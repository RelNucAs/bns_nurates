// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file opacities.h
//  \brief header files for opacity functions

#ifndef BNS_NURATES_SRC_OPACITIES_OPACITIES_H_
#define BNS_NURATES_SRC_OPACITIES_OPACITIES_H_

#include "bns_nurates.hpp"
#include "constants.hpp"
#include "functions.hpp"
#include "weak_magnetism.hpp"

#define kirchoff_flag false

/*===========================================================================*/

// nu_abs_em_beta.c
// \brief Computes emissivity and absorptivity for neutrino absorption on
// neutron
//        and for antineutrino absorption on proton and their inverse reactions
//
// Reference: Bruenn, 1985
// (https://articles.adsabs.harvard.edu/pdf/1985ApJS...58..771B)
//
// Also include:
//           Possible inclusion of phase space, recoil, weak magnetism
//           correction as in Horowitz, 2002
//           (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001)
//
//           Possible inclusion of nucleon interation correction as in
//           Hempel, 2015
//           (https://journals.aps.org/prc/abstract/10.1103/PhysRevC.91.055807)

// Constants for beta reactions
#define mu_thres 1.0e-02 // mu_hat threshold value in EtaNNAbs evaluation

// ----------------------------------------------------------------------------------

/* Inputs:
 * 	omega    [MeV] : (anti)neutrino energy
 * 	nb      [nm-3] : baryon number density
 * 	T        [MeV] : temperature
 * 	lep_mass [MeV] : lepton mass (electron/muon)
 * 	yp             : proton fraction
 * 	yn             : neutron fraction
 * 	mu_l     [MeV] : lepton chemical potential, rest mass included
 * 	mu_hat   [MeV] : neutron-proton chemical potential difference, rest mass NOT
 * included dU       [MeV] : nuclear interaction correction on mu_hat
 */

/* Generic function for nucleon phase space integration
 * Computes absorptivity from Eqn. (C14) of Bruenn
 *
 * Note:
 *
 * mu_hat [MeV] is neutron-proton chemical potential difference, rest mass NOT
 * included This function is never called directly.
 *
 *  Output: j [nm^-3]
 */
KOKKOS_INLINE_FUNCTION
bs_real EtaNNAbs(const bs_real n_in, const bs_real n_out, const bs_real mu_hat,
                 const bs_real T)
{
    // if no nucleons, enforce zero rates
    // if (n_in == 0.)
    // {
    // return 0.;
    // }

    // if mu_hat too small, neglect nucleon degeneracy as backup
    if (fabs(mu_hat) < mu_thres)
    {
        return n_in;
    }

    // Eq.(C14), [nm-3]
    const bs_real etanp = (n_out - n_in) / (SafeExp(-mu_hat / T) - 1.);

    BS_ASSERT(isfinite(etanp),
              "etanp factor in beta-process kernel is not finite.")

    // backup if etanp is negative
    if (etanp < 0.)
    {
        return n_in;
    }

    return etanp;
}

/* Perform phase space integration for:
 *    X + n -> X + p
 * Inputs:
 *      nn [MeV], np [MeV], mu_hat [MeV], T [MeV]
 * Output:
 *      Eta_np from Bruenn (C14) [nm^-3]
 */
KOKKOS_INLINE_FUNCTION
bs_real EtaNP(const bs_real nn, const bs_real np, const bs_real mu_hat,
              const bs_real T)
{
    return EtaNNAbs(nn, np, mu_hat, T);
}

/* Perform phase space integration for:
 *    X + p -> X + n
 * Inputs:
 *      nn [MeV], np [MeV], mu_hat [MeV], T [MeV]
 * Output:
 *      Eta_pn [nm^-3]
 */
KOKKOS_INLINE_FUNCTION
bs_real EtaPN(const bs_real nn, const bs_real np, const bs_real mu_hat,
              const bs_real T)
{
    return EtaNNAbs(np, nn, -mu_hat, T);
}

// @TODO: add effective mass correction to the rates

/* Compute opacities for:
 *
 * 1. Neutrino absorption on neutron (nu + n -> l- + p)
 * 2. Antineutrino absorption on neutron (anul + p -> l+ + n)
 *
 * Inputs:
 *    omega [MeV], mLep [MeV], muLep [MeV]
 * Outputs: j_x and 1/lambda_x for electron neutrino and electron antineutrino
 *    out[0]: j_nue [s^-1], out[1]: 1/lamda_nue [s^-1], out[2]: j_anue [s^-1],
 * out[3]: 1/lamda_anue [s^-1]
 */
KOKKOS_INLINE_FUNCTION
void AbsOpacitySingleLep(const bs_real omega, OpacityParams* opacity_pars,
                         MyEOSParams* eos_pars, const bs_real mass_lepton,
                         const bs_real mu_lepton, bs_real* out)
{
    const bs_real mass_lepton_squared = POW2(mass_lepton);

    bs_real Qprime, mu_np;
    bs_real etanp, etapn;
    bs_real E_e, E_p;
    bs_real E_e_squared, E_p_squared;
    bs_real cap_term = 0., dec_term = 0.;
    bs_real fd_e, fd_p;
    bs_real R = 1., Rbar = 1.;
    bs_real dU = 0., dQ = kBS_Q;

    const bs_real nb   = eos_pars->nb;   // Number baryon density [nm-3]
    const bs_real T    = eos_pars->temp; // Temperature [MeV]
    const bs_real yp   = eos_pars->yp;   // Proton fraction
    const bs_real yn   = eos_pars->yn;   // Neutron fraction
    const bs_real mu_p = eos_pars->mu_p; // Proton chemical potential [MeV]
    const bs_real mu_n = eos_pars->mu_n; // Neutron chemical potential [MeV]

    const bs_real nn = nb * yn; // Neutron number density [nm-3]
    const bs_real np = nb * yp; // Proton number density  [nm-3]

    // Mean field corrections
    if (opacity_pars->use_dU)
        dU = eos_pars->dU; // [MeV]
    if (opacity_pars->use_dm_eff)
        dQ = eos_pars->dm_eff; // [MeV]

    // Neutron minus proton chem. potentials (corrected for the mass difference)
    const bs_real mu_hat = mu_n - mu_p - dQ; // [MeV]

    Qprime = dQ + dU;     // [MeV], Eq.(79) in Hempel
    mu_np  = mu_hat - dU; // [MeV], Eq.(80,86) in Hempel

    etanp = EtaNP(nn, np, mu_np, T); // Eq. (C14)
    etapn = EtaPN(nn, np, mu_np, T);

    // Phase space, recoil and weak magnetism correction
    if (opacity_pars->use_WM_ab)
        WMAbsEm(omega, &R, &Rbar);

    // Electron/muon-type neutrino
    E_e = omega + Qprime; // Electron energy [MeV]
    E_p = -E_e;           // Positron energy [MeV]

    E_e_squared = E_e * E_e;
    E_p_squared = E_p * E_p;

    fd_e = FermiDistr(E_e, T, +mu_lepton);
    fd_p = FermiDistr(E_p, T, -mu_lepton);

    // @TODO: HeavisideTanhApprox is 0.5 instead of 1 in E = m
    if (E_e - mass_lepton >= 0.)
    {
        cap_term = E_e_squared * sqrt(1. - mass_lepton_squared / E_e_squared) *
                   R; // * HeavisideTanhApprox(E_e - mass_lepton)
    }

    if (opacity_pars->use_decay)
    {
        if (E_p - mass_lepton >= 0.)
        {
            dec_term = E_p_squared *
                       sqrt(1. - mass_lepton_squared /
                                     E_p_squared); // * HeavisideTanhApprox(E_p
                                                   // - mass_lepton)
        }
    }

    // @TODO: eventually think about a specifically designed function for
    // (1-FermiDistr)

    // Neutrino emissivity [s^-1], Eq.(C15), remove c to get output in nm^-1
    out[1] =
        kBS_Beta_Const * etapn * (cap_term * fd_e + dec_term * (1. - fd_p));
    // Neutrino absorptivity [s^-1]
    out[0] = out[1] * SafeExp((omega - (mu_p + mu_lepton - mu_n)) / T);

    // without detailed balance
    // out[0] = kAbsEmConst * etanp * (cap_term * (1. - fd_e) + dec_term *
    // fd_p); // Neutrino absorptivity [s-1], Eq.(C13) bs_real mu_nue =
    // (eos_pars->mu_e - eos_pars->mu_n + eos_pars->mu_p) / temp; out[0] =
    // kAbsEmConst * etanp * cap_term / (1. + exp(eos_pars->mu_e / temp -
    // FDI_p5(mu_nue)/FDI_p4(mu_nue)));

    cap_term = 0.;
    dec_term = 0.;

    E_p = omega - Qprime; // Positron energy [MeV]
    E_e = -E_p;           // Electron energy [MeV]

    E_e_squared = E_e * E_e;
    E_p_squared = E_p * E_p;

    fd_e = FermiDistr(E_e, T, +mu_lepton);
    fd_p = FermiDistr(E_p, T, -mu_lepton);

    if (E_p - mass_lepton >= 0.)
    {
        cap_term = E_p_squared * sqrt(1. - mass_lepton_squared / E_p_squared) *
                   Rbar; // * HeavisideTanhApprox(E_p - mass_lepton)
    }

    if (opacity_pars->use_decay)
    {
        if (E_e - mass_lepton >= 0.)
        {
            dec_term = E_e_squared *
                       sqrt(1. - mass_lepton_squared /
                                     E_e_squared); // * HeavisideTanhApprox(E_e
                                                   // - mass_lepton)
        }
    }

    // Antineutrino emissivity [s^-1], Eq.(C20), remove c to get output in nm^-1
    out[3] =
        kBS_Beta_Const * etanp * (cap_term * fd_p + dec_term * (1. - fd_e));
    // Antineutrino absorptivity [s^-1]
    out[2] = out[3] * SafeExp((omega - (mu_n - mu_p - mu_lepton)) / T);

    // without detailed balance
    // out[2] = kAbsEmConst * etapn * (cap_term * (1 - fd_p) + dec_term * fd_e);
    // // Antineutrino absorptivity [s-1], Eq.(C19)

    BS_ASSERT(isfinite(out[1]) && out[1] >= 0.,
              "Invalid beta-process nue emissivity.");
    BS_ASSERT(isfinite(out[3]) && out[3] >= 0.,
              "Invalid beta-process anue emissivity.");
    BS_ASSERT(isfinite(out[0]) && out[0] >= 0.,
              "Invalid beta-process nue absorptivity.");
    BS_ASSERT(isfinite(out[2]) && out[2] >= 0.,
              "Invalid beta-process anue absorptivity.");

    return;
}

/* Compute the absortivity and inverse mean free path
 *
 * Outputs:
 *      em_nue = j_nu [s^-1] for nue
 *      ab_nue = 1/lambda_nu [s^-1] for nue
 *      em_anue = j_nu [s^-1] for anue
 *      ab_anue = 1/lambda_nu [s^-1] for anue
 *
 * @TODO: add support for muons
 */
KOKKOS_INLINE_FUNCTION
MyOpacity AbsOpacity(const bs_real omega, OpacityParams* opacity_pars,
                     MyEOSParams* eos_pars)
{
    MyOpacity MyOut = {0}; // initialize to zero

    // Electron (anti)neutrino
    bs_real el_out[4] = {0.0};
    AbsOpacitySingleLep(omega, opacity_pars, eos_pars, kBS_Me, eos_pars->mu_e,
                        el_out);

    MyOut.abs[id_nue]  = el_out[0];
    MyOut.em[id_nue]   = el_out[1];
    MyOut.abs[id_anue] = el_out[2];
    MyOut.em[id_anue]  = el_out[3];

    // Uncomment the following when considering also muons
    // // Muon (anti)neutrino
    // bs_real mu_out[4] = {0.0};
    // AbsOpacitySingleLep(omega, opacity_pars, eos_pars, kBS_Mmu,
    // eos_pars->mu_mu, mu_out);

    // MyOut.abs[id_num] = mu_out[0];
    // MyOut.em[id_num] = mu_out[1];
    // MyOut.abs[id_anum] = mu_out[2];
    // MyOut.em[id_anum] = mu_out[3];

    return MyOut;
}

/* Compute the stimulated absorption opacities
 *
 * For all species:
 *      em_nu = j_nu [s^-1]
 *      ab_nu = j_nu + 1/lambda_nu [s^-1]
 *
 * Both opacities are energy dependent
 * @TODO: add support for muons
 */
KOKKOS_INLINE_FUNCTION
MyOpacity StimAbsOpacity(const bs_real omega, OpacityParams* opacity_pars,
                         MyEOSParams* eos_pars)
{
    MyOpacity abs_opacity = AbsOpacity(omega, opacity_pars, eos_pars);

    abs_opacity.abs[id_nue] = abs_opacity.abs[id_nue] + abs_opacity.em[id_nue];
    abs_opacity.abs[id_anue] =
        abs_opacity.abs[id_anue] + abs_opacity.em[id_anue];

    return abs_opacity;
}

KOKKOS_INLINE_FUNCTION
void BetaOpacitiesTable(MyQuadrature* quad, MyEOSParams* eos_pars,
                        OpacityParams* opacity_pars, bs_real t,
                        M1MatrixKokkos2D* out)
{
    const int n = quad->nx;

    MyOpacity beta_1, beta_2;

    for (int i = 0; i < n; ++i)
    {

        beta_1 = StimAbsOpacity(t * quad->points[i], opacity_pars, eos_pars);
        beta_2 = StimAbsOpacity(t / quad->points[i], opacity_pars, eos_pars);

        for (int idx = 0; idx < total_num_species; ++idx)
        {
            out->m1_mat_ab[idx][0][i] = beta_1.abs[idx];
            out->m1_mat_em[idx][0][i] = beta_1.em[idx];

            out->m1_mat_ab[idx][0][n + i] = beta_2.abs[idx];
            out->m1_mat_em[idx][0][n + i] = beta_2.em[idx];
        }
    }

    return;
}

/*
// (Anti)neutrino absorption on nucleons
MyOpacity AbsOpacity(const bs_real omega, OpacityParams* opacity_pars,
                     MyEOSParams* eos_pars);

// Stimulated absoption version
MyOpacity StimAbsOpacity(const bs_real omega, OpacityParams* opacity_pars,
                         MyEOSParams* eos_pars);

// Build matrix for integration
void BetaOpacitiesTable(MyQuadrature* quad, MyEOSParams* eos_pars,
                        OpacityParams* opacity_pars, bs_real t, M1Matrix* out);

*/
/*===========================================================================*/

// nu_scatt_iso.c

/*
 * @brief Implementation of iso-energetic scattering of neutrinos on nucleons.
 *
 * This file implements iso-energetic scattering on nucleons from [Bruenn
 * 1985](https://articles.adsabs.harvard.edu/pdf/1985ApJS...58..771B) with
 * optional inclusion of phase space, recoil and weak magnetism corrections from
 * [Horowitz
 * 2002](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001)
 *
 * @note The iso-energetic scattering kernel is same for production and
 * absorption and for all neutrino species.
 */


/**
 * @fn bs_real EtaNNSc(const bs_real nb, const bs_real temp, const bs_real yN)
 * @brief Computes degeneracy parameter \f$\eta_{NN}\f$ from Eqn. (C37) of
 * Bruenn.
 *
 * @param nb    baryon number density \f$[nm^{-3}]\f$
 * @param temp  temperature \f$[MeV]\f$
 * @param yN    proton/neutron fraction
 * @return      The degeneracy parameter \f$\eta_{NN} [nm^-3]\f$
 */
KOKKOS_INLINE_FUNCTION
bs_real EtaNNSc(const bs_real nb, const bs_real temp, const bs_real yN)
{
    // nucleon (neutron/proton) number density
    const bs_real nN = yN * nb; // [nm^-3]

    // Enforce zero rates if no nucleons are present
    if (nN <= 0.)
    {
        return 0.;
    }

    // Linear interpolation between degenerate and non-degnerate limit in
    // Eq.(C37)
    // Fermi energy computation
    const bs_real eFN = kBS_Iso_eF * pow(nN, 2. / 3.); // [MeV]
    const bs_real aux = 1.5 * temp / eFN;

    return nN * aux / sqrt(1. + POW2(aux)); // [nm-3]
}

/**
 * @fn bs_real IsoScattNucleon(const bs_real omega, OpacityParams *opacity_pars,
 * MyEOSParams *eos_pars, const bs_real yN, const int reacflag)
 * @brief Computes Spectral scattering opacity for iso-energetic scattering on
 * nucleons (protons/neutrons)
 * @param omega         neutrino energy \f$[MeV]\f$
 * @param opacity_pars  structure containing opacity parameters
 * @param eos_pars      structure containing equation of state parameters (needs
 * baryon number density \f$[nm^{-3}]\f$ and temperature \f$[MeV]\f$)
 * @param yN            proton/neutron fraction
 * @param reacflag      choice of nucleon (1: proton scattering 2: neutron
 * scattering)
 * @return              "Eq.(A41)" \f$[MeV nm^{3} s^{-1}]\f$
 */
KOKKOS_INLINE_FUNCTION
bs_real IsoScattNucleon(const bs_real omega, OpacityParams* opacity_pars,
                        MyEOSParams* eos_pars, const bs_real yN,
                        const int reacflag)
{
    bs_real R0 = 1., R1 = 1.;
    bs_real leg_0, leg_1;

    // Number baryon density
    const bs_real nb = eos_pars->nb; // [nm^-3]
    // Temperature
    const bs_real temp = eos_pars->temp; // [MeV]

    // Degeneracy parameter eta_NN [nm^-3] from Eqn. (C37) of Bruen
    const bs_real etaNN = EtaNNSc(nb, temp, yN);

    // Phase space, recoil and weak magnetism corrections
    // R0 (R1) is the correction to the zeroth (first) Legendre coefficient
    if (opacity_pars->use_WM_sc)
    {
        WMScatt(omega, &R0, &R1, reacflag);
    }

    if (reacflag == 1)
    {
        // Scattering on proton
        leg_0 = kBS_Iso_c0_p * R0;
        leg_1 = kBS_Iso_c1_p * R1; // [MeV nm^6 s^-1]
    }
    else if (reacflag == 2)
    {
        // Scattering on neutron
        leg_0 = kBS_Iso_c0_n * R0;
        leg_1 = kBS_Iso_c1_n * R1; // [MeV nm^6 s^-1]
    }

    // Eq.(A41)
    return etaNN * (leg_0 - leg_1 / 3.); // [MeV nm^3 s-1]
}

/**
 * @fn bs_real IsoScattProton(bs_real omega, OpacityParams *opacity_pars,
 * MyEOSParams *eos_pars)
 * @brief Computes the spectral scattering opacity for scattering of neutrinos
 * on protons
 * @param omega         neutrino energy \f$[MeV]\f$
 * @param opacity_pars  structure for opacity parameters
 * @param eos_pars      structure for equation of state parameters
 * @return              "Eq.(A41)" \f$[MeV nm{^3} s^{-1}]\f$
 */
KOKKOS_INLINE_FUNCTION
bs_real IsoScattProton(const bs_real omega, OpacityParams* opacity_pars,
                       MyEOSParams* eos_pars)
{
    return IsoScattNucleon(omega, opacity_pars, eos_pars, eos_pars->yp, 1);
}

/**
 * @fn bs_real IsoScattNeutron(bs_real omega, OpacityParams *opacity_pars,
 * MyEOSParams *eos_pars)
 * @brief Computes the spectral scattering opacity for scattering of neutrinos
 * on neutrons
 * @param omega         neutrino energy \f$[MeV]\f$
 * @param opacity_pars  structure for opacity parameters
 * @param eos_pars      structure for equation of state parameters
 * @return              "Eq.(A41)" \f$[MeV nm{^3} s^{-1}]\f$
 */
KOKKOS_INLINE_FUNCTION
bs_real IsoScattNeutron(const bs_real omega, OpacityParams* opacity_pars,
                        MyEOSParams* eos_pars)
{
    return IsoScattNucleon(omega, opacity_pars, eos_pars, eos_pars->yn, 2);
}

/**
 * @fn bs_real IsoScattTotal(bs_real omega, OpacityParams *opacity_pars,
 * MyEOSParams *eos_pars)
 * @brief Computes the total spectral scattering opacity for scattering of
 * neutrinos on protons and neutrons
 * @param omega         neutrino energy \f$[MeV]\f$
 * @param opacity_pars  structure for opacity parameters
 * @param eos_pars      structure for equation of state parameters
 * @return              "Eq.(A41)" \f$[MeV nm{^3} s^{-1}]\f$
 */
KOKKOS_INLINE_FUNCTION
bs_real IsoScattTotal(const bs_real omega, OpacityParams* opacity_pars,
                      MyEOSParams* eos_pars)
{
    // proton contribution
    const bs_real iso_nu_p = IsoScattProton(omega, opacity_pars, eos_pars);
    // neutron contribution
    const bs_real iso_nu_n = IsoScattNeutron(omega, opacity_pars, eos_pars);

    return iso_nu_p + iso_nu_n;
}


/**
 * @fn bs_real IsoScattLegCoeff(const bs_real omega, OpacityParams
 * *opacity_pars, MyEOSParams *eos_pars, const int l)
 * @brief Computes Spectral scattering opacity for iso-energetic scattering on
 * nucleons (protons/neutrons)
 * @param omega         neutrino energy \f$[MeV]\f$
 * @param opacity_pars  structure containing opacity parameters
 * @param eos_pars      structure containing equation of state parameters (needs
 * baryon number density \f$[nm^{-3}]\f$ and temperature \f$[MeV]\f$)
 * @param l             order of Legendre coefficient
 * @return              Legendre coefficient of order l of the isoscattering
 * kernel
 */
KOKKOS_INLINE_FUNCTION
bs_real IsoScattLegCoeff(const bs_real omega, OpacityParams* opacity_pars,
                         MyEOSParams* eos_pars, const int l)
{
    BS_ASSERT(l >= 0 && l <= 1); // Legendre order must be either zero or one

    bs_real R0_n = 1., R1_n = 1.;
    bs_real R0_p = 1., R1_p = 1.;
    bs_real leg;

    const bs_real nb   = eos_pars->nb;   // Number baryon density [nm^-3]
    const bs_real temp = eos_pars->temp; // Temperature [MeV]
    const bs_real yp   = eos_pars->yp;   // Proton fraction
    const bs_real yn   = eos_pars->yn;   // Neutron fraction

    // degeneracy parameter from Eqn. (C37) of Bruenn
    const bs_real eta_pp =
        EtaNNSc(nb, temp, yp); // scattering on protons [nm^-3]
    const bs_real eta_nn =
        EtaNNSc(nb, temp, yn); // scattering on neutrons [nm^-3]

    // Phase space, recoil and weak magnetism corrections
    // R0 (R1) is the correction to the zeroth (first) Legendre coefficient
    if (opacity_pars->use_WM_sc)
    {
        WMScatt(omega, &R0_p, &R1_p, 1);
        WMScatt(omega, &R0_n, &R1_n, 2);
    }

    if (l == 0)
    {
        leg = eta_pp * kBS_Iso_c0_p * R0_p +
              eta_nn * kBS_Iso_c0_n * R0_n; // [MeV nm^6 s^-1]
    }
    else
    {
        leg = eta_pp * kBS_Iso_c1_p * R1_p +
              eta_nn * kBS_Iso_c1_n * R1_n; // [MeV nm^6 s^-1]
    }

    return leg;
}

/*===========================================================================*/

// pair opacities
#ifdef GSL_INCLUDES_H_
MyKernelQuantity PairEmissivityAbsorptivityIntegrandFermi(
    bs_real var, MyEOSParams* my_eos_params, MyKernelParams* my_kernel_params);
MyKernelQuantity PairOpacitiesFermi(MyQuadrature* quad,
                                    MyEOSParams* my_eos_params,
                                    MyKernelParams* my_kernel_params);
#endif // GSL_INCLUDES_H

// bremsstrahlung opacities
#ifdef GSL_INCLUDES_H_
MyKernelQuantity
BremEmissivityAbsorptivityIntegrandFermi(bs_real omega_prime,
                                         MyEOSParams* my_eos_params,
                                         MyKernelParams* my_kernel_params);
MyKernelQuantity BremOpacitiesFermi(MyQuadrature* quad,
                                    MyEOSParams* my_eos_params,
                                    MyKernelParams* my_kernel_params);
#endif // GSL_INCLUDES_H

#endif // BNS_NURATES_SRC_OPACITIES_OPACITIES_H_
