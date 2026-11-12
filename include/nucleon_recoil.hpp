#ifndef BNS_NURATES_SRC_NUCLEON_RECOIL_HPP_
#define BNS_NURATES_SRC_NUCLEON_RECOIL_HPP_

//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  nucleon_recoil.hpp
//  \brief energy exchange in neutrino-nucleon scattering (nucleon recoil and
//         thermal motion) as grey M1 source coefficients
//
// The grey scattering opacity kappa_s treats neutrino-nucleon scattering as
// isoenergetic. This file adds the energy that scattering exchanges with the
// nucleons themselves: the neutrino loses energy to nucleon recoil and gains
// energy from the nucleons' thermal motion.
//
// Kernel (free, non-degenerate nucleons of mass m at temperature T): the ideal-gas
// dynamic structure factor
//
//   S(q0, q) = sqrt(m / (2 pi T q^2)) exp[-m (q0 - q^2/2m)^2 / (2 T q^2)],
//
// with q0 = w - w' and q^2 = w^2 + w'^2 - 2 w w' c. S(-q0) = exp(-q0/T) S(q0), so
// the exchange vanishes for a Fermi-Dirac spectrum at the gas temperature, and the
// kernel is exact to all orders in w/m. The rate of change of the fluid-frame
// energy density is
//
//   dJ/dt = (4 pi/(hc)^3) int dw w^2 f(w) int dc 4 pi (L0(w) + L1(w) c)/2
//           int dw' w'^2 [1 - f(w')] (w' - w) S(q0, q),
//
// where L0, L1 are the library's own Legendre coefficients of the isoenergetic
// kernel (as IsoScattLegCoeff). In the elastic limit (S -> delta(q0)) the
// transport rate of this kernel is 4 pi w^2 IsoScattTotal(w), i.e. exactly the
// spectral opacity behind kappa_s, so weak magnetism and the nucleon-degeneracy
// factor carry over. L0 and L1 are taken at the pair midpoint (w + w')/2, not at
// the incoming energy: the weak-magnetism factor depends on energy, and only a
// factor symmetric in w <-> w' keeps the detailed balance of S (taken at w, the
// imbalance at a Fermi-Dirac spectrum is ~1% of the gross exchange at T = 14 MeV).
// Number is conserved; only energy moves.
//
// dJ/dt is split into a gain G (up-scattering, q0 < 0) and a loss L (q0 > 0), both
// non-negative. They enter the grey coefficients as an energy emissivity G and an
// energy absorption opacity kappa_rec = L / (c J).
//
// Quadrature: the kernel is a Gaussian in q0 a few MeV wide, which the library's
// two-dimensional energy rule cannot resolve (it recovers a few percent of the
// exchange while still balancing at equilibrium). The rule used here is fitted to
// it: Gauss-Laguerre (6) in w scaled to the spectrum (J / 3n), Gauss-Legendre (8)
// in c, Gauss-Hermite (8) in q0 centred on the recoil peak q^2/2m with width
// sqrt(q^2 T/m) (elastic q). It reproduces a fine-grid evaluation to 0.25% for
// thermal, hot and hard spectra. It does not pair each scattering with its
// reverse, so for a Fermi-Dirac spectrum at the gas temperature the net is small
// but not zero: it relaxes the spectrum 0.02-0.4% too cold.

#include "bns_nurates.hpp"
#include "constants.hpp"
#include "distribution.hpp"
#include "formulas_beta_iso.hpp"

/* Species selection for OpacityFlags::use_nucleon_recoil */
#define BS_RECOIL_OFF 0
#define BS_RECOIL_HEAVY 1 // nux and anux only
#define BS_RECOIL_ALL 2   // all four species

/* Ideal-gas dynamic structure factor of free nucleons [MeV^-1] */
KOKKOS_INLINE_FUNCTION
BS_REAL RecoilStructureFactor(const BS_REAL q0, const BS_REAL q2,
                              const BS_REAL temp, const BS_REAL mass)
{
    constexpr BS_REAL two    = 2;
    constexpr BS_REAL two_pi = 2 * kBS_Pi;
    const BS_REAL d          = q0 - q2 / (two * mass);
    return Kokkos::sqrt(mass / (two_pi * temp * q2)) *
           Kokkos::exp(-mass * d * d / (two * temp * q2));
}

/* Legendre coefficients L0, L1 of the isoenergetic kernel at energy w, summed
 * over protons and neutrons, as IsoScattLegCoeff, with the energy-independent
 * degeneracy factors eta_pp, eta_nn [nm^-3] passed in; [MeV nm^6 s^-1] */
KOKKOS_INLINE_FUNCTION
void RecoilLegCoeffs(const BS_REAL w, OpacityParams* opacity_pars,
                     const BS_REAL eta_pp, const BS_REAL eta_nn, BS_REAL* leg_0,
                     BS_REAL* leg_1)
{
    BS_REAL R0_p = 1., R1_p = 1., R0_n = 1., R1_n = 1.;
    if (opacity_pars->use_WM_sc)
    {
        WMScatt(w, &R0_p, &R1_p, 1);
        WMScatt(w, &R0_n, &R1_n, 2);
    }
    *leg_0 = eta_pp * kBS_Iso_c0_p * R0_p + eta_nn * kBS_Iso_c0_n * R0_n;
    *leg_1 = eta_pp * kBS_Iso_c1_p * R1_p + eta_nn * kBS_Iso_c1_n * R1_n;
}

/* NucleonRecoilExchange
 *
 * Gain and loss of fluid-frame neutrino energy density from neutrino-nucleon
 * scattering with recoil, for species idx, using the distribution in
 * grey_pars->distr_pars and the spectral scale tau = J / (3 n).
 *
 * Output: *gain, *loss in [MeV nm^-3 s^-1] (the units of the energy emissivity).
 * Both are set to zero if n, J or T is not positive.
 */
KOKKOS_INLINE_FUNCTION
void NucleonRecoilExchange(GreyOpacityParams* grey_pars, const int idx,
                           const BS_REAL n, const BS_REAL J, BS_REAL* gain,
                           BS_REAL* loss)
{
    constexpr BS_REAL zero      = 0;
    constexpr BS_REAL half      = 0.5;
    constexpr BS_REAL one       = 1;
    constexpr BS_REAL two       = 2;
    constexpr BS_REAL three     = 3;
    constexpr BS_REAL sqrt_two  = 1.4142135623730951;
    constexpr BS_REAL two_pi    = 2 * kBS_Pi; // 4 pi * (1/2) of the angular factor
    constexpr BS_REAL q2_floor  = 1e-12;      // [MeV^2]
    constexpr BS_REAL mass      = kBS_Mn;     // [MeV]
    constexpr BS_REAL fd_step   = 1e-3;       // relative step for dL/dw
    constexpr BS_REAL fd_min    = 1e-6;       // [MeV]

    // Gauss-Laguerre, weights multiplied by exp(x)
    constexpr int n_lag              = 6;
    const BS_REAL x_lag[n_lag] = {0.2228466041792607, 1.1889321016726226,
                                  2.992736326059314,  5.77514356910451,
                                  9.83746741838259,   15.982873980601703};
    const BS_REAL w_lag[n_lag] = {0.5735355074227396, 1.3692525907123008,
                                  2.2606845933826705, 3.350524582355459,
                                  4.886826800210833,  7.849015945595829};
    // Gauss-Legendre on [-1, 1]
    constexpr int n_leg              = 8;
    const BS_REAL x_leg[n_leg] = {
        -0.9602898564975362, -0.7966664774136267, -0.525532409916329,
        -0.18343464249564978, 0.18343464249564978, 0.525532409916329,
        0.7966664774136267,   0.9602898564975362};
    const BS_REAL w_leg[n_leg] = {
        0.10122853629037669, 0.22238103445337434, 0.31370664587788705,
        0.36268378337836177, 0.36268378337836177, 0.31370664587788705,
        0.22238103445337434, 0.10122853629037669};
    // Gauss-Hermite, weights multiplied by exp(x^2)
    constexpr int n_her              = 8;
    const BS_REAL x_her[n_her] = {
        -2.930637420257244, -1.981656756695843, -1.1571937124467802,
        -0.3811869902073221, 0.3811869902073221, 1.1571937124467802,
        1.981656756695843,   2.930637420257244};
    const BS_REAL w_her[n_her] = {
        1.0719301442479818, 0.866752606563381,  0.7928900483864011,
        0.7645441286517294, 0.7645441286517294, 0.7928900483864011,
        0.866752606563381,  1.0719301442479818};

    *gain = zero;
    *loss = zero;

    const BS_REAL temp = grey_pars->eos_pars.temp;
    if (! (n > zero) || ! (J > zero) || ! (temp > zero))
    {
        return;
    }

    const BS_REAL tau = J / (three * n); // [MeV]

    // nucleon degeneracy factors (Bruenn C37), energy-independent
    const BS_REAL eta_pp =
        EtaNNSc(grey_pars->eos_pars.nb, temp, grey_pars->eos_pars.yp);
    const BS_REAL eta_nn =
        EtaNNSc(grey_pars->eos_pars.nb, temp, grey_pars->eos_pars.yn);

    BS_REAL g_sum = zero, l_sum = zero;

    for (int i = 0; i < n_lag; ++i)
    {
        const BS_REAL w   = tau * x_lag[i];
        const BS_REAL f_w = TotalNuF(w, &grey_pars->distr_pars, idx);
        if (! (f_w > zero))
        {
            continue;
        }

        const BS_REAL outer = tau * w_lag[i] * w * w * f_w;

        // L0, L1 and their slopes at w; the midpoint values below are linear in
        // q0, which differs from the exact midpoint by O(L'' q0^2) ~ 1e-5
        BS_REAL leg_0w, leg_1w, leg_0h, leg_1h;
        const BS_REAL dw = Kokkos::fmax(fd_step * w, fd_min); // [MeV]
        RecoilLegCoeffs(w, &grey_pars->opacity_pars, eta_pp, eta_nn, &leg_0w,
                        &leg_1w);
        RecoilLegCoeffs(w + dw, &grey_pars->opacity_pars, eta_pp, eta_nn,
                        &leg_0h, &leg_1h);
        const BS_REAL dleg_0 = (leg_0h - leg_0w) / dw;
        const BS_REAL dleg_1 = (leg_1h - leg_1w) / dw;

        for (int k = 0; k < n_leg; ++k)
        {
            const BS_REAL c     = x_leg[k];
            const BS_REAL q2_el = Kokkos::fmax(two * w * w * (one - c), q2_floor);
            const BS_REAL cen   = q2_el / (two * mass);
            const BS_REAL sig   = Kokkos::sqrt(q2_el * temp / mass);
            const BS_REAL jac   = sqrt_two * sig;

            for (int h = 0; h < n_her; ++h)
            {
                const BS_REAL q0 = cen + jac * x_her[h];
                const BS_REAL wp = w - q0;
                if (! (wp > zero))
                {
                    continue;
                }
                const BS_REAL q2 =
                    Kokkos::fmax(w * w + wp * wp - two * w * wp * c, q2_floor);
                const BS_REAL f_wp =
                    TotalNuF(wp, &grey_pars->distr_pars, idx);
                // angular factor (L0 + L1 c) at the pair midpoint w - q0/2
                const BS_REAL leg_0 = leg_0w - half * q0 * dleg_0;
                const BS_REAL leg_1 = leg_1w - half * q0 * dleg_1;
                const BS_REAL ang   = w_leg[k] * (leg_0 + leg_1 * c);
                const BS_REAL r   = outer * ang * w_her[h] * jac * wp * wp *
                                  (one - f_wp) *
                                  RecoilStructureFactor(q0, q2, temp, mass) * q0;
                if (q0 < zero)
                {
                    g_sum -= r; // up-scattering: the neutrino gains energy
                }
                else
                {
                    l_sum += r; // down-scattering: the neutrino loses energy
                }
            }
        }
    }

    *gain = kBS_FourPi_hc3 * two_pi * g_sum;
    *loss = kBS_FourPi_hc3 * two_pi * l_sum;

    return;
}

/* Whether species idx gets the recoil term for a given use_nucleon_recoil value */
KOKKOS_INLINE_FUNCTION
bool NucleonRecoilSelected(const int use_nucleon_recoil, const int idx)
{
    if (use_nucleon_recoil == BS_RECOIL_ALL)
    {
        return true;
    }
    if (use_nucleon_recoil == BS_RECOIL_HEAVY)
    {
        return (idx == id_nux) || (idx == id_anux);
    }
    return false;
}

#endif // BNS_NURATES_SRC_NUCLEON_RECOIL_HPP_
