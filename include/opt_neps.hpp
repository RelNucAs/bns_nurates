//! \file  opt_neps.hpp
//  \brief NES/NPS kernels with per-cell constants hoisted and Fermi integrals
//         of all orders at one argument sharing a single exponential.
//         Same algebra and branches as kernel_nes.hpp.

#ifndef BNS_NURATES_INCLUDE_OPT_NEPS_HPP_
#define BNS_NURATES_INCLUDE_OPT_NEPS_HPP_

#include "bns_nurates.hpp"
#include "constants.hpp"
#include "functions.hpp"
#include "kernel_nes.hpp"

/* Per-cell constants of one lepton (electrons for NES, positrons for NPS) */
struct NEPSLeptonTerms
{
    BS_REAL eta;
    BS_REAL f3_e, f4_e, f5_e;      // FDI_p{3,4,5}(eta)
    BS_REAL fd0, fdi0, fdi1, fdi2; // FermiDistr(0,1,eta), FDI_0, FDI_p1, FDI_p2 at eta
};

struct NEPSEOSTerms
{
    BS_REAL T, pref; // pref = kBS_NEPS_Const * T^2
    NEPSLeptonTerms lep[2]; // 0: electrons (NES), 1: positrons (NPS)
};

KOKKOS_INLINE_FUNCTION
NEPSEOSTerms NEPSPrecomputeEOS(const MyEOSParams* eos)
{
    NEPSEOSTerms c;
    c.T    = eos->temp;
    c.pref = kBS_NEPS_Const * POW2(c.T);
    for (int l = 0; l < 2; ++l)
    {
        NEPSLeptonTerms* L = &c.lep[l];
        L->eta = (l == 0) ? eos->mu_e / c.T : -eos->mu_e / c.T;
        BS_REAL f[3];
        FDI_p3to5(L->eta, f);
        L->f3_e = f[0]; L->f4_e = f[1]; L->f5_e = f[2];
        L->fd0  = FermiDistr(0., 1., L->eta);
        L->fdi0 = FDI_0(L->eta);
        L->fdi1 = FDI_p1(L->eta);
        L->fdi2 = FDI_p2(L->eta);
    }
    return c;
}

/* Kernel base (before detailed-balance factors) of NES (l=0) or NPS (l=1) for
 * the four species, with the same b-coefficient assignment as NESKernels /
 * NPSKernels. */
KOKKOS_INLINE_FUNCTION
void NEPSBaseOneLepton(const BS_REAL w, const BS_REAL wp, const NEPSLeptonTerms* L,
                       const int l, const BS_REAL pref, BS_REAL* base)
{
    constexpr BS_REAL zero = 0, one = 1, six = 6;
    const BS_REAL x   = (w > wp) ? w : wp;
    const BS_REAL y   = (w < wp) ? w : wp;
    const BS_REAL eta = L->eta;

    // (b1, b2) per species; NPS swaps the two
    const BS_REAL bA[total_num_species] = {kBS_NEPS_BPlus, kBS_NEPS_BZero, kBS_NEPS_BMinus, kBS_NEPS_BZero};
    const BS_REAL bB[total_num_species] = {kBS_NEPS_BZero, kBS_NEPS_BPlus, kBS_NEPS_BZero, kBS_NEPS_BMinus};

    if (y > eta * kTaylorSeriesEpsilon)
    {
        const int sign = 2 * signbit(wp - w) - 1;
        BS_REAL fwp[5], fw[5], fabs_[3];
        FDI_p1to5(eta - wp, fwp);
        FDI_p1to5(eta - w, fw);
        FDI_p3to5(eta - Kokkos::fabs(w - wp), fabs_);
        BS_REAL fdi_diff_w[5], fdi_diff_abs[3];
        for (int k = 0; k < 5; ++k) fdi_diff_w[k] = fwp[k] - fw[k];
        fdi_diff_abs[0] = L->f3_e - fabs_[0];
        fdi_diff_abs[1] = L->f4_e - fabs_[1];
        fdi_diff_abs[2] = L->f5_e - fabs_[2];
        for (int s = 0; s < total_num_species; ++s)
        {
            const BS_REAL b1 = (l == 0) ? bA[s] : bB[s];
            const BS_REAL b2 = (l == 0) ? bB[s] : bA[s];
            base[s] = pref * MezzacappaIntOut(w, wp, x, y, sign, b1, b2, fdi_diff_w, fdi_diff_abs);
        }
    }
    else if (x > eta * kTaylorSeriesEpsilon)
    {
        const int sign = 2 * signbit(wp - w) - 1;
        const BS_REAL fdis[5] = {
            FDI_p2(eta - x) - L->fdi2, FDI_p1(eta - x), L->fdi1,
            FDI_0(eta - x) + y * FermiDistr(zero, one, eta - x) / six,
            L->fdi0 - y * L->fd0 / six};
        for (int s = 0; s < total_num_species; ++s)
        {
            const BS_REAL b1 = (l == 0) ? bA[s] : bB[s];
            const BS_REAL b2 = (l == 0) ? bB[s] : bA[s];
            base[s] = pref * MezzacappaIntOneEnergy(x, y, sign, b1, b2, fdis);
        }
    }
    else
    {
        const BS_REAL fdis[3] = {L->fd0, L->fdi0, L->fdi1};
        for (int s = 0; s < total_num_species; ++s)
        {
            const BS_REAL b1 = (l == 0) ? bA[s] : bB[s];
            const BS_REAL b2 = (l == 0) ? bB[s] : bA[s];
            base[s] = pref * MezzacappaIntTwoEnergies(w, wp, x, y, b1, b2, fdis);
        }
    }
}

/* Positron scattering (NPS) is suppressed as exp(-eta_e): above this electron
 * degeneracy it can be skipped. On a BNS slice, a cut at 10 changes the NEPS
 * coefficients by < 1e-6. Disabled by default. */
constexpr BS_REAL kBS_NPS_EtaCut = 1e300;

KOKKOS_INLINE_FUNCTION
bool NEPSUseNPS(const NEPSEOSTerms* c)
{
    return !(c->lep[0].eta > kBS_NPS_EtaCut);
}

/* Same as InelasticScattKernels for omega = nu, omega_prime = nu_p [MeV] */
KOKKOS_INLINE_FUNCTION
MyKernelOutput InelasticScattKernelsFromEOS(const BS_REAL nu, const BS_REAL nu_p,
                                            const NEPSEOSTerms* c, const bool use_nps)
{
    constexpr BS_REAL zero = 0, one = 1;
    const BS_REAL w  = nu / c->T;
    const BS_REAL wp = nu_p / c->T;

    // NEPSExpFunc(wp - w) and NEPSExpFunc(w - wp) share exp(-|w - wp|)
    const BS_REAL d       = wp - w;
    const BS_REAL e_      = SafeExp(-Kokkos::fabs(d));
    const BS_REAL denom_d = one - e_ + (d == zero);
    const int neg_d       = signbit(d);
    const int neg_md      = signbit(w - wp);
    const BS_REAL exp_factor           = (neg_d - (!neg_d) * e_) / denom_d;
    const BS_REAL exp_factor_exchanged = (neg_md - (!neg_md) * e_) / denom_d;

    BS_REAL nes[total_num_species], nps[total_num_species] = {0, 0, 0, 0};
    NEPSBaseOneLepton(w, wp, &c->lep[0], 0, c->pref, nes);
    if (use_nps) NEPSBaseOneLepton(w, wp, &c->lep[1], 1, c->pref, nps);

    MyKernelOutput out;
    for (int s = 0; s < total_num_species; ++s)
    {
        out.abs[s] = nes[s] * exp_factor + nps[s] * exp_factor;
        out.em[s]  = -nes[s] * exp_factor_exchanged + -nps[s] * exp_factor_exchanged;
    }
    return out;
}

#endif // BNS_NURATES_INCLUDE_OPT_NEPS_HPP_
