//! \file  opt_kernels.hpp
//  \brief Per-node-cached variant of the pair kernel used by the fused 2D
//         pair/brem integral (same formulas, fewer transcendental evaluations).

#ifndef BNS_NURATES_INCLUDE_OPT_KERNELS_HPP_
#define BNS_NURATES_INCLUDE_OPT_KERNELS_HPP_

#include "bns_nurates.hpp"
#include "constants.hpp"
#include "functions.hpp"
#include "kernel_pair.hpp"

/* ------------------------------------------------------------------------ */
/* Pair process: FDIs at eta -/+ y depend on one energy only                */
/* ------------------------------------------------------------------------ */

/* Per-node Fermi integrals for the pair kernel. Only F_k(-|a|) is stored
 * (a = eta -/+ y); F_k(a) is rebuilt as in FDI_pk (reflection + polynomial for
 * a > 0), and the reflected value (-1)^k F_k(-a) used when eta > y + z is the
 * stored one. Half the storage of keeping both sets. */
struct PairNodeFDI
{
    BS_REAL am, ap;   // eta - y, eta + y
    BS_REAL nm[5];    // FDI_p{1..5}(-|eta - y|)
    BS_REAL np[5];    // FDI_p{1..5}(-|eta + y|)
};

KOKKOS_INLINE_FUNCTION
void PairNodeFDIFill(const BS_REAL eta, const BS_REAL y, PairNodeFDI* out)
{
    out->am = eta - y;
    out->ap = eta + y;
    BS_REAL f[5];
    FDI_p1to5(-Kokkos::fabs(out->am), f);
    for (int k = 0; k < 5; ++k) out->nm[k] = f[k];
    FDI_p1to5(-Kokkos::fabs(out->ap), f);
    for (int k = 0; k < 5; ++k) out->np[k] = f[k];
}

/* F_k(y) from fd = F_k(-|y|), same arithmetic as the y > 0 branch of FDI_pk */
KOKKOS_INLINE_FUNCTION
BS_REAL FDIFromNegative(const int k, const BS_REAL y, BS_REAL fd)
{
    if (!(y > zero)) return fd;
    const BS_REAL y2 = y * y;
    switch (k)
    {
    case 1: return -fd + 1.64493406684822644 + 0.5 * y * y;
    case 2: return fd + y * (3.28986813369645287 + 0.333333333333333333 * y * y);
    case 3: return -fd + 11.3643939539669510 + y2 * (4.93480220054467931 + y2 * 0.25);
    case 4: return fd + y * (45.4575758158678040 + y2 * (6.57973626739290575 + y2 * 0.2));
    default:
        return -fd + 236.532261911384425 +
               y2 * (113.643939539669510 + y2 * (8.22467033424113218 + y2 * 0.166666666666666667));
    }
}

/* Values used by PairPsiCached: direct F_k(a) or reflected (-1)^k F_k(-a) */
KOKKOS_INLINE_FUNCTION
void PairNodeValues(const PairNodeFDI* N, const bool refl, BS_REAL* m, BS_REAL* p)
{
    for (int k = 0; k < 5; ++k)
    {
        const BS_REAL sg = (k % 2 == 0) ? -1 : 1; // (-1)^(k+1) for order k+1
        m[k] = refl ? sg * N->nm[k] : FDIFromNegative(k + 1, N->am, N->nm[k]);
        p[k] = refl ? sg * N->np[k] : FDIFromNegative(k + 1, N->ap, N->np[k]);
    }
}

/* Low-energy expansion of Psi_0(y, z) to 5th order in (y, z): the O(eta^6)
 * Fermi integrals cancel analytically, leaving derivatives of the logistic
 * function f(eta) = 1/(1+exp(-eta)). */
KOKKOS_INLINE_FUNCTION
BS_REAL PairPsi0Taylor(const BS_REAL y, const BS_REAL z, const BS_REAL p, const BS_REAL q)
{
    const BS_REAL f1 = p, f2 = p * q, f3 = p * (1 - 6 * p), f4 = p * q * (1 - 12 * p);
    const BS_REAL f5 = p * (1 - 30 * p + 120 * p * p);
    const BS_REAL y2 = y * y, z2 = z * z, yz = y * z;
    return 8 * yz * (y + z) *
           (-f1 / 9 + f2 * (y - z) / 36 - f3 * (22 * y2 + 5 * yz + 22 * z2) / 2160 +
            f4 * (y - z) * (8 * y2 + 7 * yz + 8 * z2) / 4320 -
            f5 * (55 * y2 * y2 + 35 * y2 * yz + 77 * yz * yz + 35 * yz * z2 + 55 * z2 * z2) / 151200);
}

/* Same algebra as PairPsi, with per-node FDIs passed in */
KOKKOS_INLINE_FUNCTION
void PairPsiCached(const BS_REAL y, const BS_REAL z, const BS_REAL eta,
                   const BS_REAL FDI_p3_e, const BS_REAL FDI_p4_e,
                   const PairNodeFDI* Yn, const PairNodeFDI* Zn, BS_REAL* psi_out)
{
    constexpr BS_REAL four             = 4;
    constexpr BS_REAL fifteen          = 15;
    constexpr BS_REAL twenty           = 20;
    constexpr BS_REAL forty            = 40;
    constexpr BS_REAL eighty           = 80;
    constexpr BS_REAL onehundredtwenty = 120;
    constexpr BS_REAL twohundred       = 200;

    if (eta - y - z > twohundred)
    {
        psi_out[0] = zero;
        psi_out[1] = zero;
        return;
    }
    constexpr BS_REAL kPairTaylorMax = 0.5;
    if (y + z < kPairTaylorMax)
    {
        const BS_REAL em = Kokkos::exp(-Kokkos::fabs(eta));
        const BS_REAL p  = em / POW2(1 + em);
        const BS_REAL q  = -Kokkos::tanh(eta / 2);
        psi_out[0] = PairPsi0Taylor(y, z, p, q);
        psi_out[1] = PairPsi0Taylor(z, y, p, q);
        return;
    }
    // eta > y + z: all arguments positive, use F_k(a) -> (-1)^k F_k(-a)
    // (the polynomial parts cancel identically)
    const bool refl = (eta > y + z);
    struct { BS_REAL m[5], p[5]; } Yv, Zv;
    PairNodeValues(Yn, refl, Yv.m, Yv.p);
    PairNodeValues(Zn, refl, Zv.m, Zv.p);
    const auto* Y = &Yv;
    const auto* Z = &Zv;
    const BS_REAL FDI_p1_emy = Y->m[0], FDI_p1_epy = Y->p[0];
    const BS_REAL FDI_p1_emz = Z->m[0], FDI_p1_epz = Z->p[0];
    const BS_REAL FDI_p2_emy = Y->m[1], FDI_p2_epy = Y->p[1];
    const BS_REAL FDI_p2_emz = Z->m[1], FDI_p2_epz = Z->p[1];
    const BS_REAL FDI_p3_emy = Y->m[2], FDI_p3_epy = Y->p[2];
    const BS_REAL FDI_p3_emz = Z->m[2], FDI_p3_epz = Z->p[2];
    const BS_REAL FDI_p4_emy = Y->m[3], FDI_p4_epy = Y->p[3];
    const BS_REAL FDI_p4_emz = Z->m[3], FDI_p4_epz = Z->p[3];
    const BS_REAL FDI_p5_emy = Y->m[4], FDI_p5_epy = Y->p[4];
    const BS_REAL FDI_p5_emz = Z->m[4], FDI_p5_epz = Z->p[4];
    const BS_REAL emyz = eta - y - z, epyz = eta + y + z;
    const BS_REAL FDI_p3_emyz = refl ? -FDI_p3(-emyz) : FDI_p3(emyz);
    const BS_REAL FDI_p3_epyz = refl ? -FDI_p3(-epyz) : FDI_p3(epyz);
    const BS_REAL FDI_p4_emyz = refl ?  FDI_p4(-emyz) : FDI_p4(emyz);
    const BS_REAL FDI_p4_epyz = refl ?  FDI_p4(-epyz) : FDI_p4(epyz);
    const BS_REAL FDI_p5_emyz = refl ? -FDI_p5(-emyz) : FDI_p5(emyz);
    const BS_REAL FDI_p5_epyz = refl ? -FDI_p5(-epyz) : FDI_p5(epyz);
    const BS_REAL FDI_p3_e_ = refl ? -FDI_p3(-eta) : FDI_p3_e;
    const BS_REAL FDI_p4_e_ = refl ?  FDI_p4(-eta) : FDI_p4_e;

    const BS_REAL x0 = twenty * FDI_p4_emy;
    const BS_REAL x1 = twenty * FDI_p4_epz;
    const BS_REAL x2 = onehundredtwenty * FDI_p2_emy - onehundredtwenty * FDI_p2_epz;
    const BS_REAL x3  = -forty * FDI_p3_emy + forty * FDI_p3_epz;
    const BS_REAL x4  = forty * FDI_p3_e_;
    const BS_REAL x5  = forty * FDI_p3_emyz - x4;
    const BS_REAL x6  = -twenty * FDI_p4_e_;
    const BS_REAL x7  = twenty * FDI_p4_emyz + x6;
    const BS_REAL x8  = -forty * FDI_p3_epyz + x4;
    const BS_REAL x9  = twenty * FDI_p4_epyz + x6;
    const BS_REAL x10 = -four * FDI_p5_emy + four * FDI_p5_emyz - four * FDI_p5_emz +
                        four * FDI_p5_epy - four * FDI_p5_epyz + four * FDI_p5_epz;
    const BS_REAL x11 = twenty * FDI_p4_epy;
    const BS_REAL x12 = twenty * FDI_p4_emz;
    const BS_REAL x13 = -forty * FDI_p3_emz + forty * FDI_p3_epy;
    const BS_REAL x14 = onehundredtwenty * FDI_p2_emz - onehundredtwenty * FDI_p2_epy;

    const BS_REAL aux = fifteen * POW2(y * z);

    psi_out[0] =
        x10 +
        y * (-x0 + x1 + x7 +
             y * (x3 + x5 +
                  z * (x2 + z * (-onehundredtwenty * FDI_p1_emy +
                                 onehundredtwenty * FDI_p1_epz))) +
             z * (eighty * FDI_p3_emy - eighty * FDI_p3_epz - x2 * z)) +
        z * (x0 - x1 + x9 + z * (x3 + x8));

    psi_out[1] =
        x10 + y * (-x11 + x12 + x9 + y * (x13 + x8)) +
        z * (x11 - x12 + x7 + y * (eighty * FDI_p3_emz - eighty * FDI_p3_epy - x14 * y) +
             z * (x13 + x5 +
                  y * (x14 + y * (-onehundredtwenty * FDI_p1_emz +
                                  onehundredtwenty * FDI_p1_epy))));

    psi_out[0] /= aux;
    psi_out[1] /= aux;
}

/* Same as PairKernels(eos, kp, &out_for, &out_inv) for omega = nu_a,
 * omega_prime = nu_b, with cached FDIs and exponentials computed once. */
KOKKOS_INLINE_FUNCTION
void PairKernelsCached(const BS_REAL nu_a, const BS_REAL nu_b, const BS_REAL T,
                       const BS_REAL eta, const BS_REAL FDI_p3_e, const BS_REAL FDI_p4_e,
                       const PairNodeFDI* A, const PairNodeFDI* B,
                       MyKernelOutput* out_for, MyKernelOutput* out_inv)
{
    constexpr BS_REAL one  = 1;
    constexpr BS_REAL half = 0.5;
    const BS_REAL y = nu_a / T;
    const BS_REAL z = nu_b / T;

    // Detailed balance without forming exp(y+z) * emission (cf. the clamped-SafeExp
    // issue fixed for the beta absorptivity): with u = exp(-(y+z)) <= 1,
    //   1 / (1 - exp(y+z))         = -u / (1 - u)   (emission)
    //   exp(y+z) / (1 - exp(y+z))  = -1 / (1 - u)   (absorption)
    // u only underflows to 0 for large y+z, so no clamp is needed.
    const BS_REAL u       = Kokkos::exp(-(y + z));
    const BS_REAL abs_fac = -one / (one - u);
    const BS_REAL aux     = kBS_Pair_Phi * POW2(T) * u * abs_fac;
    const BS_REAL aux_abs = kBS_Pair_Phi * POW2(T) * abs_fac;

    BS_REAL psi[2];
    PairPsiCached(y, z, eta, FDI_p3_e, FDI_p4_e, A, B, psi);

    const BS_REAL c[total_num_species] = {
        POW2(kBS_Pair_Alpha1_0) * psi[0] + POW2(kBS_Pair_Alpha2_0) * psi[1],
        POW2(kBS_Pair_Alpha1_0) * psi[1] + POW2(kBS_Pair_Alpha2_0) * psi[0],
        POW2(kBS_Pair_Alpha1_1) * psi[0] + POW2(kBS_Pair_Alpha2_1) * psi[1],
        POW2(kBS_Pair_Alpha1_1) * psi[1] + POW2(kBS_Pair_Alpha2_1) * psi[0]};

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        out_for->em[idx]  = half * c[idx] * aux;
        out_for->abs[idx] = half * c[idx] * aux_abs;
    }

    out_inv->em[id_nue]   = out_for->em[id_anue];
    out_inv->em[id_anue]  = out_for->em[id_nue];
    out_inv->em[id_nux]   = out_for->em[id_anux];
    out_inv->em[id_anux]  = out_for->em[id_nux];
    out_inv->abs[id_nue]  = out_for->abs[id_anue];
    out_inv->abs[id_anue] = out_for->abs[id_nue];
    out_inv->abs[id_nux]  = out_for->abs[id_anux];
    out_inv->abs[id_anux] = out_for->abs[id_nux];
}

#endif // BNS_NURATES_INCLUDE_OPT_KERNELS_HPP_
