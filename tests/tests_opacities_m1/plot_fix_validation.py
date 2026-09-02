#!/usr/bin/env python3
"""
Validate the proposed AthenaK fix: does  eta_eq * (n_field/n_eq)  reproduce the
reconstructed (M1rec) NEPS emissivity? Run the driver in NEPS-only mode; it now
also outputs the corrected eq emissivity (eta0_eqc, eta_eqc).

Column layout (0-based): r rho T + 12 groups x 4 species
  eta0_re=3 eta_re=7 kap0_re=11 kap_re=15 eta0_eq=19 eta_eq=23 kap0_eq=27
  kap_eq=31 kaps_re=35 kaps_eq=39 eta0_eqc=43 eta_eqc=47     (nux = +2)

Usage: python plot_fix_validation.py neps.txt
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INUX = 2
C = {"eta0_re": 3, "eta_re": 7, "eta0_eq": 19, "eta_eq": 23,
     "eta0_eqc": 43, "eta_eqc": 47}


def col(d, name):
    return np.abs(d[:, C[name] + INUX])


def main():
    d = np.loadtxt(sys.argv[1], comments="#")
    o = np.argsort(d[:, 1]); d = d[o]           # sort by density
    rho, T = d[:, 1], d[:, 2]

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
    for a, (n_re, n_eq, n_eqc, ttl) in zip(
            ax,
            [("eta_re", "eta_eq", "eta_eqc", "energy emissivity"),
             ("eta0_re", "eta0_eq", "eta0_eqc", "number emissivity")]):
        a.plot(rho, col(d, n_re),  "-",  color="k",  lw=2.2, label="M1rec (truth)")
        a.plot(rho, col(d, n_eq),  "--", color="C3", lw=1.6, label="eq (uncorrected)")
        a.plot(rho, col(d, n_eqc), "-",  color="C0", lw=1.4, label="eq + fix")
        a.axvline(1e12, color="grey", ls=":", lw=1)
        a.set_xscale("log"); a.set_yscale("log")
        a.set_xlabel(r"$\rho$ [g cm$^{-3}$]"); a.set_ylabel(f"nux NEPS {ttl}")
        a.grid(True, which="both", alpha=0.3); a.legend(fontsize=9)
        a.set_title(f"nux NEPS {ttl}: does eq+fix match M1rec?")
    fig.tight_layout(); fig.savefig("fix_validation.png", dpi=140)
    print("wrote fix_validation.png")

    # numeric: how close is eq+fix to M1rec vs uncorrected eq (energy emiss)
    re, eq, eqc = col(d, "eta_re"), col(d, "eta_eq"), col(d, "eta_eqc")
    m = re > 0
    print(f"\n{'rho[g/cc]':>11} {'T':>5} {'eq/re':>7} {'(eq+fix)/re':>12}")
    for j in range(0, len(rho), max(1, len(rho) // 16)):
        if not m[j]:
            continue
        print(f"{rho[j]:11.2e} {T[j]:5.1f} {eq[j]/re[j]:7.2f} {eqc[j]/re[j]:12.2f}")
    # summary error metric over 1e10-1e13 (the atmosphere where the bump lives)
    sel = m & (rho > 1e10) & (rho < 1e13)
    err_eq  = np.abs(np.log10(eq[sel] / re[sel]))
    err_eqc = np.abs(np.log10(eqc[sel] / re[sel]))
    print(f"\nmedian |log10(ratio to M1rec)| over 1e10-1e13 g/cc:")
    print(f"  uncorrected eq : {np.median(err_eq):.3f} dex")
    print(f"  eq + fix       : {np.median(err_eqc):.3f} dex")


if __name__ == "__main__":
    main()
