#!/usr/bin/env python3
"""
Check whether the NEPS eq-vs-recon emissivity "bump" (~1e12 g/cc, x-axis) tracks
a temperature feature. Consumes the driver output which now begins with
    r[km]  rho[g/cc]  T[MeV]  then 8 groups x 4 species.
Group starts (0-based) after the 3 leading cols:
    eta0_re=3 eta_re=7 kap0_re=11 kap_re=15 eta0_eq=19 eta_eq=23 kap0_eq=27 kap_eq=31
species order within a group: nue, anue, nux, anux  (nux index = 2)

Usage: python plot_neps_bump.py neps.txt [all.txt]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INUX = 2
G = {"eta0": 3, "eta": 7, "eta0_eq": 19, "eta_eq": 23}  # nux columns = start+INUX


def load(fn):
    d = np.loadtxt(fn, comments="#")
    return d


def main():
    neps = load(sys.argv[1])
    r, rho, T = neps[:, 0], neps[:, 1], neps[:, 2]
    # sort by radius for line plots
    o = np.argsort(r)
    r, rho, T, neps = r[o], rho[o], T[o], neps[o]

    eta0_re = neps[:, G["eta0"] + INUX]
    eta_re = neps[:, G["eta"] + INUX]
    eta0_eq = neps[:, G["eta0_eq"] + INUX]
    eta_eq = neps[:, G["eta_eq"] + INUX]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))

    # ---- Panel A: vs radius, with T overlaid ----
    ax = axes[0]
    ax.plot(r, np.abs(eta_eq), "-", color="C3", label="energy emiss  eq")
    ax.plot(r, np.abs(eta_re), "--", color="C3", label="energy emiss  recon")
    ax.plot(r, np.abs(eta0_eq), "-", color="C0", label="number emiss  eq")
    ax.plot(r, np.abs(eta0_re), "--", color="C0", label="number emiss  recon")
    ax.set_yscale("log"); ax.set_xlabel("radius [km]")
    ax.set_ylabel("NEPS nux emissivity  (scaled)")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8, loc="upper right")
    axT = ax.twinx()
    axT.plot(r, T, ":", color="k", lw=1.6, label="T [MeV]")
    axT.set_ylabel("T [MeV]"); axT.legend(fontsize=8, loc="lower left")
    ax.set_title("NEPS nux emissivity & T vs radius")

    # ---- Panel B: vs density (Albino's axis), with T overlaid ----
    ax = axes[1]
    od = np.argsort(rho)
    ax.plot(rho[od], np.abs(eta_eq[od]), "-", color="C3", label="energy emiss  eq")
    ax.plot(rho[od], np.abs(eta_re[od]), "--", color="C3", label="energy emiss  recon")
    ax.plot(rho[od], np.abs(eta0_eq[od]), "-", color="C0", label="number emiss  eq")
    ax.plot(rho[od], np.abs(eta0_re[od]), "--", color="C0", label="number emiss  recon")
    ax.axvline(1e12, color="grey", ls=":", lw=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\rho$ [g cm$^{-3}$]"); ax.set_ylabel("NEPS nux emissivity  (scaled)")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8, loc="upper left")
    axT = ax.twinx()
    axT.plot(rho[od], T[od], ":", color="k", lw=1.6)
    axT.set_ylabel("T [MeV]")
    ax.set_title(r"NEPS nux emissivity & T vs $\rho$  (dotted @ $10^{12}$)")

    fig.tight_layout()
    fig.savefig("neps_bump_vs_T.png", dpi=140)
    print("wrote neps_bump_vs_T.png")

    # ---- quick numeric: locate the eq/recon energy-emissivity peak & the T there ----
    ratio = np.where(np.abs(eta_re) > 0, np.abs(eta_eq) / np.abs(eta_re), np.nan)
    k = np.nanargmax(ratio)
    print(f"peak eq/recon (energy emiss, nux): {ratio[k]:.2f} at "
          f"r={r[k]:.1f} km, rho={rho[k]:.2e} g/cc, T={T[k]:.1f} MeV")
    # print a table around 1e11-1e13 g/cc
    m = (rho > 3e10) & (rho < 3e13)
    print(f"\n{'rho[g/cc]':>11} {'T[MeV]':>7} {'eta_eq/re':>10} {'eta0_eq/re':>11}")
    idx = np.argsort(rho[m])[::-1]
    rr, TT = rho[m][idx], T[m][idx]
    e_eq, e_re = np.abs(eta_eq[m][idx]), np.abs(eta_re[m][idx])
    n_eq, n_re = np.abs(eta0_eq[m][idx]), np.abs(eta0_re[m][idx])
    for j in range(0, len(rr), max(1, len(rr) // 18)):
        er = e_eq[j] / e_re[j] if e_re[j] > 0 else np.nan
        nr = n_eq[j] / n_re[j] if n_re[j] > 0 else np.nan
        print(f"{rr[j]:11.2e} {TT[j]:7.1f} {er:10.2f} {nr:11.2f}")


if __name__ == "__main__":
    main()
