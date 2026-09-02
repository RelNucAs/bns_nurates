#!/usr/bin/env python3
"""
Overlay, on the same axes, the heavy-lepton (nux) NEPS absorption opacity kappa_a
and the iso-scattering opacity kappa_s, vs density -- to locate the scattering
atmosphere edge (where kappa_s overtakes kappa_a). eq (solid) vs recon (dashed).

Driver column layout (0-based): r rho T then 10 groups x 4 species:
  eta0_re=3 eta_re=7 kap0_re=11 kap_re=15 eta0_eq=19 eta_eq=23 kap0_eq=27
  kap_eq=31 kaps_re=35 kaps_eq=39      (species order nue,anue,nux,anux; nux=+2)

Usage: python plot_kappa_atmosphere.py neps.txt iso.txt
  neps.txt -> supplies kappa_a (NEPS absorption)
  iso.txt  -> supplies kappa_s (iso-energetic scattering)
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INUX = 2
KAP_A_RE, KAP_A_EQ = 15 + INUX, 31 + INUX   # energy absorption (nux)
KAPS_RE, KAPS_EQ = 35 + INUX, 39 + INUX     # scattering (nux)


def main():
    neps = np.loadtxt(sys.argv[1], comments="#")
    iso = np.loadtxt(sys.argv[2], comments="#")
    # both share the same profile/zone ordering; sort by density
    rho = neps[:, 1]; T = neps[:, 2]
    o = np.argsort(rho)
    rho, T = rho[o], T[o]
    ka_re = np.abs(neps[o, KAP_A_RE]); ka_eq = np.abs(neps[o, KAP_A_EQ])
    ks_re = np.abs(iso[o, KAPS_RE]);  ks_eq = np.abs(iso[o, KAPS_EQ])

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(rho, ka_eq, "-",  color="C3", label=r"$\kappa_a$ NEPS  eq")
    ax.plot(rho, ka_re, "--", color="C3", label=r"$\kappa_a$ NEPS  recon")
    ax.plot(rho, ks_eq, "-",  color="C0", label=r"$\kappa_s$ iso  eq")
    ax.plot(rho, ks_re, "--", color="C0", label=r"$\kappa_s$ iso  recon")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\rho$ [g cm$^{-3}$]")
    ax.set_ylabel(r"nux opacity  [cm$^{-1}$, scaled $\times 10^{7}$]")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9, loc="lower right")
    axT = ax.twinx(); axT.plot(rho, T, ":", color="k", lw=1.4); axT.set_ylabel("T [MeV]")
    ax.set_title(r"nux: NEPS $\kappa_a$ vs iso-scattering $\kappa_s$  (scattering atmosphere)")
    fig.tight_layout(); fig.savefig("kappa_atmosphere.png", dpi=140)
    print("wrote kappa_atmosphere.png")

    # crossing (kappa_s == kappa_a) => scattering-atmosphere edge, eq case
    diff = ks_eq - ka_eq
    sign = np.sign(diff)
    cross = np.where(np.diff(sign) != 0)[0]
    print("kappa_s = kappa_a crossings (eq), scattering dominates at LOWER rho below:")
    for c in cross:
        print(f"  rho ~ {rho[c]:.2e} g/cc  (T~{T[c]:.1f} MeV)")
    print(f"\n{'rho[g/cc]':>11} {'T[MeV]':>7} {'ka_eq':>10} {'ks_eq':>10} {'ks/ka_eq':>9}")
    for j in range(0, len(rho), max(1, len(rho) // 16)):
        ratio = ks_eq[j] / ka_eq[j] if ka_eq[j] > 0 else np.nan
        print(f"{rho[j]:11.2e} {T[j]:7.1f} {ka_eq[j]:10.2e} {ks_eq[j]:10.2e} {ratio:9.2f}")


if __name__ == "__main__":
    main()
