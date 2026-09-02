#!/usr/bin/env python3
"""
Plot the equilibrium-vs-reconstructed neutrino opacity comparison produced by
test_opacities_bns_eq_vs_recon (the Chiesa+25 / arXiv:2412.04570 Fig. 8 check).

Pipeline:
    ./test_opacities_bns_eq_vs_recon inputs/BNS/1d/DD2_..._x_axis.txt > eqrecon.txt
    python plot_eq_vs_recon.py eqrecon.txt

Produces:
    eq_vs_recon_opacities.png  -- 2x2 panels: number emissivity, energy
                                  emissivity, stimulated number opacity,
                                  stimulated energy opacity; reconstructed
                                  (solid) vs equilibrium (dashed), per species.
    eq_vs_recon_ratio.png      -- eq/recon ratio for each quantity, nux only
                                  (the ~2x luminosity question).

Column layout written by the driver (0-based; col 0 = radius [km]):
    eta0_re[4] eta_re[4] kap0_re[4] kap_re[4] eta0_eq[4] eta_eq[4] kap0_eq[4] kap_eq[4]
    species order within each group of 4: nue, anue, nux, anux
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SPECIES = ["nue", "anue", "nux", "anux"]
# quantity groups and their starting column (after col 0 = r)
GROUPS = {  # label -> (recon_start, eq_start)
    "Number emissivity  $\\eta_0$": (1, 17),
    "Energy emissivity  $\\eta$":   (5, 21),
    "Stim. opacity (number)  $\\kappa_{0,a}$": (9, 25),
    "Stim. opacity (energy)  $\\kappa_a$":     (13, 29),
}
COLORS = {"nue": "C0", "anue": "C1", "nux": "C2", "anux": "C3"}


def main():
    infile = sys.argv[1] if len(sys.argv) > 1 else "eqrecon.txt"
    d = np.loadtxt(infile, comments="#")
    r = d[:, 0]
    order = np.argsort(r)
    r = r[order]
    d = d[order]

    # ---- Figure 1: 2x2 panels, recon (solid) vs eq (dashed) ----
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for ax, (label, (re0, eq0)) in zip(axes.flat, GROUPS.items()):
        for s, sp in enumerate(SPECIES):
            if sp == "anux":       # anux == nux for this profile; skip clutter
                continue
            yre = d[:, re0 + s]
            yeq = d[:, eq0 + s]
            ax.plot(r, np.abs(yre), "-", color=COLORS[sp], label=f"{sp} recon")
            ax.plot(r, np.abs(yeq), "--", color=COLORS[sp], label=f"{sp} eq")
        ax.set_yscale("log")
        ax.set_title(label)
        ax.grid(True, which="both", alpha=0.3)
    for ax in axes[-1]:
        ax.set_xlabel("radius [km]")
    axes[0, 0].legend(fontsize=8, ncol=3)
    fig.suptitle("bns_nurates: reconstructed (solid) vs equilibrium (dashed)")
    fig.tight_layout()
    fig.savefig("eq_vs_recon_opacities.png", dpi=140)
    print("wrote eq_vs_recon_opacities.png")

    # ---- Figure 2: eq/recon ratio for nux (the ~2x question) ----
    inux = SPECIES.index("nux")
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    for label, (re0, eq0) in GROUPS.items():
        yre = d[:, re0 + inux]
        yeq = d[:, eq0 + inux]
        ratio = np.where(np.abs(yre) > 0, yeq / yre, np.nan)
        ax2.plot(r, ratio, label=label)
    ax2.axhline(1.0, color="k", lw=0.8, ls=":")
    ax2.axhline(2.0, color="grey", lw=0.8, ls=":")
    ax2.set_xlabel("radius [km]")
    ax2.set_ylabel("equilibrium / reconstructed  (nux)")
    ax2.set_title("nux: eq/recon ratio  (emissivity ~2, opacity ~1  => physical enhancement)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8)
    fig2.tight_layout()
    fig2.savefig("eq_vs_recon_ratio.png", dpi=140)
    print("wrote eq_vs_recon_ratio.png")


if __name__ == "__main__":
    main()
