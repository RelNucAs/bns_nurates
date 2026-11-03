"""
Plot the output of test_opacities_slice_timing: fraction of the time spent
by each reaction in the computation of the M1 emissivities and opacities
(m1_opacities.hpp) on a 2D slice of a BNS simulation.

Usage:
    python plot_opacities_cost.py TIMING_FILE [TIMING_FILE ...]
                                          [--out FIGURE]

Each TIMING_FILE (e.g. one per number of quadrature points) is shown in a
separate panel.
"""

import argparse
import os

import numpy as np

LABELS = {
    "none": "overhead\n(no reactions)",
    "beta": r"$\beta$ processes",
    "iso": "isoenergetic\nscattering",
    "pair": "pair\nprocesses",
    "brem": "NN brem\n(HR98)",
    "inelastic": "inelastic\nscattering",
}


def read_timing(fname):
    with open(fname) as fp:
        info = fp.readline().lstrip("# ").strip()
    data = np.genfromtxt(fname, comments="#", dtype=None, encoding="utf-8")
    t_mean = {str(row[0]): float(row[1]) for row in data}
    t_std = {str(row[0]): float(row[2]) for row in data}
    return info, t_mean, t_std


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("timing_files", nargs="+")
    parser.add_argument("--out", default=None,
                        help="output figure (default: next to the first file)")
    args = parser.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nfiles = len(args.timing_files)
    fig, axes = plt.subplots(1, nfiles, figsize=(7.5 * nfiles, 5.5),
                             constrained_layout=True, squeeze=False)

    for ax, fname in zip(axes[0], args.timing_files):
        info, t_mean, t_std = read_timing(fname)
        t_all = t_mean["all"]
        t_none = t_mean["none"]

        keys = ["none", "beta", "iso", "pair", "brem", "inelastic"]
        # net cost of each reaction (the overhead is shown separately)
        net = np.array([t_none if k == "none" else t_mean[k] - t_none for k in keys])
        err = np.array([t_std[k] if k == "none" else np.hypot(t_std[k], t_std["none"])
                        for k in keys])
        frac = 100. * net / t_all
        frac_err = 100. * err / t_all

        print(f"{fname}: {info}")
        print(f"  {'reactions':>10s} {'net time [s]':>13s} {'% of all':>9s}")
        for k, t, f in zip(keys, net, frac):
            print(f"  {k:>10s} {t:13.5e} {f:8.2f}%")
        print(f"  {'all':>10s} {t_all:13.5e} {100.:8.2f}%")
        print(f"  {'sum':>10s} {net.sum():13.5e} {100. * net.sum() / t_all:8.2f}%"
              " (overhead + net single reactions)")

        colors = ["0.6"] + [f"C{i}" for i in range(len(keys) - 1)]
        x = np.arange(len(keys))
        bars = ax.bar(x, frac, yerr=frac_err, color=colors, capsize=3)
        for b, f in zip(bars, frac):
            ax.annotate(f"{f:.1f}%", (b.get_x() + b.get_width() / 2, b.get_height()),
                        xytext=(0, 4), textcoords="offset points", ha="center",
                        fontsize=9)
        ax.set_xticks(x)
        ax.set_xticklabels([LABELS[k] for k in keys], fontsize=9)
        ax.set_ylabel(r"$(t_{\rm reaction} - t_{\rm none}) / t_{\rm all}$  [%]")
        ax.set_ylim(0, max(frac.max() * 1.18, 5))
        ax.set_title(f"{info}\n" + rf"$t_{{\rm all}}$ = {t_all:.3e} s, "
                     f"sum of bars = {100. * net.sum() / t_all:.1f}%", fontsize=9)

    out = args.out or os.path.join(os.path.dirname(os.path.abspath(args.timing_files[0])),
                                   "opacities_slice_timing.png")
    fig.savefig(out, dpi=150)
    print(f"Written {out}")


if __name__ == "__main__":
    main()
