#!/usr/bin/env python3
"""
plot_cuts.py
============
Load `events_gold` and `events_uranium` from coincidences.root and produce
clear, publication-style plots of the cut variables, split by neutron-energy
region (mirroring the step table in getCutsNominal, cuts.h).

For each sample and each energy region it draws:
  1. (amp0+amp1) vs (tof1-tof0)                  -- amplitude-sum cut lines + ROI
  2. (amp1-amp0)/(amp0+amp1) vs (tof1-tof0)       -- ratio cut lines + ROI
  3. 1-D histogram of tof1-tof0                   -- ROI / U-peak windows shaded

Requires: uproot, numpy, matplotlib
    pip install uproot numpy matplotlib

Usage:
    python plot_cuts.py --infile coincidences.root --outdir plots
"""

import argparse
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import uproot

# halo used behind every cut line so it reads clearly against any part of
# the density colormap (dark or light)
_HALO = [pe.withStroke(linewidth=4.2, foreground="white")]

# ---------------------------------------------------------------------------
# Nominal cut table -- mirrors getCutsNominal() in cuts.h exactly.
# (label, e_lo, e_hi, roi_min, roi_max, amp_min, amp_max,
#  ratio_max, ratio_min, upeak_min, upeak_max)
# ---------------------------------------------------------------------------
URANIUM_CUTS = [
    ("E \u2265 1000 MeV",  1000, 2000, -5.0,  5.0,  8e3, 35e3, 0.1, -0.9, 0.0, 0.0),
    ("500 \u2013 1000 MeV", 500, 1000, -4.0,  4.0, 13e3, 35e3, 0.0, -0.9, 0.0, 0.0),
    ("100 \u2013 500 MeV",  100,  500, -4.0,  5.0, 13e3, 35e3, 0.1, -0.9, 0.0, 0.0),
    ("10 \u2013 100 MeV",    10,  100, -5.5,  7.0, 10e3, 35e3, 0.1, -0.9, 0.0, 0.0),
    ("E < 10 MeV",            0,   10, -5.0,  7.0,  8e3, 38e3, 0.1, -1.0, 0.0, 0.0),
]

GOLD_CUTS = [
    ("E \u2265 1000 MeV",  1000, 2000, -3.0,  3.0, 18e3, 37e3, 0.30, -0.7, -15.0, -3.0),
    ("600 \u2013 1000 MeV", 600, 1000, -3.0,  2.5, 20e3, 37e3, 0.40, -0.7, -15.0, -3.0),
    ("300 \u2013 600 MeV",  300,  600, -2.5,  2.5, 20e3, 35e3, 0.30, -0.6, -15.0, -3.0),
    ("150 \u2013 300 MeV",  150,  300, -2.0,  2.3, 20e3, 35e3, 0.15, -0.7, -15.0, -2.5),
    ("E < 150 MeV",          40,  150, -2.0,  2., 21e3, 35e3, 0.10, -0.6, -15.0, -2.5),
]

SAMPLES = {
    "uranium": dict(tree="events_uranium", cuts=URANIUM_CUTS,
                     cmap="Blues", accent="#1d6fa5", has_upeak=False,
                     title="Uranium"),
    "gold":    dict(tree="events_gold",    cuts=GOLD_CUTS,
                     cmap="Oranges", accent="#a45112", has_upeak=True,
                     title="Gold"),
}

# saturated, mutually distinct cut-line colors -- always drawn with a white
# halo (see _HALO) so they stay legible against the grayscale density map
ROI_COLOR   = "#2b008f"   # teal
AMP_COLOR   = "#e8630a"   # orange
UPEAK_COLOR = "#24aa24"   # violet
RATIO_COLOR = "#e8630a"   # orange (ratio_max)
RATIO2_COLOR = "#c62828"  # red    (ratio_min)


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
def set_style():
    mpl.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "#fbfbfd",
        "axes.edgecolor": "#444444",
        "axes.linewidth": 0.9,
        "axes.labelcolor": "#222222",
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelsize": 10.5,
        "axes.grid": True,
        "grid.color": "#e2e2e6",
        "grid.linewidth": 0.6,
        "grid.alpha": 0.7,
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "xtick.color": "#333333",
        "ytick.color": "#333333",
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.frameon": False,
        "legend.fontsize": 9.5,
        "figure.dpi": 120,
        "savefig.dpi": 220,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
    })


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_sample(infile, treename):
    with uproot.open(infile) as f:
        tree = f[treename]
        arrs = tree.arrays(
            ["amp0", "amp1", "tof0", "tof1", "neutron_energy",
             "cos_theta", "cos_theta_det"],
            library="np",
        )

    amp0 = arrs["amp0"].astype(float)
    amp1 = arrs["amp1"].astype(float)
    dt   = arrs["tof1"] - arrs["tof0"]
    E    = arrs["neutron_energy"]

    amp_sum   = amp0 + amp1
    amp_ratio = (amp1 - amp0) / (amp0 + amp1)

    # basic geometric quality cut, mirrors fillHistograms()
    good = (
        (arrs["cos_theta_det"] >= 0.0)
        & (np.abs(arrs["cos_theta_det"]) <= 1.0)
        & (np.abs(arrs["cos_theta"]) <= 1.0) & (amp_sum>4000)
    )

    return dict(
        amp_sum=amp_sum[good],
        amp_ratio=amp_ratio[good],
        dt=dt[good],
        E=E[good],
    )


# ---------------------------------------------------------------------------
# Per-region plotting helpers
# ---------------------------------------------------------------------------
def _region_mask(E, e_lo, e_hi):
    return (E >= e_lo) & (E < e_hi)


def plot_amp_vs_time(ax, data, cut, cmap, has_upeak):
    (label, e_lo, e_hi, roi_min, roi_max, amp_min, amp_max,
     ratio_max, ratio_min, up_min, up_max) = cut

    m = _region_mask(data["E"], e_lo, e_hi)
    dt, amp = data["dt"][m], data["amp_sum"][m]

    if dt.size == 0:
        ax.text(0.5, 0.5, "no events", ha="center", va="center",
                transform=ax.transAxes, color="#999999")
        ax.set_title(f"{label}  (n=0)")
        return

    xr = np.percentile(amp, [0.5, 99.5])
    yr = np.percentile(dt, [0.5, 99.5])
    h = ax.hist2d(amp, dt, bins=[110, 110], range=[xr, yr],
                   cmap=cmap, norm=mpl.colors.LogNorm(), cmin=1)
    cb = ax.figure.colorbar(h[3], ax=ax, pad=0.015, fraction=0.046)
    cb.ax.tick_params(labelsize=7.5)
    cb.set_label("counts", fontsize=8)

    ax.axhline(roi_min, color=ROI_COLOR, lw=1.4, ls="-", path_effects=_HALO, zorder=5)
    ax.axhline(roi_max, color=ROI_COLOR, lw=1.4, ls="-", path_effects=_HALO,
               zorder=5, label="ROI")

    ax.axvline(amp_min, color=AMP_COLOR, lw=1.4, ls="--", path_effects=_HALO, zorder=5)
    ax.axvline(amp_max, color=AMP_COLOR, lw=1.4, ls="--", path_effects=_HALO,
               zorder=5, label="amp cut")

    if has_upeak and (up_min != 0.0 or up_max != 0.0):
        ax.axhline(up_min, color=UPEAK_COLOR, lw=1.2, ls=":", path_effects=_HALO, zorder=5)
        ax.axhline(up_max, color=UPEAK_COLOR, lw=1.2, ls=":", path_effects=_HALO,
                   zorder=5, label="U-peak")

    ax.set_title(f"{label}  (n={dt.size:})")
    ax.set_xlabel(r"amp$_0$+amp$_1$")
    ax.set_ylabel(r"$t_1-t_0$ [ns]")
    return h[3]  # QuadMesh, for shared colorbar if desired


def plot_ratio_vs_time(ax, data, cut, cmap):
    (label, e_lo, e_hi, roi_min, roi_max, amp_min, amp_max,
     ratio_max, ratio_min, up_min, up_max) = cut

    m = _region_mask(data["E"], e_lo, e_hi)
    dt, ratio = data["dt"][m], data["amp_ratio"][m]

    if dt.size == 0:
        ax.text(0.5, 0.5, "no events", ha="center", va="center",
                transform=ax.transAxes, color="#999999")
        ax.set_title(f"{label}  (n=0)")
        return

    yr = np.percentile(dt, [0.5, 99.5])
    h = ax.hist2d(ratio, dt, bins=[110, 110], range=[[-1, 1], yr],
                   cmap=cmap, norm=mpl.colors.LogNorm(), cmin=1)
    cb = ax.figure.colorbar(h[3], ax=ax, pad=0.015, fraction=0.046)
    cb.ax.tick_params(labelsize=7.5)
    cb.set_label("counts", fontsize=8)

    ax.axhline(roi_min, color=ROI_COLOR, lw=1.0, ls="-", path_effects=_HALO, zorder=5)
    ax.axhline(roi_max, color=ROI_COLOR, lw=1.0, ls="-", path_effects=_HALO,
               zorder=5, label="ROI")

    ax.axvline(ratio_max, color=RATIO_COLOR, lw=1.0, ls="--", path_effects=_HALO,
               zorder=5, label="ratio_max")
    ax.axvline(ratio_min, color=RATIO2_COLOR, lw = 1.0, ls=":", path_effects=_HALO,
               zorder=5, label="ratio_min")

    ax.set_title(f"{label}  (n={dt.size:})")
    ax.set_xlabel(r"(amp$_1$-amp$_0$)/(amp$_0$+amp$_1$)")
    ax.set_ylabel(r"$t_1-t_0$ [ns]")


def plot_dt_hist(ax, data, cut, accent, has_upeak):
    (label, e_lo, e_hi, roi_min, roi_max, amp_min, amp_max,
     ratio_max, ratio_min, up_min, up_max) = cut

    m = _region_mask(data["E"], e_lo, e_hi)
    dt = data["dt"][m]

    if dt.size == 0:
        ax.text(0.5, 0.5, "no events", ha="center", va="center",
                transform=ax.transAxes, color="#999999")
        ax.set_title(f"{label}  (n=0)")
        return

    xr = np.percentile(dt, [0.2, 99.8])
    bins = np.linspace(xr[0], xr[1], 140)
    ax.hist(dt, bins=bins, color=accent, alpha=0.55, histtype="stepfilled",
             linewidth=0.0, zorder=1)
    ax.hist(dt, bins=bins, color=accent, histtype="step",
             linewidth=1.4, zorder=2)

    ax.axvline(roi_min, color=ROI_COLOR, lw=2.4, ls="-", path_effects=_HALO, zorder=5)
    ax.axvline(roi_max, color=ROI_COLOR, lw=2.4, ls="-", path_effects=_HALO,
               zorder=5, label="ROI")
    if has_upeak and (up_min != 0.0 or up_max != 0.0):
        ax.axvline(up_min, color=UPEAK_COLOR, lw=2.2, ls=":", path_effects=_HALO, zorder=5)
        ax.axvline(up_max, color=UPEAK_COLOR, lw=2.2, ls=":", path_effects=_HALO,
                   zorder=5, label="U-peak")

    ax.set_title(f"{label}  (n={dt.size:})")
    ax.set_xlabel(r"$t_1-t_0$ [ns]")
    ax.set_ylabel("events")
    ax.set_yscale("log")


# ---------------------------------------------------------------------------
# Figure assembly: one row-of-5 (energy regions) figure per variable
# ---------------------------------------------------------------------------
def make_grid_figure(data, cuts, plot_fn, suptitle, outpath, ncols=5, **kw):
    n = len(cuts)  # normalmente 5
    nrows = 2      # fijo: 2 filas
    ncols = 3      # fijo: 3 columnas
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.3 * ncols, 3.9 * nrows),
                              squeeze=False)
    axes_flat = axes.flatten()

    for ax, cut in zip(axes_flat, cuts):
        plot_fn(ax, data, cut, **kw)
        h, l = ax.get_legend_handles_labels()
        if h:
            leg = ax.legend(loc="best", handlelength=1.6)
        for text in leg.get_texts():
            text.set_fontweight("bold")

    for ax in axes_flat[len(cuts):]:
        ax.axis("off")

    fig.suptitle(suptitle, fontsize=15, fontweight="bold", x=0.51, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(outpath)
    plt.close(fig)
    print(f"  saved {outpath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--infile", default="/Users/nico/Desktop/Tese/Analysis/cross_section/data/events_selection.root",
                     help="path to coincidences.root")
    ap.add_argument("--outdir", default="plots", help="output directory")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    set_style()

    for sample_key, cfg in SAMPLES.items():
        print(f"[{sample_key}] loading tree '{cfg['tree']}' ...")
        data = load_sample(args.infile, cfg["tree"])
        print(f"[{sample_key}] {data['dt'].size:,} events after quality cuts")

        make_grid_figure(
            data, cfg["cuts"], plot_amp_vs_time,
            suptitle=f"{cfg['title']}" ,
            outpath=os.path.join(args.outdir, f"{sample_key}_amp_vs_dt.pdf"),
            cmap=cfg["cmap"], has_upeak=cfg["has_upeak"],
        )

        make_grid_figure(
            data, cfg["cuts"], plot_ratio_vs_time,
            suptitle=f"{cfg['title']} ",
            outpath=os.path.join(args.outdir, f"{sample_key}_ratio_vs_dt.pdf"),
            cmap=cfg["cmap"],
        )

        make_grid_figure(
            data, cfg["cuts"], plot_dt_hist,
            suptitle=f"{cfg['title']}",
            outpath=os.path.join(args.outdir, f"{sample_key}_dt_hist.pdf"),
            accent=cfg["accent"], has_upeak=cfg["has_upeak"],
        )

    print("done.")


if __name__ == "__main__":
    main()