#!/usr/bin/env python3
# ---------------------------------------------------------------------------
#  plot_tof_cuts.py
#
#  Reproduces the selection in fillHistograms() and plots the TOF /
#  amplitude cuts, one panel per energy region.
#
#  Usage:
#     python3 plot_tof_cuts.py
#     python3 plot_tof_cuts.py --infile /path/to.root --outdir plots
# ---------------------------------------------------------------------------
import argparse
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import uproot


# ===========================================================================
#  VISUAL LANGUAGE
# ===========================================================================
CUT_LINES_INK = True

if CUT_LINES_INK:
    C_ROI   = "#111111"; LS_ROI   = "-";  LW_ROI   = 2.0
    C_AMP   = "#111111"; LS_AMP   = "--"; LW_AMP   = 1.5
    C_UPEAK = "#111111"; LS_UPEAK = ":";  LW_UPEAK = 1.5
    C_RAT   = "#111111"; LS_RAT   = "--"; LW_RAT   = 1.5
else:
    C_ROI   = "#4b0082"; LS_ROI   = "-";  LW_ROI   = 2.0
    C_AMP   = "#e8630a"; LS_AMP   = "--"; LW_AMP   = 1.5
    C_UPEAK = "#1f8a3b"; LS_UPEAK = ":";  LW_UPEAK = 1.5
    C_RAT   = "#e8630a"; LS_RAT   = "--"; LW_RAT   = 1.5

HALO = [pe.withStroke(linewidth=3.0, foreground="white", alpha=0.9)]

FILL_U  = "#4a90c4"
FILL_AU = "#d08030"


# ===========================================================================
#  CUT TABLE — from getCutsNominal() in cuts.h
#
#  (label, e_lo, e_hi, roi_min, roi_max, amp_min, amp_max,
#   ratio_min, ratio_max, upeak_min, upeak_max)
# ===========================================================================
URANIUM_CUTS = [
    ("E #geq 1000 MeV", 1000, 1000,   -6.0,   6.0,  14000,  40000,  -1,   0.0,   0.0,  0.0),
    ("500 -- 1000 MeV",  500, 1000,   -6.0,   6.0,  12000,  40000,  -1,   0.0,   0.0,  0.0),
    ("100 -- 500 MeV",   100,  500,   -6.0,   6.0,  10000,  40000,  -1,   0.0,   0.0,  0.0),
    ("10 -- 100 MeV",     10,  100,   -6.0,   6.0,  10000,  40000,  -1,   0.0,   0.0,  0.0),
    ("E < 10 MeV",         0,   10,   -6.0,  10.0,   7000,  40000,  -1,   0.0,   0.0,  0.0),
]

GOLD_CUTS = [
    ("E #geq 1000 MeV", 1000, 1000, -3.0, 4.2, 20000, 37000, -0.6, 0.4, -14.5, -3.5),
    ("600 -- 1000 MeV",  600, 1000, -3.5, 4.0, 20000, 37000, -0.6, 0.4, -14.5, -3.5),
    ("300 -- 600 MeV",   300,  600, -3.5, 4.0, 20000, 37000, -0.6, 0.4, -14.5, -3.5),
    ("150 -- 300 MeV",   150,  300, -3.5, 4.0, 20000, 37000, -0.6, 0.4, -14.0, -3.5),
    ("E < 150 MeV",       40,  150, -3.0, 4.0, 20000, 37000, -0.6, 0.4, -14.0, -3.5),
]

SAMPLES = {
    "uranium": dict(
        tree="events_uranium", cuts=URANIUM_CUTS,
        cmap="Blues", accent=FILL_U, fill_alpha=0.18,
        has_upeak=False, title="Uranium (n,f)",
    ),
    "gold": dict(
        tree="events_gold", cuts=GOLD_CUTS,
        cmap="Oranges", accent=FILL_AU, fill_alpha=0.18,
        has_upeak=True, title="Gold (n,#gamma)",
    ),
}

BAD_RUNS = {118771, 118789}


# ===========================================================================
#  Guard: verify cmap names before doing anything else
# ===========================================================================
def check_cmaps():
    valid = set(mpl.colormaps)
    for key, cfg in SAMPLES.items():
        name = cfg["cmap"]
        if name not in valid:
            raise ValueError(
                f"SAMPLES['{key}']['cmap'] = '{name}' is not a matplotlib "
                f"colormap. Try one of: Blues, Blues_r, Oranges, Oranges_r, "
                f"viridis, magma, inferno, plasma, cividis, Greys, Greys_r."
            )


# ===========================================================================
#  Style
# ===========================================================================
def set_style():
    mpl.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "#fcfcfd",
        "axes.edgecolor": "#3a3a3a",
        "axes.linewidth": 0.9,
        "axes.labelcolor": "#1f1f1f",
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelsize": 10,
        "axes.grid": True,
        "grid.color": "#e6e6ea",
        "grid.linewidth": 0.55,
        "grid.alpha": 0.8,
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "xtick.color": "#333333",
        "ytick.color": "#333333",
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.frameon": True,
        "legend.framealpha": 0.85,
        "legend.facecolor": "white",
        "legend.edgecolor": "#dcdcdc",
        "legend.fontsize": 8.5,
        "legend.handlelength": 1.6,
        "figure.dpi": 120,
        "savefig.dpi": 220,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
    })


# ===========================================================================
#  compute_angles — same as the C++ version
# ===========================================================================
def compute_angles(x0, y0, x1, y1, offx0, offx1, offy0, offy1):
    dx = (2 * x1 + 2.5) - (2 * x0 - 2.5)
    dy =  2 * y1 - 2 * y0 - (2 * offy1 - 2 * offy0)
    dz =  5.0
    nd = np.sqrt(dx * dx + dy * dy + dz * dz)

    cos_theta_det = np.full_like(dx, -999.0, dtype=float)
    cos_theta     = np.full_like(dx, -999.0, dtype=float)

    ok = nd > 0.0
    sth                = np.zeros_like(dx, dtype=float)
    sth[ok]            = np.sqrt(dx[ok] ** 2 + dy[ok] ** 2) / nd[ok]
    cos_theta_det[ok]  = dz / nd[ok]
    phi_det            = np.arctan2(dy[ok], dx[ok])

    nx = (-sth[ok] * np.cos(phi_det) + cos_theta_det[ok]) / np.sqrt(2.0)
    ny =  sth[ok] * np.sin(phi_det)
    nz = ( sth[ok] * np.cos(phi_det) + cos_theta_det[ok]) / np.sqrt(2.0)
    nb = np.sqrt(nx * nx + ny * ny + nz * nz)

    cos_theta[ok] = nz / nb
    return cos_theta_det, cos_theta


# ===========================================================================
#  Nominal cuts as a vectorised function of (sample, E)
# ===========================================================================
def get_cuts_nominal(sample, E):
    E = np.asarray(E, dtype=float)

    if sample == "uranium":
        roi_min = np.full_like(E, -6.0)
        roi_max = np.where(E < 10.0, 10.0, 6.0)
        amp_min = np.where(E >= 1000, 14000.0,
                  np.where(E >= 100,  10000.0,
                  np.where(E >=  10,  10000.0,  7000.0)))
        amp_max = np.full_like(E, 40000.0)
        rat_min = np.where(E >= 1000, -0.9, -0.8)
        rat_max = np.zeros_like(E)
        up_min  = np.zeros_like(E)
        up_max  = np.zeros_like(E)

    elif sample == "gold":
        roi_min = np.where(E >= 1000, -3.0, -3.5)
        roi_max = np.where(E >= 1000,  4.2,  4.0)
        amp_min = np.full_like(E, 20000.0)
        amp_max = np.full_like(E, 37000.0)
        rat_min = np.full_like(E, -0.6)
        rat_max = np.full_like(E,  0.4)
        up_min  = np.where(E >= 150, -14.5, -14.0)
        up_max  = np.full_like(E, -3.5)

    else:  # uranium_mc
        roi_min = np.where(E >= 500, -4.0,
                  np.where(E >= 100, -5.0,
                  np.where(E >=  10, -4.5, -7.0)))
        roi_max = np.where(E >= 500,  4.0,
                  np.where(E >= 100,  5.0,
                  np.where(E >=  10,  9.0,  9.0)))
        amp_min = np.where(E >= 500, 10000.0,
                  np.where(E >= 100, 11000.0,
                  np.where(E >=  10, 10000.0, 8000.0)))
        amp_max = np.where(E >= 500, 35000.0,
                  np.where(E >= 100, 35000.0,
                  np.where(E >=  10, 34000.0, 38000.0)))
        rat_min = np.full_like(E, 0.3)
        rat_max = np.where(E >= 500, 0.8, 0.9)
        up_min  = np.zeros_like(E)
        up_max  = np.zeros_like(E)

    return dict(roi_min=roi_min, roi_max=roi_max,
                amp_min=amp_min, amp_max=amp_max,
                ratio_min=rat_min, ratio_max=rat_max,
                upeak_min=up_min, upeak_max=up_max)


# ===========================================================================
#  Loading — mirrors fillHistograms()
# ===========================================================================
def load_sample(infile, treename, sample):
    with uproot.open(infile) as f:
        if treename not in f:
            raise RuntimeError(f"tree '{treename}' not found in {infile}")
        tree = f[treename]
        arrs = tree.arrays(
            ["amp0", "amp1", "tof0", "tof1", "neutron_energy",
             "x0", "y0", "x1", "y1",
             "x1_amp0", "x2_amp0", "y1_amp1", "y2_amp1",
             "RunNumber", "full_position"],
            library="np",
        )

    E      = arrs["neutron_energy"].astype(float)
    a0     = arrs["amp0"].astype(float)
    a1     = arrs["amp1"].astype(float)
    dt     = arrs["tof1"].astype(float) - arrs["tof0"].astype(float)
    amp_sum = a0 + a1
    amp_rat = np.where(amp_sum > 0, (a1 - a0) / np.maximum(amp_sum, 1e-9), 0.0)

    # --- run / position / gold quality ------------------------------------
    good  = ~np.isin(arrs["RunNumber"], list(BAD_RUNS))
    good &= (arrs["full_position"] == 1)

    if sample == "gold":
        good &= ~((arrs["x1_amp0"] - a0 + arrs["x2_amp0"] - a0) > -9000)
        good &= ~((arrs["y1_amp1"] - a1 + arrs["y2_amp1"] - a1) > -6000)

    # --- per-event amplitude + ratio cut ----------------------------------
    cuts     = get_cuts_nominal(sample, E)
    pass_sum = (amp_sum >= cuts["amp_min"]) & (amp_sum <= cuts["amp_max"])
    pass_rat = (amp_rat >= cuts["ratio_min"]) & (amp_rat <= cuts["ratio_max"])
    good    &= pass_sum & pass_rat

    # --- angular cuts (recomputed, not read from branches) ----------------
    x0 = arrs["x0"].astype(float); y0 = arrs["y0"].astype(float)
    x1 = arrs["x1"].astype(float); y1 = arrs["y1"].astype(float)

    if np.any(good):
        mean_x0 = np.mean(x0[good]); mean_x1 = np.mean(x1[good])
        mean_y0 = np.mean(y0[good]); mean_y1 = np.mean(y1[good])
    else:
        mean_x0 = mean_x1 = mean_y0 = mean_y1 = 0.0

    cos_det, cos_beam = compute_angles(x0, y0, x1, y1,
                                        mean_x0, mean_x1, mean_y0, mean_y1)
    good &= (cos_det >= 0.0)
    good &= (np.abs(cos_det)  <= 1.0)
    good &= (np.abs(cos_beam) <= 1.0)

    # --- analysis hard cut: nothing above 1 GeV enters the histograms -----
    good &= (E <= 1000.0)

    return dict(
        dt        = dt[good],
        E         = E[good],
        amp_sum   = amp_sum[good],
        amp_ratio = amp_rat[good],
    )


# ===========================================================================
#  Panel helpers
# ===========================================================================
def _region_mask(E, lo, hi):
    return (E >= lo) & (E < hi)


def _no_events(ax, label):
    ax.text(0.5, 0.5, "no events in this window",
            ha="center", va="center", transform=ax.transAxes,
            color="#aaaaaa", fontsize=10, style="italic")
    ax.set_title(f"{label}   (n=0)")


def _panel_2d(ax, x, y, xlabel, ylabel, cmap, xrange, yrange, nbin=140):
    h = ax.hist2d(x, y, bins=[nbin, nbin], range=[xrange, yrange],
                   cmap=cmap, norm=mpl.colors.LogNorm(), cmin=1)
    cb = ax.figure.colorbar(h[3], ax=ax, pad=0.015, fraction=0.046)
    cb.ax.tick_params(labelsize=7)
    cb.set_label("counts", fontsize=7.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def _draw_vline(ax, x, label=None, color="#111", ls="-", lw=1.5):
    ax.axvline(x, color=color, ls=ls, lw=lw, path_effects=HALO,
               zorder=5, label=label)


def _draw_hline(ax, y, label=None, color="#111", ls="-", lw=1.5):
    ax.axhline(y, color=color, ls=ls, lw=lw, path_effects=HALO,
               zorder=5, label=label)


def _add_legend(ax, loc="upper right"):
    h, l = ax.get_legend_handles_labels()
    if h:
        ax.legend(loc=loc)


# ===========================================================================
#  Panel renderers
# ===========================================================================
def plot_amp_vs_time(ax, data, cut, cmap, has_upeak):
    (label, e_lo, e_hi, roi_min, roi_max, amp_min, amp_max,
     rat_min, rat_max, up_min, up_max) = cut

    m = _region_mask(data["E"], e_lo, e_hi)
    amp, dt = data["amp_sum"][m], data["dt"][m]
    if dt.size == 0:
        _no_events(ax, label); return

    xr = np.percentile(amp, [0.2, 99.8])
    yr = np.percentile(dt,  [0.2, 99.8])

    _panel_2d(ax, amp, dt, r"$A_{\Sigma}=a_0+a_1$",
              r"$\Delta t = t_1 - t_0$ [ns]", cmap, xr, yr)

    _draw_vline(ax, amp_min, "amp window", C_AMP, LS_AMP, LW_AMP)
    _draw_vline(ax, amp_max, None,         C_AMP, LS_AMP, LW_AMP)
    _draw_hline(ax, roi_min, "ROI",        C_ROI, LS_ROI, LW_ROI)
    _draw_hline(ax, roi_max, None,         C_ROI, LS_ROI, LW_ROI)
    if has_upeak:
        _draw_hline(ax, up_min, "U contamination", C_UPEAK, LS_UPEAK, LW_UPEAK)
        _draw_hline(ax, up_max, None,              C_UPEAK, LS_UPEAK, LW_UPEAK)

    ax.set_title(f"{label}   (n={dt.size:,})")
    _add_legend(ax)


def plot_ratio_vs_time(ax, data, cut, cmap):
    (label, e_lo, e_hi, roi_min, roi_max, amp_min, amp_max,
     rat_min, rat_max, up_min, up_max) = cut

    m = _region_mask(data["E"], e_lo, e_hi)
    ratio, dt = data["amp_ratio"][m], data["dt"][m]
    if dt.size == 0:
        _no_events(ax, label); return

    yr = np.percentile(dt, [0.2, 99.8])

    _panel_2d(ax, ratio, dt,
              r"$R=(a_1-a_0)/(a_0+a_1)$",
              r"$\Delta t = t_1 - t_0$ [ns]",
              cmap, (-1.0, 1.0), yr)

    _draw_vline(ax, rat_min, "ratio window", C_RAT, LS_RAT, LW_RAT)
    _draw_vline(ax, rat_max, None,           C_RAT, LS_RAT, LW_RAT)
    _draw_hline(ax, roi_min, "ROI",          C_ROI, LS_ROI, LW_ROI)
    _draw_hline(ax, roi_max, None,           C_ROI, LS_ROI, LW_ROI)

    ax.set_title(f"{label}   (n={dt.size:,})")
    _add_legend(ax)


def plot_dt_hist(ax, data, cut, accent, fill_alpha, has_upeak):
    (label, e_lo, e_hi, roi_min, roi_max, amp_min, amp_max,
     rat_min, rat_max, up_min, up_max) = cut

    m = _region_mask(data["E"], e_lo, e_hi)
    dt = data["dt"][m]
    if dt.size == 0:
        _no_events(ax, label); return

    xr = np.percentile(dt, [0.1, 99.9])
    bins = np.linspace(xr[0], xr[1], 160)

    ax.axvspan(roi_min, roi_max, color=accent, alpha=fill_alpha,
               lw=0, zorder=0, label="ROI")

    ax.hist(dt, bins=bins, color=accent, alpha=0.55,
            histtype="stepfilled", linewidth=0.0, zorder=1)
    ax.hist(dt, bins=bins, color="#1f1f1f", histtype="step",
            linewidth=1.3, zorder=2)

    _draw_vline(ax, roi_min, None, C_ROI, LS_ROI, LW_ROI)
    _draw_vline(ax, roi_max, None, C_ROI, LS_ROI, LW_ROI)
    if has_upeak:
        _draw_vline(ax, up_min, "U contamination", C_UPEAK, LS_UPEAK, LW_UPEAK)
        _draw_vline(ax, up_max, None,              C_UPEAK, LS_UPEAK, LW_UPEAK)

    ax.set_yscale("log")
    ax.set_title(f"{label}   (n={dt.size:,})")
    ax.set_xlabel(r"$\Delta t = t_1 - t_0$ [ns]")
    ax.set_ylabel("events")
    _add_legend(ax)


# ===========================================================================
#  Figure assembly
# ===========================================================================
def make_grid_figure(data, cuts, plot_fn, suptitle, outpath, **kw):
    nrows, ncols = 2, 3
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4.4 * ncols, 3.9 * nrows),
                             squeeze=False)
    flat = axes.flatten()

    for ax, cut in zip(flat, cuts):
        plot_fn(ax, data, cut, **kw)

    for ax in flat[len(cuts):]:
        ax.axis("off")

    fig.suptitle(suptitle, fontsize=14, fontweight="bold", x=0.03, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(outpath)
    plt.close(fig)
    print(f"  saved {outpath}")


# ===========================================================================
#  Main
# ===========================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--infile",
        default="/Users/nico/Desktop/Tese/Analysis/cross_section/"
                "data/coincidences_final.root",
    )
    ap.add_argument("--outdir", default="plots")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    set_style()
    check_cmaps()

    for key, cfg in SAMPLES.items():
        print(f"[{key}] loading tree '{cfg['tree']}' ...")
        try:
            data = load_sample(args.infile, cfg["tree"], key)
        except Exception as e:
            print(f"[{key}] FAILED: {e}")
            continue
        print(f"[{key}] {data['dt'].size:,} events after all selections")

        make_grid_figure(
            data, cfg["cuts"], plot_amp_vs_time,
            suptitle=cfg["title"],
            outpath=os.path.join(args.outdir, f"{key}_amp_vs_dt.pdf"),
            cmap=cfg["cmap"], has_upeak=cfg["has_upeak"],
        )
        make_grid_figure(
            data, cfg["cuts"], plot_ratio_vs_time,
            suptitle=cfg["title"],
            outpath=os.path.join(args.outdir, f"{key}_ratio_vs_dt.pdf"),
            cmap=cfg["cmap"],
        )
        make_grid_figure(
            data, cfg["cuts"], plot_dt_hist,
            suptitle=cfg["title"],
            outpath=os.path.join(args.outdir, f"{key}_dt_hist.pdf"),
            accent=cfg["accent"], fill_alpha=cfg["fill_alpha"],
            has_upeak=cfg["has_upeak"],
        )

    print("done.")


if __name__ == "__main__":
    main()