"""Extended Data Fig. 1b: relative brightness of HaloTag vs cpHalo vs cpWHalo across 17 JF
dyes in live HEK293 cells.

Relative brightness (rel_brightness) is the dye-channel intensity normalized to
the co-expressed mTurquoise2 reference: for each FOV both the dye channel
(HyD X 4) and the mTq2 channel (HyD X 2) are normalized as
Sum_Intensity_kCnts * speed_Hz / power / accumulation, and
rel_brightness = dye_int / mTq2_int. See data/ExtDataFig1/ExtDataFig1b.csv (values
pre-computed from the source spreadsheet).

Per-FOV points (jittered) with mean +/- SD, faceted per dye (ordered by dye
excitation maximum), log y axis. Reads data/ExtDataFig1/ExtDataFig1b.csv and writes
figures/ExtDataFig1b.svg and .png.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR  # noqa: E402

apply_style()

ORDER = ["HaloTag", "cpHalo", "cpWHalo"]
VALUE = "rel_brightness"
YLABEL = "relative brightness (dye / mTq2)"


def dye_order(dyes):
    return sorted(dyes, key=lambda d: int(d[-3:]))


# Per-subplot axes size and the gaps between them, measured from the published
# panels; margins are absolute so the grid can be swapped in directly.
LEFT_IN, RIGHT_IN, TOP_IN, BOT_IN = 0.47, 0.05, 0.05, 0.45


def facet_panel(csv_name, value, ylabel, stem, logy,
                ax_w=0.637, ax_h=0.983, wgap=0.3345, hgap=0.4305):
    df = pd.read_csv(DATA_DIR / "ExtDataFig1" / csv_name)
    dyes = dye_order(df["dye"].unique())
    ncol = 6
    nrow = int(np.ceil(len(dyes) / ncol))
    fig_w = LEFT_IN + ncol * ax_w + (ncol - 1) * wgap + RIGHT_IN
    fig_h = TOP_IN + nrow * ax_h + (nrow - 1) * hgap + BOT_IN
    fig, axes = plt.subplots(nrow, ncol, figsize=(fig_w, fig_h),
                             squeeze=False, sharey=True)
    rng = np.random.default_rng(0)

    for k, dye in enumerate(dyes):
        ax = axes[k // ncol][k % ncol]
        sub = df[df.dye == dye]
        for i, g in enumerate(ORDER):
            y = sub.loc[sub.construct == g, value].dropna().to_numpy()
            jx = i + (rng.random(len(y)) - 0.5) * 0.3
            ax.scatter(jx, y, s=10, color="blue", alpha=0.6,
                       linewidths=0, zorder=3)
            if len(y):
                m = np.mean(y)
                sd = np.std(y, ddof=1) if len(y) > 1 else 0.0
                ax.plot([i - 0.22, i + 0.22], [m, m], color="k", lw=0.9,
                        zorder=4)
                ax.errorbar(i, m, yerr=sd, color="k", lw=0.8, capsize=2.5,
                            capthick=0.8, fmt="none", zorder=4)
        if logy:
            ax.set_yscale("log")
        # dye name inside the panel, bottom-left, as published (a title would
        # eat into the row gap and push the top row off the canvas)
        ax.text(0.04, 0.03, dye, transform=ax.transAxes, ha="left",
                va="bottom", fontsize=7)
        ax.set_xticks(range(len(ORDER)))
        last_in_col = k + ncol >= len(dyes)
        ax.set_xticklabels(ORDER if last_in_col else [], rotation=45,
                           ha="right", fontsize=6)
        ax.set_xlim(-0.6, len(ORDER) - 0.4)
        ax.spines[["top", "right"]].set_visible(False)
        if k % ncol != 0:
            ax.tick_params(axis="y", labelleft=False)

    for k in range(len(dyes), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")

    fig.subplots_adjust(left=LEFT_IN / fig_w, right=1 - RIGHT_IN / fig_w,
                        top=1 - TOP_IN / fig_h, bottom=BOT_IN / fig_h,
                        wspace=wgap / ax_w, hspace=hgap / ax_h)
    fig.supylabel(ylabel, fontsize=7, x=0.012)
    save(fig, stem)


def main():
    facet_panel("ExtDataFig1b.csv", VALUE, YLABEL, "ExtDataFig1b", logy=True)


if __name__ == "__main__":
    main()
