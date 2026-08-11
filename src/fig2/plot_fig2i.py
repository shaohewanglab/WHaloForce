"""Figure 2i: Vin-cpH representative FOV (LV186-003) dual-axis time trace of
bright-costamere lifetime (tau_brightFA_ns, construct color, left axis) and
beating/motion (grey, right axis).

Reads data/Fig2/Fig2i.csv and writes figures/Fig2i.svg and .png.
"""

import sys
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR  # noqa: E402

apply_style()

CPH_C, MOT_C = "#DD8452", "0.45"
CONSTRUCT = "Vin-cpH"
YLIM = (3.52, 3.60)


def main():
    df = pd.read_csv(DATA_DIR / "Fig2" / "Fig2i.csv")
    t = df["t_s"].to_numpy()
    tau = df["tau_brightFA_ns"].to_numpy()
    mot = df["motion_au"].to_numpy()

    fig, ax = plt.subplots(figsize=(3.3, 1.0))
    ax2 = ax.twinx()
    ax2.plot(t, mot, color=MOT_C, lw=0.8, alpha=0.8, zorder=2)
    ax.plot(t, tau, color=CPH_C, lw=1.0, alpha=0.8, zorder=3)
    ax.set_zorder(ax2.get_zorder() + 1)      # lifetime drawn on top
    ax.patch.set_visible(False)

    tmax = float(t.max())
    ax.set_xlim(0, tmax)
    ax.set_ylim(*YLIM)
    ax.yaxis.set_major_locator(MultipleLocator(0.04))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.tick_params(axis="y", colors=CPH_C, labelsize=6.5)
    ax.set_title(f"{CONSTRUCT} lifetime (ns)", loc="left", fontsize=6.5,
                 color=CPH_C, pad=2)

    ax2.set_ylabel("Beating (a.u.)", fontsize=7, color=MOT_C, rotation=270,
                   labelpad=9)
    ax2.tick_params(axis="y", colors=MOT_C, labelsize=6.5)
    ax2.yaxis.set_major_locator(MultipleLocator(4))

    ticks = list(range(0, int(round(tmax)) + 1, 10))
    labels = [str(v) for v in ticks]
    labels[-1] = f"{ticks[-1]} s"
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.tick_params(axis="x", labelsize=6.5, pad=1.5)

    ax.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax2.spines["right"].set_visible(True)
    ax2.spines["right"].set_color(MOT_C)
    ax.spines["left"].set_color(CPH_C)

    fig.subplots_adjust(left=0.11, right=0.88, top=0.86, bottom=0.20)
    save(fig, "Fig2i")


if __name__ == "__main__":
    main()
