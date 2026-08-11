"""Figure S4i: radial profile of junction lifetime (intensity-weighted mean
fluorescence lifetime, tau_int, ns) from the rosette center for the single
representative E-cadherin FOV shown in Figure S4f (190-15-apical-extrusion).

Each point is the mean lifetime of junctional pixels within a 2-um-thick radial
ring; point color encodes distance from the rosette center (plasma colormap),
matching the color coding used in S4f.

Reads data/FigS4/FigS4i.csv and writes figures/FigS4i.svg and .png.
"""

import sys
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR, TAU  # noqa: E402

apply_style()

# Axes size measured from the published panel; margins are absolute so the
# canvas can be swapped in directly.
AX_W, AX_H = 2.745, 1.327
LEFT_IN, RIGHT_IN, TOP_IN, BOT_IN = 0.34, 0.05, 0.05, 0.30
FIG_W, FIG_H = AX_W + LEFT_IN + RIGHT_IN, AX_H + TOP_IN + BOT_IN


def main():
    df = pd.read_csv(DATA_DIR / "FigS4" / "FigS4i.csv").sort_values("ring_distance_um")
    r = df["ring_distance_um"].to_numpy()
    tau = df["tau_int_ns"].to_numpy()

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.plot(r, tau, "-", color="0.6", lw=0.8, zorder=1)
    sc = ax.scatter(r, tau, c=r, cmap="plasma", s=16, edgecolors="black",
                    linewidths=0.3, zorder=3)
    # no colorbar here: the published panel keys colour to the one under Fig. S4f
    ax.set_xlabel("Distance from rosette center (µm)", labelpad=1.5)
    ax.set_ylabel(f"Junction {TAU} (ns)", labelpad=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="both", labelsize=7, pad=1.5)
    fig.subplots_adjust(left=LEFT_IN / FIG_W, right=1 - RIGHT_IN / FIG_W,
                        top=1 - TOP_IN / FIG_H, bottom=BOT_IN / FIG_H)
    save(fig, "FigS4i")


if __name__ == "__main__":
    main()
