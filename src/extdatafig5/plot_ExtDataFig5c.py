"""
Extended Data Fig. 5c: representative single-cell junction lifetime profile around the cell
periphery (walking-box sampling), for one stretched aCat-WF cell.

Reads data/ExtDataFig5/ExtDataFig5c.csv by a relative path and regenerates the line plot:
the raw per-sample profile plus a circular moving-average smooth line, against
**arc position in degrees** (0-360, clockwise from the top-point reference), the
same convention Fig. 3k plots on its x axis. The CSV already carries the angle
as `arc_clock_deg`, so nothing has to be re-derived here. Dotted guides mark the
R/B/L quarter positions at 90/180/270 deg. No statistics (single cell,
illustrative panel).

SMOOTHING. The smooth line is a circular moving average with a window defined in
**degrees of arc**, not in samples, so it is reproducible for any cell regardless
of how densely the walking box sampled its perimeter: the sample window is
derived from the median angular step of the data. This cell is sampled every
~4.4 deg, so the 30 deg window is 7 samples wide. The smooth curve is a display
aid computed here, not a stored column.

Run:
    python plot_ExtDataFig5c.py

Writes ExtDataFig5c.svg and ExtDataFig5c.png to ../figures.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR, TAU  # noqa: E402

apply_style()

DATA_CSV = DATA_DIR / "ExtDataFig5" / "ExtDataFig5c.csv"

RAW_C = "#f4a09c"
SMOOTH_C = "#8b1a1a"
SMOOTH_DEG = 30.0          # circular moving-average width, in degrees of arc

# Quarter positions, clockwise from the top reference: right, bottom, left.
QUARTERS = {90: "R", 180: "B", 270: "L"}

# Canvas fixed to the figure slot (the published axes measured 2.027 x 1.308 in,
# but SW asked for a larger panel). Margins are absolute, so the plot region is
# whatever is left: 2.76 x 1.68 in.
FIG_W, FIG_H = 3.2, 2.2
LEFT_IN, RIGHT_IN, TOP_IN, BOT_IN = 0.35, 0.09, 0.17, 0.35   # right fits the "360" label


def circular_moving_average(y, window):
    """Moving average that wraps around the closed contour."""
    half = window // 2
    ypad = np.concatenate([y[-half:], y, y[:half]])
    return np.convolve(ypad, np.ones(window) / window, mode="valid")


def main():
    df = pd.read_csv(DATA_CSV).sort_values("arc_clock_deg")
    x = df["arc_clock_deg"].to_numpy()
    y = df["tau_int_ns"].to_numpy()

    # window in samples, from the requested width in degrees (odd, >= 3)
    step_deg = float(np.median(np.diff(x)))
    window = max(3, int(round(SMOOTH_DEG / step_deg)) | 1)
    y_smooth = circular_moving_average(y, window)
    print(f"[ExtDataFig5c] median step {step_deg:.3f} deg -> {SMOOTH_DEG:.0f} deg "
          f"smoothing window = {window} samples ({len(x)} samples total)")

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    for q in QUARTERS:
        ax.axvline(q, color="0.7", ls=":", lw=0.7, zorder=1)
    ax.plot(x, y, color=RAW_C, lw=0.7, marker="o", markersize=2.2, zorder=2,
            label="raw")
    ax.plot(x, y_smooth, color=SMOOTH_C, lw=1.3, zorder=3, label="smooth")

    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xlabel("Arc position (degrees)", labelpad=1.5)
    ax.set_ylim(1.8, 3.0)
    ax.set_yticks([2.0, 2.5, 3.0])
    ax.set_ylabel(f"{TAU} (ns)", labelpad=1.5)
    ax.tick_params(axis="both", labelsize=7, pad=1.5)
    for sp in ax.spines.values():      # framed, as published
        sp.set_visible(True)

    # R / B / L markers just inside the bottom of the frame
    for q, lab in QUARTERS.items():
        ax.annotate(lab, xy=(q, 0.0), xycoords=ax.get_xaxis_transform(),
                    xytext=(0.0, 2.0), textcoords="offset points",
                    ha="center", va="bottom", fontsize=6, color="0.35")

    ax.legend(frameon=False, fontsize=6.5, loc="upper right", handlelength=1.4,
              labelspacing=0.3, borderaxespad=0.2, handletextpad=0.5)
    ax.set_title(f"Orientation-aware {TAU} profile", fontsize=7, pad=3)

    fig.subplots_adjust(left=LEFT_IN / FIG_W, right=1 - RIGHT_IN / FIG_W,
                        top=1 - TOP_IN / FIG_H, bottom=BOT_IN / FIG_H)
    save(fig, "ExtDataFig5c")


if __name__ == "__main__":
    main()
