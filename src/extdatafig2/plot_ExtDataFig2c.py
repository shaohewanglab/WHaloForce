"""Extended Data Fig. 2c: representative cpWHalo-JFX554 reversibility trace.

Force (left, black) and fluorescence intensity (right, red) versus time for a
single tether cycled through force, showing reversible quenching/unquenching.
Smoothed-only display, matching the original force_trace_simple_cli.py output.
Reads data/ExtDataFig2/ExtDataFig2c.csv and writes figures/ExtDataFig2c.svg and .png.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR  # noqa: E402

apply_style()

# Axes size measured from the published panel; margins are absolute so the
# canvas can be swapped in directly.
# The published panel also carries a Quenched/Unquenched state bar and a red
# dashed intensity threshold. Both were set empirically during figure
# assembly and are not reproducible from the committed CSV, so only the
# plot itself is regenerated here (agreed with SW).
AX_W, AX_H = 5.280, 1.288
LEFT_IN, RIGHT_IN, TOP_IN, BOT_IN = 0.36, 0.46, 0.05, 0.32
FIG_W, FIG_H = AX_W + LEFT_IN + RIGHT_IN, AX_H + TOP_IN + BOT_IN

FORCE_SMOOTH_S = 0.05
INTENSITY_SMOOTH_PTS = 15


def _smooth(x, window):
    if window <= 1:
        return np.asarray(x, dtype=float)
    minp = max(1, window // 3)
    return pd.Series(x).rolling(window, center=True, min_periods=minp).mean().to_numpy()


def plot_trace(csv_name, stem, smoothed_only=True):
    df = pd.read_csv(DATA_DIR / "ExtDataFig2" / csv_name)
    t = df["time_s"].to_numpy(float)
    force = df["force_pN"].to_numpy(float)
    inten = df["intensity_raw"].to_numpy(float)

    dt = float(np.median(np.diff(t)))
    w_force = max(1, int(round(FORCE_SMOOTH_S / max(dt, 1e-9))))
    force_sm = _smooth(force, w_force)
    inten_sm = _smooth(inten, INTENSITY_SMOOTH_PTS)

    fig, axf = plt.subplots(figsize=(FIG_W, FIG_H))
    if not smoothed_only:
        axf.plot(t, force, color="black", alpha=0.25, lw=0.7)
    axf.plot(t, force_sm, color="black", lw=1.0)
    axf.set_ylabel("Force (pN)", labelpad=1.5)
    axf.set_xlabel("Time (s)", labelpad=1.5)
    axf.set_xlim(min(0.0, float(np.nanmin(t))), float(np.nanmax(t)))
    axf.spines["top"].set_visible(False)

    axi = axf.twinx()
    if not smoothed_only:
        axi.plot(t, inten, color="red", alpha=0.25, lw=0.7)
    axi.plot(t, inten_sm, color="red", lw=1.0)
    axi.set_ylabel("Intensity (a.u.)", color="red", labelpad=1.5)
    axi.tick_params(axis="y", colors="red")
    axi.spines["top"].set_visible(False)
    axi.spines["right"].set_color("red")

    for a in (axf, axi):
        a.tick_params(axis="both", labelsize=7, pad=1.5)
    fig.subplots_adjust(left=LEFT_IN / FIG_W, right=1 - RIGHT_IN / FIG_W,
                        top=1 - TOP_IN / FIG_H, bottom=BOT_IN / FIG_H)
    save(fig, stem)


def main():
    plot_trace("ExtDataFig2c.csv", "ExtDataFig2c")


if __name__ == "__main__":
    main()
