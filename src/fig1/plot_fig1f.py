"""Figure 1f: representative cpHalo-JFX554 force-ramp trace.

Smoothed force (left, black) and fluorescence intensity (right, red) versus time for a
single tether, reproducing the dual-axis layout of the original
force_trace_simple_cli.py. Reads data/Fig1/Fig1f.csv and writes
figures/Fig1f.svg and .png.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR  # noqa: E402

apply_style()

# Display smoothing matches force_trace_simple_cli.py defaults.
FORCE_SMOOTH_S = 0.05          # force rolling-mean window (seconds)
INTENSITY_SMOOTH_PTS = 15      # intensity rolling-mean window (points)


def _smooth(x, window):
    if window <= 1:
        return np.asarray(x, dtype=float)
    minp = max(1, window // 3)
    return pd.Series(x).rolling(window, center=True, min_periods=minp).mean().to_numpy()


def plot_trace(csv_name, stem, smoothed_only=True,
               force_ticks=None, force_lim=None,
               inten_ticks=None, inten_lim=None):
    df = pd.read_csv(DATA_DIR / "Fig1" / csv_name)
    t = df["time_s"].to_numpy(float)
    force = df["force_pN"].to_numpy(float)
    inten = df["intensity_raw"].to_numpy(float)

    dt = float(np.median(np.diff(t)))
    w_force = max(1, int(round(FORCE_SMOOTH_S / max(dt, 1e-9))))
    force_sm = _smooth(force, w_force)
    inten_sm = _smooth(inten, INTENSITY_SMOOTH_PTS)

    fig, axf = plt.subplots(figsize=(2.5, 1.2), constrained_layout=True)
    if not smoothed_only:
        axf.plot(t, force, color="black", alpha=0.25, lw=0.7)
    axf.plot(t, force_sm, color="black", lw=1.0, zorder=3)
    axf.set_ylabel("Force (pN)")
    axf.set_xlabel("Time (s)")
    axf.set_xlim(min(0.0, float(np.nanmin(t))), float(np.nanmax(t)))
    if force_ticks is not None:
        axf.set_yticks(force_ticks)
        axf.set_ylim(force_lim)
    # Boxed frame with light grid behind the traces, as in the published panel.
    for sp in axf.spines.values():
        sp.set_visible(True)
    axf.grid(True, which="major", color="0.85", lw=0.5, zorder=0)
    axf.set_axisbelow(True)

    axi = axf.twinx()
    if not smoothed_only:
        axi.plot(t, inten, color="red", alpha=0.25, lw=0.7)
    axi.plot(t, inten_sm, color="red", lw=1.0, zorder=3)
    axi.set_ylabel("Intensity (a.u.)", color="red")
    axi.tick_params(axis="y", colors="red")
    if inten_ticks is not None:
        axi.set_yticks(inten_ticks)
        axi.set_ylim(inten_lim)
    axi.spines["right"].set_color("red")
    axi.grid(False)          # gridlines come from the force axis only

    save(fig, stem)


def main():
    plot_trace("Fig1f.csv", "Fig1f",
               force_ticks=[0, 5, 10], force_lim=(-1.2, 11.2),
               inten_ticks=[0, 20, 40], inten_lim=(-2.0, 44.0))


if __name__ == "__main__":
    main()
