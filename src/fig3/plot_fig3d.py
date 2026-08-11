"""Figure 3d: E-cadherin radial profile of junction lifetime change (Delta tau)
vs distance from the extrusion (rosette) center.

Delta tau is computed per FOV as the per-ring mean junction lifetime
(intensity-weighted mean tau_int, ns) minus that FOV's baseline lifetime,
where the baseline is the mean of rings at r >= 10 um (falling back to the
outer third of rings if fewer than two qualify). This mirrors the original
plot_per_group_dtau_small.py; the committed CSV stores the raw per-ring
tau_int_ns, and the baseline / Delta tau are computed here. Per-construct
mean +/- SD bands are drawn across FOVs (rings shown only where >=2 FOVs
contribute).

Reads data/Fig3/Fig3d.csv and writes figures/Fig3d.svg and .png.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR  # noqa: E402

apply_style()

ORDER = ["Ecad-WF", "Ecad-cpH", "Ecad-C-WF"]
COLORS = {"Ecad-WF": "#2c7fb8", "Ecad-cpH": "#e6842e", "Ecad-C-WF": "#d62728"}
BASELINE_FROM_UM = 10.0
MIN_FOVS_PER_RING = 5


def baseline_per_fov(df):
    out = {}
    for img, g in df.groupby("fov_id"):
        outer = g[g["ring_distance_um"] >= BASELINE_FROM_UM]
        if len(outer) < 2:
            outer = g[g["ring_distance_um"] >= g["ring_distance_um"].quantile(0.66)]
        out[img] = outer["tau_int_ns"].mean() if len(outer) else np.nan
    return pd.Series(out)


def main():
    df = pd.read_csv(DATA_DIR / "Fig3" / "Fig3d.csv")
    base = baseline_per_fov(df)
    df["dtau"] = df["tau_int_ns"] - df["fov_id"].map(base)
    half_min = max(2, MIN_FOVS_PER_RING // 2)

    fig, ax = plt.subplots(figsize=(1.9, 1.9))
    for c in ORDER:
        gg = df[df.construct == c]
        a = (gg.groupby("ring_distance_um")["dtau"]
             .agg(["mean", "std", "count"]).reset_index())
        a = a[a["count"] >= half_min]
        col = COLORS[c]
        ax.fill_between(a["ring_distance_um"], a["mean"] - a["std"],
                        a["mean"] + a["std"], color=col, alpha=0.22, lw=0)
        ax.plot(a["ring_distance_um"], a["mean"], "-", lw=0.9, color=col,
                label=f"{c} (n={gg.fov_id.nunique()})")
    ax.axhline(0, color="black", lw=0.5, ls="--", alpha=0.6)
    ax.set_xlabel("Distance from center (µm)")
    ax.set_ylabel(r"Junction $\Delta\tau$ (ns)")
    ax.legend(frameon=False, fontsize=6, handlelength=1.2, handletextpad=0.4,
              borderpad=0.2, labelspacing=0.25, loc="best")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(pad=0.4)
    save(fig, "Fig3d")


if __name__ == "__main__":
    main()
