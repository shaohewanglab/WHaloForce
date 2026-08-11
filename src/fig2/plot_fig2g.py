"""Figure 2g: baseline-normalized in-FA lifetime change (delta tau) during cell
spreading for Vin-WF, comparing DMSO vs Y-27632 (5 uM),
over the first 125 min (0-25 frames at 5 min/frame).

Reads data/Fig2/Fig2g.csv and writes figures/Fig2g.svg and .png. Shows the
per-condition mean trajectory +/- SEM band across FOVs. Per-FOV trapezoidal AUC
of delta tau vs time is compared Control vs Y-27632 by a single two-sided
Mann-Whitney U test (p printed).
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR, TAU, DELTA  # noqa: E402

apply_style()

# np.trapz was renamed np.trapezoid in NumPy 2.0 and removed in 2.4. Same
# function, same result; this keeps the panel working on either.
trapezoid = getattr(np, "trapezoid", None) or np.trapz

CONDITIONS = ["DMSO", "Y-27632"]
COLORS = {"DMSO": "#1f77b4", "Y-27632": "#9467bd"}


def main():
    df = pd.read_csv(DATA_DIR / "Fig2" / "Fig2g.csv")

    fig, ax = plt.subplots(figsize=(2.2, 1.7))
    n_fov = {}
    for cond in CONDITIONS:
        sub = df[df.condition == cond]
        n_fov[cond] = sub.fov_id.nunique()
        wide = sub.pivot_table(index="fov_id", columns="time_min",
                               values="dtau_ns")
        times = wide.columns.to_numpy(dtype=float)
        mat = wide.to_numpy()
        n = np.sum(~np.isnan(mat), axis=0)
        mean = np.nanmean(mat, axis=0)
        sd = np.nanstd(mat, axis=0, ddof=1)
        sem = np.where(n > 0, sd / np.sqrt(np.maximum(n, 1)), np.nan)
        th = times / 60.0
        ax.fill_between(th, mean - sem, mean + sem, color=COLORS[cond],
                        alpha=0.22, lw=0, zorder=1)
        ax.plot(th, mean, color=COLORS[cond], lw=1.0, zorder=2,
                label=f"{cond} (n={n_fov[cond]})")

    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel(f"{DELTA}{TAU}, in FA (ns)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=7, loc="best")
    fig.tight_layout(pad=0.3)
    save(fig, "Fig2g")

    # ---- per-FOV trapezoidal AUC, DMSO vs Y-27632 (Mann-Whitney U) ----
    auc = {c: [] for c in CONDITIONS}
    for (cond, fov), g in df.groupby(["condition", "fov_id"]):
        g = g.sort_values("time_min")
        auc[cond].append(trapezoid(g.dtau_ns.to_numpy(), g.time_min.to_numpy()))
    c = np.asarray(auc["DMSO"])
    y = np.asarray(auc["Y-27632"])
    u, p = mannwhitneyu(c, y, alternative="two-sided")
    print(f"[Fig2g] AUC DMSO (n={len(c)}, mean={c.mean():.2f}) vs "
          f"Y-27632 (n={len(y)}, mean={y.mean():.2f}): MWU U={u:.1f}, p={p:.4g}")


if __name__ == "__main__":
    main()
