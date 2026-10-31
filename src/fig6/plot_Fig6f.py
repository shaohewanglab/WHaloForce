"""
Fig. 6f: spermathecal local lifetime contrast (tau_Spt(norm)) versus number of
ovulations per spermatheca, in young adults (CeLam-WF). One dot per
spermatheca (53 spermathecae from 42 worms).

Reads Fig6f.csv (relative path) and regenerates the scatter with a linear fit
line and Pearson correlation.

The two spermathecae of one worm are not independent (11 worms contribute
both), so no independent-point Pearson p value is shown. Uncertainty in r is
a 95% percentile confidence interval from 50,000 whole-worm bootstrap
resamples (seed 20261029): the 42 worm IDs are resampled with replacement and
every measurement of a sampled worm is kept together. The interval is printed
here and reported in the figure legend.

Run:
    python plot_Fig6f.py

Writes Fig6f.svg and Fig6f.png next to this file.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import linregress, pearsonr

# --- house style ---
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR, TAU  # noqa: E402

apply_style()


X_COL = "mean_ovulations_per_spermatheca"
Y_COL = "tau_Spt_norm_ns"
WORM_COL = "worm_id"
N_BOOT = 50_000
SEED = 20261029


def worm_bootstrap_r(x, y, worms, n_boot=N_BOOT, seed=SEED):
    """95% percentile CI of Pearson r, resampling whole worms with replacement
    (all spermathecae of a sampled worm kept together)."""
    ids = np.unique(worms)
    rows = [np.flatnonzero(worms == w) for w in ids]
    rng = np.random.default_rng(seed)
    rs = np.empty(n_boot)
    for b in range(n_boot):
        idx = np.concatenate([rows[k] for k in rng.integers(0, ids.size, ids.size)])
        rs[b] = np.corrcoef(x[idx], y[idx])[0, 1]
    lo, hi = np.percentile(rs, [2.5, 97.5])
    return lo, hi, ids.size


def main():
    df = pd.read_csv(DATA_DIR / "Fig6" / "Fig6f.csv")
    x = df[X_COL].to_numpy(dtype=float)
    y = df[Y_COL].to_numpy(dtype=float)
    worms = df[WORM_COL].astype(str).to_numpy()

    r, _ = pearsonr(x, y)
    fit = linregress(x, y)
    lo, hi, n_worms = worm_bootstrap_r(x, y, worms)

    fig, ax = plt.subplots(figsize=(1.55, 2.1))
    # translucent dots: at n = 53 the low-ovulation columns overlap heavily, and
    # transparency shows where points stack up
    ax.scatter(x, y, s=18, color="0.3", alpha=0.55, zorder=3, linewidths=0)

    xs = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs, fit.intercept + fit.slope * xs, color="black", linewidth=0.8, zorder=2)

    # zero reference, as on the published panel
    ax.axhline(0.0, color="0.6", linewidth=0.6, linestyle=(0, (4, 3)), zorder=0)

    ax.set_xlabel("Number of ovulations", labelpad=1.5)
    ax.set_ylabel(f"{TAU}, spermatheca (norm) (ns)", labelpad=1.5)
    ax.set_xticks([0, 3, 6, 9])
    ax.set_ylim(-0.22, 0.62)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6])
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="both", labelsize=7, pad=1.5)

    # n / r in the empty upper-left corner (the worm-bootstrap CI is in the legend)
    ax.text(0.03, 0.99, f"n = {len(x)}   r = {r:.2f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=7)

    fig.subplots_adjust(left=0.30, right=0.985, top=0.985, bottom=0.135)
    save(fig, "Fig6f")
    print(f"Pearson r = {r:.4f}, n = {len(x)} spermathecae from {n_worms} worms; "
          f"95% CI {lo:.3f}-{hi:.3f} ({N_BOOT:,} whole-worm bootstrap resamples, "
          f"seed {SEED})")


if __name__ == "__main__":
    main()
