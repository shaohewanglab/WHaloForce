"""
Fig. 6f: spermathecal local lifetime contrast (tau_Spt(norm)) versus number of
ovulations per spermatheca, in young adults (CeLam-WF). One dot per
spermatheca (53 spermathecae from 42 worms).

Reads Fig6f.csv (relative path) and regenerates the scatter with a linear fit
line and Pearson correlation.

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
from style import apply_style, save, sup_text, format_p, DATA_DIR, TAU  # noqa: E402

apply_style()


X_COL = "mean_ovulations_per_spermatheca"
Y_COL = "tau_Spt_norm_ns"


def main():
    df = pd.read_csv(DATA_DIR / "Fig6" / "Fig6f.csv")
    x = df[X_COL].to_numpy(dtype=float)
    y = df[Y_COL].to_numpy(dtype=float)

    r, p = pearsonr(x, y)
    fit = linregress(x, y)

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

    # n / r / p in the empty upper-left corner. Two lines: on a 1.55 in panel the
    # single published line does not fit. The p exponent is built from plain text
    # runs so the whole string stays editable in the SVG.
    base, exp = format_p(p)
    ax.text(0.03, 0.99, f"n = {len(x)}   r = {r:.2f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=7)
    if exp is None:
        ax.text(0.03, 0.895, base, transform=ax.transAxes, ha="left", va="top",
                fontsize=7)
    else:
        sup_text(ax, 0.03, 0.895, base, exp, fontsize=7, ha="left", va="top",
                 xycoords="axes fraction", pad_pt=0.0)

    fig.subplots_adjust(left=0.30, right=0.985, top=0.985, bottom=0.135)
    save(fig, "Fig6f")
    print(f"Pearson r = {r:.4f}, p = {p:.4g}, n = {len(x)}")


if __name__ == "__main__":
    main()
