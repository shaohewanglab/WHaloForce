"""Shared plotting helpers for Extended Data Fig. 6d and Extended Data Fig. 6e (house style)."""
import sys
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.transforms import blended_transform_factory
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR, FIG_DIR  # noqa: E402

apply_style()

STRAIN_ORDER = [
    "N2", "CeEcad-WF", "CeEcad-cpH", "CeEcad-C-WF",
    "CeVin-WF", "CeVin-cpH", "CeVin-C-WF",
    "CeLam-WF", "CeLam-cpH", "CeLam-C-WF",
]
CONTROL = "N2"


def sd(vals):
    import numpy as np
    vals = np.asarray(vals, dtype=float)
    if len(vals) < 2:
        return 0.0
    return vals.std(ddof=1)


def fmt_p(p):
    if p < 0.0001:
        return "p<0.0001"
    if p < 0.001:
        return f"p={p:.4f}"
    return f"p={p:.3f}"


def welch_holm_vs_control(df, value_col, group_col="strain", control=CONTROL, order=STRAIN_ORDER):
    """Two-sided Welch's t-test of each non-control group vs control, with
    Holm-Bonferroni adjustment across the (n_groups - 1) comparisons.

    WHY WELCH AND NOT MANN-WHITNEY U: the Fig. 6 panels were unified onto
    Mann-Whitney U with Holm correction; Extended Data Fig. 6d and ED6e deliberately were
    not, and keep Welch's t-test. Both carry bounded, ceiling-limited
    measurements -- brood size and embryonic viability (%) -- and the N2
    reference piles up against that ceiling: 17 of 24 N2 worms in ExtDataFig6e.csv
    have viability of exactly 100.00%. Because a rank test sees only order, a
    strain sitting consistently at 98-99% ranks below almost every N2 worm and
    is reported as a strong shift: Mann-Whitney U with Holm gives CeVin-WF
    p = 0.0027 and CeLam-C-WF p = 0.0048, for median differences from N2 of
    1.4 and 1.6 percentage points respectively. Those p-values describe a real
    but biologically negligible offset, and on a panel whose purpose is to show
    that the sensor strains are healthy they would invite the opposite reading.
    Welch's t-test on the same comparisons returns p = 0.58 and p = 0.36.

    Neither test is ideal for ceiling-bounded proportions. Welch is retained
    because it is what the published panel reports and because it does not
    manufacture significance from the tie structure at 100%. If these panels
    are ever revisited, a test built for proportions (e.g. a beta-binomial or
    logistic model on the hatched/total counts, which are in the CSV as
    n_hatched and n_embryos) would be a better answer than either.
    """
    others = [g for g in order if g != control and g in set(df[group_col])]
    ctrl_vals = df.loc[df[group_col] == control, value_col].to_numpy()
    raw = []
    for g in others:
        vals = df.loc[df[group_col] == g, value_col].to_numpy()
        _, p = ttest_ind(vals, ctrl_vals, equal_var=False)
        raw.append(p)
    _, p_adj, _, _ = multipletests(raw, method="holm")
    return dict(zip(others, p_adj))


def bar_with_sd_panel(df, value_col, ylabel, ymax, tick_step, out_stub, p_adj,
                      ax_w=1.75, ax_h=1.11):
    """House-style bar + SD panel for Extended Data Fig. 6d/e: gray bars (mean),
    black error bars (SD; SEM before 2026-10-31), Holm-adjusted p (or 'ref')
    printed inside each bar at 7 pt, a horizontal gray n row under the axis
    sharing one "n:" at the left, and rotated strain labels below it.
    """
    groups = [g for g in STRAIN_ORDER if g in set(df["strain"])]
    means = [df.loc[df["strain"] == g, value_col].mean() for g in groups]
    sds = [sd(df.loc[df["strain"] == g, value_col]) for g in groups]
    ns = [int((df["strain"] == g).sum()) for g in groups]

    # Margins are absolute. The axes are wider than the original 1.18 in panel
    # (and the bars wider, 0.62 -> 0.78) so the 7 pt p labels fit inside the
    # bars. The bottom carries the n row plus the 90-degree strain labels and is
    # sized to the longest label after drawing.
    left_in, right_in, top_in, bot_in = 0.34, 0.04, 0.04, 0.70
    fig_w, fig_h = ax_w + left_in + right_in, ax_h + top_in + bot_in
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    x = range(len(groups))
    ax.bar(x, means, yerr=sds, width=0.78, color="0.85", edgecolor="0.4",
           linewidth=0.8, error_kw=dict(elinewidth=0.8, capsize=2, capthick=0.8, ecolor="black"))

    for i, g in enumerate(groups):
        label = "ref" if g == CONTROL else fmt_p(p_adj[g])
        ax.text(i, means[i] * 0.5, label, ha="center", va="center", rotation=90, fontsize=7)

    ax.set_xticks(list(x))
    ax.set_xticklabels(groups, rotation=90)
    # n row: horizontal, gray, one shared "n:" left of the first bar
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for i, n in enumerate(ns):
        ax.annotate(f"{n}", xy=(i, 0.0), xycoords=trans, xytext=(0.0, -8.0),
                    textcoords="offset points", ha="center", va="center",
                    fontsize=7, color="0.45", annotation_clip=False)
    ax.annotate("n:", xy=(0, 0.0), xycoords=trans, xytext=(-9.0, -8.0),
                textcoords="offset points", ha="right", va="center", fontsize=7,
                color="0.45", annotation_clip=False)
    ax.set_ylabel(ylabel, labelpad=1.5)
    # SD bars can pass a bounded ceiling (viability, 100%); extend the axis so
    # no bar or cap is clipped, but keep the ticks at or below ymax
    top = max(ymax, max(m + s for m, s in zip(means, sds)) * 1.01)
    ax.set_ylim(0, top)
    ax.set_yticks(range(0, ymax + 1, tick_step))
    ax.spines[["top", "right"]].set_visible(False)

    ax.tick_params(axis="both", labelsize=7, pad=1.5)
    ax.tick_params(axis="x", pad=14.0)          # strain labels start below the n row
    fig.subplots_adjust(left=left_in / fig_w, right=1 - right_in / fig_w,
                        top=1 - top_in / fig_h, bottom=bot_in / fig_h)
    # grow the bottom margin to reach the lowest strain label (+ 0.03 in)
    fig.canvas.draw()
    ax_y0 = ax.get_window_extent().y0
    lab_y0 = min(t.get_window_extent().y0 for t in ax.get_xticklabels())
    bot_in = (ax_y0 - lab_y0) / fig.dpi + 0.03
    fig_h = ax_h + top_in + bot_in
    fig.set_size_inches(fig_w, fig_h)
    fig.subplots_adjust(top=1 - top_in / fig_h, bottom=bot_in / fig_h)
    save(fig, out_stub)
