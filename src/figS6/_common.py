"""Shared plotting helpers for Fig. S6d and Fig. S6e (house style)."""
import sys
from pathlib import Path

from matplotlib import pyplot as plt
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


def sem(vals):
    import numpy as np
    vals = np.asarray(vals, dtype=float)
    if len(vals) < 2:
        return 0.0
    return vals.std(ddof=1) / len(vals) ** 0.5


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
    Mann-Whitney U with Holm correction; Fig. S6d and S6e deliberately were
    not, and keep Welch's t-test. Both carry bounded, ceiling-limited
    measurements -- brood size and embryonic viability (%) -- and the N2
    reference piles up against that ceiling: 17 of 24 N2 worms in FigS6e.csv
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


def bar_with_sem_panel(df, value_col, ylabel, ymax, tick_step, out_stub, p_adj,
                       ax_w=1.18, ax_h=1.11):
    """House-style bar + SEM panel matching the published Fig. S6d/E layout:
    gray bars, black error bars, Holm-adjusted p (or 'ref') printed inside
    each bar, rotated strain labels with n= below the axis.
    """
    groups = [g for g in STRAIN_ORDER if g in set(df["strain"])]
    means = [df.loc[df["strain"] == g, value_col].mean() for g in groups]
    sems = [sem(df.loc[df["strain"] == g, value_col]) for g in groups]
    ns = [int((df["strain"] == g).sum()) for g in groups]

    # Axes size measured from the published panel; margins are absolute. The
    # bottom has to carry the 90-degree strain labels, which are long.
    left_in, right_in, top_in, bot_in = 0.34, 0.04, 0.04, 0.70
    fig_w, fig_h = ax_w + left_in + right_in, ax_h + top_in + bot_in
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    x = range(len(groups))
    ax.bar(x, means, yerr=sems, width=0.62, color="0.85", edgecolor="0.4",
           linewidth=0.8, error_kw=dict(elinewidth=0.8, capsize=2, capthick=0.8, ecolor="black"))

    for i, g in enumerate(groups):
        label = "ref" if g == CONTROL else fmt_p(p_adj[g])
        ax.text(i, means[i] * 0.5, label, ha="center", va="center", rotation=90, fontsize=6)

    ax.set_xticks(list(x))
    ax.set_xticklabels([f"{g}\n(n={n})" for g, n in zip(groups, ns)], rotation=90)
    ax.set_ylabel(ylabel, labelpad=1.5)
    ax.set_ylim(0, ymax)
    ax.set_yticks(range(0, ymax + 1, tick_step))
    ax.spines[["top", "right"]].set_visible(False)

    ax.tick_params(axis="both", labelsize=7, pad=1.5)
    fig.subplots_adjust(left=left_in / fig_w, right=1 - right_in / fig_w,
                        top=1 - top_in / fig_h, bottom=bot_in / fig_h)
    save(fig, out_stub)
