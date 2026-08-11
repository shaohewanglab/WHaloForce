"""
Shared statistics for the Fig. 6 panels: two-sided Mann-Whitney U with
Holm-Bonferroni control of family-wise error.

Every two-group comparison in Fig. 6 uses the same test, so it lives here
rather than being re-implemented per panel. The Holm implementation matches
src/fig2/plot_fig2c.py.
"""

import numpy as np
from scipy.stats import mannwhitneyu


def holm(pvals):
    """Holm-Bonferroni step-down adjusted p-values, order-preserving.

    Sort ascending, scale the i-th smallest by (m - i), enforce monotonicity
    with a running maximum, clip at 1.
    """
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    adj = np.empty(p.size)
    run = 0.0
    for i, idx in enumerate(order):
        run = max(run, (p.size - i) * p[idx])
        adj[idx] = min(run, 1.0)
    return adj


def mwu(a, b):
    """Two-sided Mann-Whitney U on samples a, b. Returns (U, p, effect).

    `effect` is the common-language effect size (probability of superiority),
    U / (n_a * n_b): the chance a randomly drawn value from a exceeds one from
    b, ties counted as half. 0.5 is no effect, 1.0 complete separation.
    """
    res = mannwhitneyu(a, b, alternative="two-sided")
    return res.statistic, res.pvalue, res.statistic / (a.size * b.size)


def mwu_holm(pairs):
    """Run mwu() on each (a, b) in `pairs` and Holm-adjust across all of them.

    Returns a list of (U, p_raw, p_holm, effect), one per input pair, in the
    input order. The caller decides which comparisons belong in one family by
    choosing what to pass in a single call.
    """
    stats = [mwu(a, b) for a, b in pairs]
    p_adj = holm([p for _, p, _ in stats])
    return [(U, p, pa, eff) for (U, p, eff), pa in zip(stats, p_adj)]
