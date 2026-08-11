"""Figure 5c: per-FOV in-FA pooled-NNLS fluorescence lifetime under cytoskeletal
inhibitors in Vinculin-WHaloForce (Vcl-WF) MEFs.

Each FOV contributes one pooled-decay NNLS lifetime (tau_inFA_nnls_ns): every
in-FA pixel of the image is summed into a single high-SNR decay and fit once.
Treatments are pooled across three acquisitions (Y5uM-12h, drug-panel,
drug-panel2); the CK666 arms of drug-panel2 are excluded. Each inhibitor is
compared with pooled DMSO by a two-sided Mann-Whitney U test, Holm-Bonferroni
adjusted across the six inhibitor comparisons.

Reads data/Fig5/Fig5c.csv and writes figures/Fig5c.svg and .png. Exact
Holm-adjusted p values are printed.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from matplotlib import pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, p_text, DATA_DIR, TAU  # noqa: E402

apply_style()

# category structure: (label, [group codes]); DMSO is the pooled vehicle.
CATEGORIES = [
    ("vehicle", ["DMSO"]),
    ("ROCK inhib. (Y-27632)", ["Y5", "Y20"]),
    ("myosin-II inhib. (blebbistatin)", ["Bleb5", "Bleb20"]),
    ("actin depoly. (latrunculin A)", ["LatA0.5", "LatA2"]),
]
ORDER = [g for _, gs in CATEGORIES for g in gs]
INHIBITORS = [g for g in ORDER if g != "DMSO"]

XLABEL = {
    "DMSO": "DMSO", "Y5": "Y 5 µM", "Y20": "Y 20 µM",
    "Bleb5": "Bleb 5 µM", "Bleb20": "Bleb 20 µM",
    "LatA0.5": "LatA 0.5 µM", "LatA2": "LatA 2 µM",
}

BOX_FILL = "0.88"   # light fill so the dots read on top of it
BOX_LINE = "black"  # box, whiskers, caps and median all black
DOT_GRAY = "0.35"   # dots darker than the fill, lighter than the box line
ANNOT_GRAY = "0.25" # p values and n labels
CAT_GAP = 0.45


def holm(pvals):
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    adj = np.empty(p.size)
    run = 0.0
    for i, idx in enumerate(order):
        run = max(run, (p.size - i) * p[idx])
        adj[idx] = min(run, 1.0)
    return adj


def positions():
    pos, x = {}, 1.0
    for ci, (_, groups) in enumerate(CATEGORIES):
        if ci:
            x += CAT_GAP
        for g in groups:
            pos[g] = x
            x += 1.0
    return pos


def main():
    df = pd.read_csv(DATA_DIR / "Fig5" / "Fig5c.csv")
    df = df[df.treatment.isin(ORDER)].copy()

    data = {g: df.loc[df.treatment == g, "tau_inFA_ns"].dropna().to_numpy()
            for g in ORDER}
    ref = data["DMSO"]

    raw = {}
    for g in INHIBITORS:
        x = data[g]
        p = mannwhitneyu(x, ref, alternative="two-sided").pvalue
        raw[g] = dict(p=p, n=len(x), median=float(np.median(x)))
    padj = holm([raw[g]["p"] for g in INHIBITORS])
    for g, q in zip(INHIBITORS, padj):
        raw[g]["p_holm"] = float(q)
    dmso_med = float(np.median(ref))

    print("=== Fig. 5c: per-FOV in-FA pooled-NNLS tau, each inhibitor vs pooled "
          "DMSO (Mann-Whitney two-sided, Holm-Bonferroni) ===")
    print("%-9s %4s %8s %10s %10s" % ("group", "n", "median", "p", "p_holm"))
    print("%-9s %4d %8.3f %10s %10s" % ("DMSO", len(ref), dmso_med, "-", "-"))
    for g in INHIBITORS:
        r = raw[g]
        print("%-9s %4d %8.3f %10.3g %10.3g" %
              (g, r["n"], r["median"], r["p"], r["p_holm"]))

    pos = positions()
    fig, ax = plt.subplots(figsize=(3.0, 2.2))
    xs = [pos[g] for g in ORDER]
    boxdata = [data[g] for g in ORDER]

    ax.boxplot(
        boxdata, positions=xs, widths=0.55, showfliers=False, patch_artist=True,
        medianprops=dict(color=BOX_LINE, lw=1.1),
        boxprops=dict(facecolor=BOX_FILL, edgecolor=BOX_LINE, lw=0.8),
        whiskerprops=dict(color=BOX_LINE, lw=0.8),
        capprops=dict(color=BOX_LINE, lw=0.8))

    rng = np.random.default_rng(0)
    for g in ORDER:
        v = data[g]
        jx = pos[g] + (rng.random(len(v)) - 0.5) * 0.26
        ax.scatter(jx, v, s=6, color=DOT_GRAY, alpha=0.85,
                   linewidths=0.3, edgecolor="white", zorder=3)

    ax.axhline(dmso_med, color="0.6", lw=0.6, ls=(0, (4, 3)), zorder=0)

    ax.set_ylim(2.05, 3.72)
    ax.set_yticks([2.0, 2.5, 3.0])
    ax.set_ylabel(f"Per-FOV in-FA {TAU} (ns)")
    ax.set_xlim(min(xs) - 0.8, max(xs) + 0.8)
    ax.spines[["top", "right"]].set_visible(False)

    p_hi, p_lo = 3.38, 3.08
    for i, g in enumerate(INHIBITORS):
        y = p_hi if i % 2 == 0 else p_lo
        p_text(ax, pos[g], y, raw[g]["p_holm"], fontsize=7, ha="center",
               va="bottom", pad_pt=0.0, stacked=True, color=ANNOT_GRAY)

    ax.set_xticks(xs)
    ax.set_xticklabels([XLABEL[g] for g in ORDER], rotation=35, ha="right")
    ax.tick_params(axis="x", pad=20)

    trans = ax.get_xaxis_transform()
    for g in ORDER:
        ax.text(pos[g], -0.075, f"{len(data[g])}", ha="center", va="top",
                fontsize=7, color=ANNOT_GRAY, transform=trans)
    ax.text(min(xs) - 0.78, -0.075, "n:", ha="left", va="top", fontsize=7,
            color=ANNOT_GRAY, transform=trans)

    fig.tight_layout(pad=0.3)
    save(fig, "Fig5c")


if __name__ == "__main__":
    main()
