"""Shared drawing for the staged-labeling panels (Fig. 5d, Extended Data Fig.
7c,d,f-h): per-FOV in-FA lifetime under DMSO or 20 uM blebbistatin (Bleb) for
12 h, with the dye-HTL labeling (1 h) either at the end of treatment
(treat->label) or before it (label->treat).

Each panel's CSV has one row per FOV (different FOVs contain different cells),
pooled across the independent experiments in its `experiment_date` column and
drawn with a single symbol. Statistics: two-sided Mann-Whitney U on the pooled
FOVs; Holm-Bonferroni over the comparisons drawn in a panel, with the two Bleb
vs DMSO comparisons (one per schedule) as one family and the two schedule
comparisons (same treatment, treat->label vs label->treat) as a separate family.

Style is the Fig. 2e house style (box + jittered dots, vehicle lighter than
drug), with label->treat drawn in the same hue mixed 35% toward black. Fixed
figure size (default 2.0 x 2.2 in); the y range is tight per panel (lowest FOV
- 0.08 ns to the top of the highest p label), so the ns-per-pt scale differs
between panels.
"""

import matplotlib.colors
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import mannwhitneyu

from style import p_text, TAU, VIN_COLORS, DOT_EDGE, save

SCHEDULES = ["treat->label", "label->treat"]
SCHEDULE_LABEL = {"treat->label": "treat→label", "label->treat": "label→treat"}
TREATMENTS = ["DMSO", "Bleb"]
DARKEN = {"treat->label": 0.0, "label->treat": 0.35}
BOX_ALPHA = {"DMSO": 0.25, "Bleb": 0.45}
DOT_ALPHA = {"DMSO": 0.60, "Bleb": 0.90}

# geometry in pt
LEFT, RIGHT, TOP, BOTTOM = 28.3, 3.0, 4.0, 33.0
LABEL_PT = 11.0          # p label height incl. pad + exponent
PAD_LO = 0.08            # ns below the lowest FOV


def holm(pvals):
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    adj = np.empty(p.size)
    run = 0.0
    for i, idx in enumerate(order):
        run = max(run, (p.size - i) * p[idx])
        adj[idx] = min(run, 1.0)
    return adj


def shade(hue, schedule):
    f = DARKEN[schedule]
    return tuple(c * (1 - f) for c in matplotlib.colors.to_rgb(hue))


def groups(df):
    """{(schedule, treatment): array of per-FOV tau}, in CSV row order."""
    return {(s, t): df.loc[(df.schedule == s) & (df.treatment == t),
                           "tau_inFA_ns"].to_numpy()
            for s in SCHEDULES for t in TREATMENTS}


def stats(data):
    """Bleb vs DMSO per schedule, and schedule per treatment; Holm within each
    family of two. Returns {(kind, key): (p_raw, p_holm)}."""
    out = {}
    fams = {"drug": [(s, (data[(s, "DMSO")], data[(s, "Bleb")])) for s in SCHEDULES],
            "schedule": [(t, (data[(SCHEDULES[0], t)], data[(SCHEDULES[1], t)]))
                         for t in TREATMENTS]}
    for kind, fam in fams.items():
        praw = [mannwhitneyu(a, b, alternative="two-sided").pvalue for _, (a, b) in fam]
        for (key, _), pr, pa in zip(fam, praw, holm(praw)):
            out[(kind, key)] = (pr, pa)
    return out


def positions():
    return {(s, t): 1 + j + k * (len(TREATMENTS) + 0.6)
            for k, s in enumerate(SCHEDULES) for j, t in enumerate(TREATMENTS)}


def layout(data, pos, yl, ax_h):
    """[(x0, x1, y, kind, key)]: Bleb vs DMSO bars 0.11 span above the data they
    span, then the schedule bars stacked above them."""
    span = yl[1] - yl[0]
    step, lab = 0.11 * span, LABEL_PT / ax_h * span
    out = []
    for s in SCHEDULES:
        top = max(data[(s, t)].max() for t in TREATMENTS)
        out.append((pos[(s, "DMSO")], pos[(s, "Bleb")], top + step, "drug", s))
    clear = lab + 0.045 * span
    y = max([max(v.max() for v in data.values()) + step]
            + [o[2] + clear for o in out])
    for t in TREATMENTS:
        out.append((pos[(SCHEDULES[0], t)], pos[(SCHEDULES[1], t)], y, "schedule", t))
        y += clear
    return out


def plot_panel(df, stem, w=2.0, h=2.2, tag=""):
    data = groups(df)
    st = stats(data)
    for (kind, key), (pr, pa) in st.items():
        a, b = ((data[(key, "DMSO")], data[(key, "Bleb")]) if kind == "drug"
                else (data[(SCHEDULES[0], key)], data[(SCHEDULES[1], key)]))
        print(f"[{tag or stem}] {kind} {key}: n={a.size}/{b.size}, MWU p={pr:.3g}, "
              f"Holm p={pa:.3g}")

    W, H = w * 72, h * 72
    ax_h = H - TOP - BOTTOM
    pos = positions()
    # bracket spacing is a fraction of the span, so iterate the top to a fixed
    # point (top of the highest p label = top of the axes)
    lo = min(v.min() for v in data.values()) - PAD_LO
    yl = (lo, max(v.max() for v in data.values()) + 1.0)
    for _ in range(100):
        lay = layout(data, pos, yl, ax_h)
        yl = (lo, max(o[2] for o in lay) + LABEL_PT / ax_h * (yl[1] - yl[0]))

    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([LEFT / W, BOTTOM / H, (W - LEFT - RIGHT) / W, ax_h / H])
    hue = VIN_COLORS["Vin-WF"]      # one hue for Vin-WF and endogenous Vcl-WF
    keys = [(s, t) for s in SCHEDULES for t in TREATMENTS]
    xs = [pos[k] for k in keys]
    bp = ax.boxplot([data[k] for k in keys], positions=xs, widths=0.55,
                    patch_artist=True, showfliers=False, zorder=2,
                    medianprops=dict(color="black", linewidth=1.1),
                    whiskerprops=dict(color="black", linewidth=0.8),
                    capprops=dict(color="black", linewidth=0.8))
    for patch, (s, t) in zip(bp["boxes"], keys):
        c = shade(hue, s)
        patch.set_facecolor(c); patch.set_edgecolor(c); patch.set_linewidth(0.8)
        patch.set_alpha(BOX_ALPHA[t])
    rng = np.random.default_rng(0)
    for x, (s, t) in zip(xs, keys):
        v = data[(s, t)]
        jx = x + (rng.random(v.size) - 0.5) * 0.25
        ax.scatter(jx, v, s=8, color=shade(hue, s), alpha=DOT_ALPHA[t],
                   edgecolors=DOT_EDGE, linewidths=0.3, zorder=3, marker="o")
    step = 0.11 * (yl[1] - yl[0])
    for x0, x1, y, kind, key in lay:
        yt = y - 0.30 * step
        ax.plot([x0, x0, x1, x1], [yt, y, y, yt], lw=0.8, color="black",
                clip_on=False, zorder=5)
        p_text(ax, (x0 + x1) / 2, y, st[(kind, key)][1])
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{t}\n(n={data[(s, t)].size})" for s, t in keys], fontsize=7)
    for s in SCHEDULES:
        xc = (pos[(s, TREATMENTS[0])] + pos[(s, TREATMENTS[-1])]) / 2
        ax.annotate(SCHEDULE_LABEL[s], (xc, 0), xycoords=("data", "axes fraction"),
                    xytext=(0, -23), textcoords="offset points",
                    ha="center", va="top", fontsize=7)
    ax.set_xlim(0.4, max(xs) + 0.6)
    ax.set_ylim(*yl)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel(f"Per-FOV in-FA {TAU} (ns)")
    save(fig, stem)
    return st


def read(csv):
    df = pd.read_csv(csv)
    assert set(df.schedule) <= set(SCHEDULES) and set(df.treatment) <= set(TREATMENTS)
    return df
