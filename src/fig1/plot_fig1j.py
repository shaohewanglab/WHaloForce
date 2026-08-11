"""Figure 1j: Kaplan-Meier survival of the quenched-state dwell time per force.

Empirical survival S(t)=P(T>t) of the quenched state for each force condition,
on a log x-axis. QU events cause the step drops; END/BL observations are
right-censored (no drop, marked with a tick). Under-sampled 3, 10, and 15 pN
conditions are dropped, matching plot_QU_force_dependence.py defaults. Reads
data/Fig1/Fig1j.csv and writes figures/Fig1j.svg and .png.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR  # noqa: E402

apply_style()

EXCLUDE = {3, 10, 15}
# Okabe-Ito palette. 8 pN is the reddish-purple #CC79A7, sampled from the
# published figure -- paler than the #CC0066 used previously.
COLORS = {3: "#E6AB02", 4: "#000000", 5: "#0072B2",
          6: "#D55E00", 7: "#009E73", 8: "#CC79A7"}


def km_estimate(times, observed):
    times = np.asarray(times, dtype=float)
    observed = np.asarray(observed, dtype=bool)
    event_times = np.unique(times[observed])
    surv_after = np.empty(event_times.size)
    s = 1.0
    for i, ut in enumerate(event_times):
        at_risk = int(np.sum(times >= ut))
        deaths = int(np.sum((times == ut) & observed))
        if at_risk > 0:
            s *= 1.0 - deaths / at_risk
        surv_after[i] = s
    return event_times, surv_after


def _surv_at(event_times, surv_after, t):
    idx = int(np.searchsorted(event_times, t, side="right")) - 1
    return float(surv_after[idx]) if idx >= 0 else 1.0


def main():
    df = pd.read_csv(DATA_DIR / "Fig1" / "Fig1j.csv")
    df = df[~df["force_pN"].isin(EXCLUDE)]
    forces = sorted(df["force_pN"].unique())

    all_times = df["dwell_time_s"].to_numpy(float)
    t_max = float(np.nanmax(all_times))
    pos = all_times[all_times > 0]
    t_min_pos = float(np.nanmin(pos)) if pos.size else 1.0
    left, right = t_min_pos * 0.6, t_max * 1.4

    fig, ax = plt.subplots(figsize=(1.8, 2.0))
    for f in forces:
        sub = df[df["force_pN"] == f]
        times = sub["dwell_time_s"].to_numpy(float)
        observed = sub["observed"].to_numpy(bool)
        event_times, surv_after = km_estimate(times, observed)
        color = COLORS[f]
        x_step = np.concatenate(([left], event_times, [right]))
        if surv_after.size:
            y_step = np.concatenate(([1.0], surv_after, [surv_after[-1]]))
        else:
            y_step = np.ones_like(x_step)
        ax.step(x_step, y_step, where="post", color=color, lw=1.0, label=f"{f} pN")
        if event_times.size:
            ax.scatter(event_times, surv_after, s=7, color=color,
                       edgecolors="none", alpha=0.9, zorder=3)
        cens = np.sort(times[~observed])
        if cens.size:
            y_c = [_surv_at(event_times, surv_after, tc) for tc in cens]
            ax.scatter(cens, y_c, s=18, marker="|", color=color,
                       linewidths=0.8, zorder=4)

    ax.set_xscale("log")
    ax.set_xlabel("Quenched state dwell time (s)")
    ax.set_ylabel("1 - CDF  (survival probability)")
    ax.set_xlim(left, right)
    ax.set_ylim(0, 1.05)
    ax.grid(True, linestyle="-", linewidth=0.3, color="0.85", alpha=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower left", ncol=1, frameon=False, fontsize=6)

    fig.tight_layout(pad=0.4)
    save(fig, "Fig1j")


if __name__ == "__main__":
    main()
