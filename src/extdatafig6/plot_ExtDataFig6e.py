"""
Extended Data Fig. 6e: embryonic viability (%) of embryos laid 24-48 h post-L4 per
strain, mean +/- SD, Holm-Bonferroni-adjusted two-sided Welch's t-test vs
N2.

viability_percent per worm = 100 * n_hatched / n_embryos, where n_embryos =
hatched + unhatched + UFO (unfertilized) eggs on that worm's plate.

NOTE ON STATISTICS: this panel keeps Welch's t-test and is deliberately not
unified with the Mann-Whitney U tests used in Fig. 6. Viability is a
ceiling-bounded percentage with most N2 worms at exactly 100%, where a rank
test reports strong significance for ~1 percentage-point shifts. See
welch_holm_vs_control() in _common.py for the full rationale and the numbers
both tests give.

Reads data/ExtDataFig6/ExtDataFig6e.csv. Run:
    python plot_ExtDataFig6e.py
Writes figures/ExtDataFig6e.svg and .png.
"""
import pandas as pd

from _common import welch_holm_vs_control, bar_with_sd_panel, DATA_DIR

VALUE_COL = "viability_percent"


def main():
    df = pd.read_csv(DATA_DIR / "ExtDataFig6" / "ExtDataFig6e.csv")
    p_adj = welch_holm_vs_control(df, VALUE_COL)
    for g, p in p_adj.items():
        print(f"{g}: Holm-adjusted p = {p:.4g}")
    bar_with_sd_panel(
        df, VALUE_COL,
        ylabel="Embryonic viability (%)",
        ymax=100, tick_step=20,
        out_stub="ExtDataFig6e",
        p_adj=p_adj,
        ax_w=1.75, ax_h=1.11,
    )


if __name__ == "__main__":
    main()
