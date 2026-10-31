"""
Extended Data Fig. 6d: brood size (eggs laid 24-48 h post-L4) per strain, mean +/- SEM,
Holm-Bonferroni-adjusted two-sided Welch's t-test vs N2.

NOTE ON STATISTICS: this panel keeps Welch's t-test and is deliberately not
unified with the Mann-Whitney U tests used in Fig. 6, for consistency with its
companion panel Extended Data Fig. 6e. See welch_holm_vs_control() in _common.py for the
rationale.

Reads data/ExtDataFig6/ExtDataFig6d.csv. Run:
    python plot_ExtDataFig6d.py
Writes figures/ExtDataFig6d.svg and .png.
"""
import pandas as pd

from _common import welch_holm_vs_control, bar_with_sem_panel, DATA_DIR

VALUE_COL = "brood_size"


def main():
    df = pd.read_csv(DATA_DIR / "ExtDataFig6" / "ExtDataFig6d.csv")
    p_adj = welch_holm_vs_control(df, VALUE_COL)
    for g, p in p_adj.items():
        print(f"{g}: Holm-adjusted p = {p:.4g}")
    bar_with_sem_panel(
        df, VALUE_COL,
        ylabel="Brood size (eggs/plate)",
        ymax=180, tick_step=20,
        out_stub="ExtDataFig6d",
        p_adj=p_adj,
        ax_w=1.18, ax_h=1.107,
    )


if __name__ == "__main__":
    main()
