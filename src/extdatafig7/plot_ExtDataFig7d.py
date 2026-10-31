"""Extended Data Figure 7d: staged labeling of endogenous Vcl-WF in homozygous knock-in MEFs, JF571-HTL.
Per-FOV in-FA lifetime (pooled-NNLS) after 12 h DMSO or 20 uM blebbistatin,
labeled for 1 h at the end of treatment (treat->label) or before it
(label->treat). Two independent experiments (2026-10-24, 2026-10-27) pooled, one symbol.

Reads data/ExtDataFig7/ExtDataFig7d.csv and writes figures/ExtDataFig7d.svg
and .png. Drawing and statistics (two-sided Mann-Whitney U; Holm over the two
Bleb vs DMSO comparisons, and separately over the two schedule comparisons) are
in src/common/staged.py.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, DATA_DIR  # noqa: E402
import staged  # noqa: E402

apply_style()


def main():
    staged.plot_panel(staged.read(DATA_DIR / "ExtDataFig7" / "ExtDataFig7d.csv"),
                      "ExtDataFig7d")


if __name__ == "__main__":
    main()
