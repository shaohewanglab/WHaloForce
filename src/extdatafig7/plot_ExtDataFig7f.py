"""Extended Data Figure 7f: staged labeling of ectopically expressed Vin-WF in NIH-3T3 cells, JFX554-HTL.
Per-FOV in-FA lifetime (pooled-NNLS) after 12 h DMSO or 20 uM blebbistatin,
labeled for 1 h at the end of treatment (treat->label) or before it
(label->treat). Three independent experiments (2026-10-13, 2026-10-17, 2026-10-24)
pooled, one symbol. One FOV (3T3-185-stain-treat-DMSO-8, 2026-10-24) was
excluded upstream and is not in the CSV.

Reads data/ExtDataFig7/ExtDataFig7f.csv and writes figures/ExtDataFig7f.svg
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
    staged.plot_panel(staged.read(DATA_DIR / "ExtDataFig7" / "ExtDataFig7f.csv"),
                      "ExtDataFig7f")


if __name__ == "__main__":
    main()
