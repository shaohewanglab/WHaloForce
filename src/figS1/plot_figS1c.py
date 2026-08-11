"""Figure S1c: intensity-weighted mean fluorescence lifetime of HaloTag vs
cpHalo vs cpWHalo across 17 JF dyes in live HEK293 cells.

Per-FOV intensity-weighted mean tau (ns) from a three-component fit of the dye
channel, points (jittered) with mean +/- SD, faceted per dye (ordered by dye
excitation maximum). Reads data/FigS1/FigS1c.csv and writes figures/FigS1c.svg
and .png.
"""

import sys
from pathlib import Path

# reuse the faceted panel from the S1b script
sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_figS1b import facet_panel  # noqa: E402
from style import TAU  # noqa: E402


def main():
    facet_panel(  # S1c subplots are slightly wider and more tightly stacked
        "FigS1c.csv", "tau_int_ns",
                f"Mean {TAU}, intensity weighted (ns)", "FigS1c", logy=False,
        ax_w=0.667, ax_h=0.983, wgap=0.328, hgap=0.235)


if __name__ == "__main__":
    main()
