"""Extended Data Fig. 2a: normalized absorption (solid) and emission (dashed) spectra of
cpHalo-JFX554-HTL, cpWHalo-JFX554-HTL, and JFX554-HTL.

Absorption is normalized to the global peak absorbance across the three
samples; emission is normalized to the global peak emission intensity — the
same normalization as the original run_analysis.py combined abs/em panel.
Reads data/ExtDataFig2/ExtDataFig2a.csv and writes figures/ExtDataFig2a.svg and .png.
"""

import sys
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR  # noqa: E402

apply_style()

# Axes size measured from the published panel; margins are absolute so the
# canvas can be swapped in directly.
AX_W, AX_H = 5.388, 2.015
LEFT_IN, RIGHT_IN, TOP_IN, BOT_IN = 0.36, 0.09, 0.05, 0.32
FIG_W, FIG_H = AX_W + LEFT_IN + RIGHT_IN, AX_H + TOP_IN + BOT_IN

# Sampled from the published Extended Data Fig. 2a: cpHalo is black, cpWHalo blue, free dye
# orange. The earlier mapping was shifted by one (and used a green that does not
# appear in the figure at all).
COLORS = {
    "cpHalo-JFX554-HTL": "black",
    "cpWHalo-JFX554-HTL": "#1f77b4",
    "JFX554-HTL": "#ff7f0e",
}
ORDER = ["cpHalo-JFX554-HTL", "cpWHalo-JFX554-HTL", "JFX554-HTL"]
STYLE = {"absorption": "-", "emission": "--"}


def main():
    df = pd.read_csv(DATA_DIR / "ExtDataFig2" / "ExtDataFig2a.csv")

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    for sample in ORDER:
        for stype in ("absorption", "emission"):
            sub = df[(df["sample"] == sample) & (df["spectrum_type"] == stype)]
            sub = sub.sort_values("wavelength_nm")
            ax.plot(sub["wavelength_nm"], sub["normalized_value"],
                    color=COLORS[sample], ls=STYLE[stype], lw=1.0)

    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Normalized Abs / Em")
    ax.set_xlim(450, 700)
    ax.set_ylim(0, 1.08)
    ax.spines[["top", "right"]].set_visible(False)

    handles = [Line2D([0], [0], color=COLORS[s], lw=1.2, label=s) for s in ORDER]
    handles += [
        Line2D([0], [0], color="black", ls="-", lw=1.0, label="Absorption"),
        Line2D([0], [0], color="black", ls="--", lw=1.0, label="Emission"),
    ]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=6,
              handletextpad=0.5, labelspacing=0.25)

    for sp in ax.spines.values():      # framed, as published
        sp.set_visible(True)
    ax.tick_params(axis="both", labelsize=7, pad=1.5)
    fig.subplots_adjust(left=LEFT_IN / FIG_W, right=1 - RIGHT_IN / FIG_W,
                        top=1 - TOP_IN / FIG_H, bottom=BOT_IN / FIG_H)
    save(fig, "ExtDataFig2a")


if __name__ == "__main__":
    main()
