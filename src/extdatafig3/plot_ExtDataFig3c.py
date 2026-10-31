"""Extended Data Fig. 3c: representative pooled in-FA photon-decay with the NNLS fit and
weighted residuals underneath (Vin-WF, one DMSO field of view).

Reads data/ExtDataFig3/ExtDataFig3c.csv (the exported decay curve) and ExtDataFig3c_meta.csv (the
fit summary for the title), and writes figures/ExtDataFig3c.svg and .png.

The decay curve was exported once from the raw photon data (see
docs/NON_SCRIPTABLE_STEPS.md); this script only renders the committed CSV.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR, FIG_DIR  # noqa: E402

apply_style()

# Axes size measured from the published panel; margins are absolute so the
# canvas can be swapped in directly.
# Two stacked axes: the decay and its weighted residual.
AX_W, AX_H0, AX_H1, GAP_IN = 5.383, 1.210, 0.565, 0.04
LEFT_IN, RIGHT_IN, TOP_IN, BOT_IN = 0.48, 0.04, 0.24, 0.30
FIG_W = AX_W + LEFT_IN + RIGHT_IN
FIG_H = AX_H0 + AX_H1 + GAP_IN + TOP_IN + BOT_IN


def main():
    d = pd.read_csv(DATA_DIR / "ExtDataFig3" / "ExtDataFig3c.csv")
    m = pd.read_csv(DATA_DIR / "ExtDataFig3" / "ExtDataFig3c_meta.csv").iloc[0]

    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(FIG_W, FIG_H), sharex=True,
        height_ratios=[AX_H0, AX_H1], gridspec_kw=dict(hspace=0.08))

    ax0.semilogy(d.time_ns, np.clip(d.counts, 1, None), ".", ms=1.8, color="0.35",
                 label="pooled in-FA decay", rasterized=True)
    ax0.semilogy(d.time_ns, np.clip(d.fit, 1, None), "-", lw=1.0, color="#C44E52",
                 label="NNLS fit")
    ax0.set_ylabel("Photon counts\n(log)")
    ax0.legend(loc="upper right", fontsize=6, frameon=False, handlelength=1.5)
    ax0.set_title(
        f"Pooled in-FA NNLS fit - {m.construct}, {m.image} ({m.treatment})   "
        rf"$\tau_{{int}}$={m.tau_int_ns:.2f} ns, {m.n_photons/1e3:.0f}k photons, "
        rf"$n_{{FA}}$={int(m.n_fa)}, reduced $\chi^2$={m.chi2_red:.0f}",
        fontsize=6.5)
    ax0.spines[["top", "right"]].set_visible(False)

    ax1.plot(d.time_ns, d.residual, "-", lw=0.5, color="#4C72B0")
    ax1.axhline(0, color="k", lw=0.5)
    ax1.set_ylabel("Weighted\nresidual")
    ax1.set_xlabel("time (ns)")
    ax1.set_xlim(d.time_ns.min(), d.time_ns.max())
    ax1.spines[["top", "right"]].set_visible(False)

    fig.align_ylabels([ax0, ax1])
    for a in (ax0, ax1):
        a.tick_params(axis="both", labelsize=7, pad=1.5)
    fig.subplots_adjust(left=LEFT_IN / FIG_W, right=1 - RIGHT_IN / FIG_W,
                        top=1 - TOP_IN / FIG_H, bottom=BOT_IN / FIG_H,
                        hspace=GAP_IN / ((AX_H0 + AX_H1) / 2))
    save(fig, "ExtDataFig3c")


if __name__ == "__main__":
    main()
