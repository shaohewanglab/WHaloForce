"""Figure 2j: normalized bright-FA lifetime vs beating (motion) power spectra for
representative Vin-WF (LV185-019) and Vin-cpH (LV186-003) FOVs. Two columns (one
per construct); within each column the lifetime power is stacked on top and the
beating power on the bottom, sharing the frequency axis (mirrors the original
fig_brightFA_spectra.py layout). The beat fundamental f0 is drawn as a dashed
line and the surrogate coupling r / p are annotated on the lifetime panel.

Reads data/Fig2/Fig2j.csv (+ Fig2j_meta.csv) and writes figures/Fig2j.svg/.png.
"""

import sys
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import apply_style, save, DATA_DIR  # noqa: E402

apply_style()

COLOR = {"Vin-WF": "#4C72B0", "Vin-cpH": "#DD8452"}
MOT_C = "0.45"
ORDER = ["Vin-WF", "Vin-cpH"]
FMAX = 2.5


def fmt_p(p):
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3f}" if p < 0.01 else f"p = {p:.2f}"


def main():
    df = pd.read_csv(DATA_DIR / "Fig2" / "Fig2j.csv")
    meta = pd.read_csv(DATA_DIR / "Fig2" / "Fig2j_meta.csv").set_index("construct")

    fig, axes = plt.subplots(2, 2, figsize=(3.6, 1.3), sharex=True,
                             gridspec_kw=dict(hspace=0.25, wspace=0.35))

    for j, construct in enumerate(ORDER):
        color = COLOR[construct]
        d = df[df["construct"] == construct]
        f = d["freq_hz"].to_numpy()
        sel = (f > 0) & (f <= FMAX)
        fv = f[sel]
        Pt = d["P_tau_brightFA_au"].to_numpy()[sel]
        Pm = d["P_motion_au"].to_numpy()[sel]
        f0 = float(meta.loc[construct, "f0_hz"])
        r = float(meta.loc[construct, "coupling_r"])
        p = float(meta.loc[construct, "coupling_p"])

        axl, axb = axes[0, j], axes[1, j]      # lifetime top, beating bottom
        axl.plot(fv, Pt / (Pt.max() + 1e-12), color=color, lw=0.9)
        axb.plot(fv, Pm / (Pm.max() + 1e-12), color=MOT_C, lw=0.9)
        for ax in (axl, axb):
            ax.axvline(f0, color="0.7", ls="--", lw=0.5, zorder=0)
            ax.set_ylim(0, 1.32)
            ax.set_yticks([0, 1])
            ax.set_xlim(0, FMAX)
            ax.set_xticks([0, 1, 2])
            ax.tick_params(labelsize=6.5)
            ax.spines[["top", "right"]].set_visible(False)
        axl.tick_params(labelbottom=False)
        axb.tick_params(axis="x", pad=1.5)
        axb.set_xlabel("Frequency (Hz)", fontsize=7, labelpad=1.5)

        axl.set_title(construct, loc="center", fontsize=7, color=color, pad=2)
        axl.text(0.025, 0.97, "lifetime", transform=axl.transAxes, va="top",
                 ha="left", fontsize=6, color=color)
        axb.text(0.025, 0.97, "beating", transform=axb.transAxes, va="top",
                 ha="left", fontsize=6, color=MOT_C)
        axl.text(0.975, 0.97, f"r = {r:.2f}, {fmt_p(p)}", transform=axl.transAxes,
                 va="top", ha="right", fontsize=6)

    fig.subplots_adjust(left=0.12, right=0.97, top=0.88, bottom=0.24)
    fig.text(0.03, 0.58, "Norm. power", rotation=90, va="center", ha="center",
             fontsize=7)
    save(fig, "Fig2j")


if __name__ == "__main__":
    main()
