"""Shared figure style and paths for the WHaloForce plotting scripts.

House style: Arial 7 pt, editable text in SVG, 0.8 pt axes, hidden top and
right spines, no in-figure legend by default. See docs and FIGURE_STYLE notes.
"""

from pathlib import Path

from matplotlib import rcParams

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
FIG_DIR = REPO_ROOT / "figures"

# --------------------------------------------------------------------------- #
# Colors used in the published figures.
#
# These are transcribed from the scripts that produced the final panels, so the
# regenerated panels match the manuscript rather than falling back to a
# monochrome or matplotlib-default scheme. The source script is named for each.
# --------------------------------------------------------------------------- #

# Vinculin constructs (3T3 / MEF panels). Source: 3T3-vinculin
# plot_dmso_bleb20_box.py and 3T3-distribution fa_replot_panels.py.
VIN_COLORS = {
    "Vin-WF": "#4C72B0",     # blue
    "Vin-cpH": "#DD8452",    # orange
    "Vin-N-WF": "#C44E52",   # red
    "Vin-WF-TL": "#55A868",  # green (construct 187, excluded from all panels)
}

# Cadherin-catenin constructs (MDCKII junction panels). Source:
# MDCKII-extrusion plot_radial_profiles.py GROUP_COLORS_LIGHT.
JUNCTION_COLORS = {"WF": "#2c7fb8", "cpH": "#e6842e", "C-WF": "#d62728"}

# NOTE: several panels are deliberately monochrome in the published figures and
# must stay that way -- Fig. 5c (gray boxes and dots), Fig. 6d/6f/6h, Fig. 4e-g
# data, and Fig. S6d/S6e (gray bars). Colored palettes exist in the working-folder
# scripts for some of these (a Dark2 by-dataset scheme for the MEF inhibitor
# panel, per-strain hues for the worm panels), but those were exploratory and did
# NOT make it into the figures. Verified against the published PDF; do not
# "restore" color to those panels.

# Category swatches drawn under the x-axis of Fig. 4e-g, and the matching
# schematic labels. Sampled from the published figure: the swatches are these
# Okabe-Ito hues at alpha 0.7 over white. Perpendicular pairs with pharynx, and
# parallel with body-wall muscle, so one pair of colours serves all three panels.
CATEGORY_COLORS = {"perpendicular": "#E69F00", "parallel": "#009E73",
                   "pharynx": "#E69F00", "body_wall_muscle": "#009E73",
                   "bodywall": "#009E73"}
CATEGORY_SWATCH_ALPHA = 0.7


def add_category_swatches(ax, keys, width=0.44, height=0.035, pad=0.135):
    """Draw the small colour swatches under the x tick labels of Fig. 4e-g.

    The published panels mark each category with a short colour bar below the
    axis (gold = perpendicular / pharynx, teal-green = parallel / body-wall
    muscle) so the plot keys back to the schematic above it. The data themselves
    stay monochrome, as in the figure.
    """
    from matplotlib import pyplot as plt
    from matplotlib.transforms import blended_transform_factory
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    for xi, key in enumerate(keys):
        ax.add_patch(plt.Rectangle(
            (xi - width / 2, -pad), width, height, transform=tr, clip_on=False,
            facecolor=CATEGORY_COLORS[key], alpha=CATEGORY_SWATCH_ALPHA,
            edgecolor="none", zorder=5))

# Alpha conventions for the colored box+dot panels (fa_replot_panels.render and
# plot_dmso_bleb20_box.py): translucent boxes, more opaque dots, white dot edge.
BOX_ALPHA = 0.30
DOT_ALPHA = 0.90
DOT_EDGE = "white"


def apply_style():
    rcParams.update({
        "svg.fonttype": "none",     # keep text editable in SVG
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        "font.size": 7,
        "axes.titlesize": 7,
        "axes.labelsize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        # Math text (the tau symbols, and the 10^n in p-value labels) must match
        # the surrounding Arial. The default "stixsans" set renders them in STIX
        # Sans, which reads as a different font next to plain Arial text -- most
        # visibly in "p = 3.3x10^-3", where "p =" was Arial and "10^-3" was not.
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "mathtext.default": "rm",
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.8,
        "ytick.minor.width": 0.8,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
    })


# --------------------------------------------------------------------------- #
# Text that stays editable in the SVG.
#
# Matplotlib renders any $...$ mathtext by positioning each glyph individually,
# so it lands in the SVG as one <text> holding a dozen <tspan>s with explicit
# coordinates. Illustrator imports that as a scatter of separate characters
# rather than an editable string -- the "fragmented" text. Plain strings, by
# contrast, come out as a single <text> you can click and retype.
#
# Convention:
#   * Axis labels and p-values: NO mathtext. Use the literal Greek character
#     (Arial has tau and Delta), and build p-value exponents from two plain text
#     runs via p_text(), so both stay editable.
#   * Log-axis tick labels: leave matplotlib's default mathtext alone. Hand-built
#     replacements throw off the axis-label spacing for little gain, and the
#     decade labels are rarely edited.
#
# Unicode superscripts are NOT an option for the p-values -- macOS Arial only
# carries the superscripts for 1, 2 and 3, and has no superscript minus, so the
# rest would silently fall back to another font.
# --------------------------------------------------------------------------- #

TAU = "τ"        # Arial has these; "$\\tau$" would fragment the label
DELTA = "Δ"


def format_p(p, spell_out_above=0.01):
    """p-value label, split into (base, exponent) plain-text runs.

    Returns (text, None) when the value is spelled out, or ("p = 1.8x10", "-8")
    when it needs scientific notation. Values above `spell_out_above` are spelled
    out; below it they are given to one decimal as m x 10^e.
    """
    import numpy as np
    if not np.isfinite(p):
        return "n/a", None
    if p <= 0:                      # underflow to exactly 0
        return "p < 1×10", "-300"
    if p > spell_out_above:
        return f"p = {p:.3g}", None
    exp = int(np.floor(np.log10(p)))
    return f"p = {p / 10 ** exp:.1f}×10", f"{exp}"


# Gap in points between a significance bar and the p-value label above it.
P_LABEL_PAD = 2.5

# How far the exponent's baseline sits above the base text's, as a fraction of
# the base font size. Typographic superscripts sit around 0.35-0.45 em; matching
# the tops of the two boxes (the obvious approach) leaves them looking too low.
SUP_RISE = 0.45


def sup_text(ax, x, y, base, exp, fontsize=7, ha="center", va="bottom",
             zorder=6, xycoords="data", pad_pt=0.0, stacked=False, **kw):
    """Text with a superscript, kept editable in the SVG.

    Draws two plain text runs -- e.g. "p = 1.8x10" and a smaller, raised "-8" --
    packed side by side. Both stay single <text> elements in Arial, unlike
    mathtext, which matplotlib positions glyph by glyph.

    The exponent is raised by padding *below* it and bottom-aligning the pair,
    which gives direct control over the rise (SUP_RISE, in em of the base font).
    """
    from matplotlib.offsetbox import (AnnotationBbox, DrawingArea, HPacker,
                                      TextArea, VPacker)
    tp = dict(fontsize=fontsize, **kw)
    rise = SUP_RISE * fontsize                      # points
    stack = VPacker(children=[TextArea(exp, textprops=dict(tp, fontsize=fontsize * 0.72)),
                              DrawingArea(0.0, rise, 0.0, 0.0)],
                    align="left", pad=0.0, sep=0.0)
    if stacked:
        # Narrow panels cannot take "p = 1.8x10^-8" on one line, so break after
        # "p =" and put the value underneath, as the published Fig. 4 panels do.
        prefix, value = base.split("=", 1)
        row = HPacker(children=[TextArea(value.strip(), textprops=tp), stack],
                      align="bottom", pad=0.0, sep=0.0)
        box = VPacker(children=[TextArea(prefix.strip() + "=", textprops=tp), row],
                      align="left", pad=0.0, sep=1.0)
    else:
        box = HPacker(children=[TextArea(base, textprops=tp), stack],
                      align="bottom", pad=0.0, sep=0.0)
    ab = AnnotationBbox(
        box, (x, y), xybox=(0.0, pad_pt), xycoords=xycoords,
        boxcoords="offset points", frameon=False, pad=0.0,
        box_alignment=({"center": 0.5, "left": 0.0, "right": 1.0}[ha],
                       {"bottom": 0.0, "center": 0.5, "top": 1.0}[va]),
        annotation_clip=False, zorder=zorder)
    ax.add_artist(ab)
    return ab


def p_text(ax, x, y, p, fontsize=7, ha="center", va="bottom", zorder=6,
           spell_out_above=0.01, pad_pt=P_LABEL_PAD, stacked=False, **kw):
    """Draw a p-value label that stays editable in the SVG.

    `pad_pt` lifts the label off the significance bar so the two do not touch.
    `stacked` breaks the label after "p =" for narrow panels.
    """
    base, exp = format_p(p, spell_out_above)
    if exp is None:
        if stacked:
            prefix, value = base.split("=", 1)
            base = f"{prefix.strip()}=\n{value.strip()}"
        return ax.annotate(base, (x, y), xytext=(0.0, pad_pt),
                           textcoords="offset points", ha=ha, va=va,
                           fontsize=fontsize, zorder=zorder,
                           annotation_clip=False, **kw)
    return sup_text(ax, x, y, base, exp, fontsize=fontsize, ha=ha, va=va,
                    zorder=zorder, pad_pt=pad_pt, stacked=stacked, **kw)


# --------------------------------------------------------------------------- #
# Dot styling for the per-observation scatter panels (Fig. 6d/6f/6h).
# --------------------------------------------------------------------------- #

DOT_SIZE = 18          # scatter `s`, i.e. marker area in pt^2
DOT_COLOR = "0.3"
DOT_ALPHA = 0.55


def dot_diameter_in(size=DOT_SIZE):
    """Marker diameter in inches for a scatter `s` (which is an area in pt^2)."""
    return (size ** 0.5) / 72.0


def swarm_offsets(values, dx, dy, max_off, step_frac=0.25):
    """Beeswarm x offsets: nudge overlapping dots sideways instead of jittering.

    Random jitter still lets dots land on top of each other, which is very
    visible once the markers are large. This walks the points in y order and
    gives each the smallest sideways offset that clears every already-placed
    point within one marker height, so a dense cluster fans out into columns
    while an isolated point stays centred.

    `dx`/`dy` are the marker diameter expressed in x/y data units, and `max_off`
    caps the spread so neighbouring groups never run into each other. Points
    that cannot be placed within the cap fall back to centred.
    """
    import numpy as np
    v = np.asarray(values, dtype=float)
    step = dx * step_frac
    n_steps = int(np.ceil(max_off / step)) if step > 0 else 0
    offs = np.zeros(v.size)
    placed = []
    for i in np.argsort(v, kind="stable"):
        near = [o for o, y in placed if abs(y - v[i]) < dy]
        chosen = 0.0
        for k in range(n_steps + 1):
            cands = (0.0,) if k == 0 else (k * step, -k * step)
            hit = next((c for c in cands
                        if all(abs(c - o) >= dx for o in near)), None)
            if hit is not None:
                chosen = hit
                break
        offs[i] = chosen
        placed.append((chosen, v[i]))
    return offs


def save(fig, stem):
    """Save a figure as SVG plus a 300-dpi PNG preview into figures/."""
    FIG_DIR.mkdir(exist_ok=True)
    fig.savefig(FIG_DIR / f"{stem}.svg")
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=300)
    print(f"wrote figures/{stem}.svg and .png")
