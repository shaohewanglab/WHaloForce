# WHaloForce reproducibility package

Analysis code and source data to regenerate the plotted figure panels of
Wu et al., 2026, WHaloForce (a FLIM-based molecular tension sensor).

Each plotted panel has one tidy source-data CSV under `data/` and one plotting
script under `src/`. The scripts read only the committed CSVs (relative paths)
and write panels to `figures/`. Raw microscopy data are large and live in an
external repository (see `docs/RAW_DATA.md`); they are not needed to regenerate
the plots.

## Requirements

The analysis reported in the paper was performed on Python 3.11; this package was
last verified on 3.12.11 (see below). Either works. Install the minimal plotting
dependencies:

```
pip install -r requirements.txt
```

or build the fuller environment (adds the upstream FLIM example steps):

```
conda env create -f environment.yml && conda activate whaloforce
```

## Regenerate all figures

```
python make_all_figures.py
```

This runs every `src/<fig>/plot_*.py` and writes each panel as SVG plus a
300-dpi PNG into `figures/`. To smoke-test the whole package -- every panel
regenerating from its committed CSV, plus the upstream analysis-step checks:

```
python examples/run_examples.py
```

Last run: 41/41 panels and 2/2 step checks passed on Python 3.12.11 with the
pinned versions in `requirements.txt`.

## Figure panel to script and source data

| Panel | Plotted quantity | Script | Source data |
|---|---|---|---|
| 1d | Mean lifetime (int/amp weighted), cpHalo/cpWHalo/dye | src/fig1/plot_fig1d.py | data/Fig1/Fig1d.csv |
| 1f | cpHalo-554 force-ramp trace | src/fig1/plot_fig1f.py | data/Fig1/Fig1f.csv |
| 1g | cpWHalo-554 force-ramp trace | src/fig1/plot_fig1g.py | data/Fig1/Fig1g.csv |
| 1h | cpWHalo-554 force-clamp trace | src/fig1/plot_fig1h.py | data/Fig1/Fig1h.csv |
| 1i | Quenched-state dwell times (box) | src/fig1/plot_fig1i.py | data/Fig1/Fig1i.csv |
| 1j | Dwell-time survival | src/fig1/plot_fig1j.py | data/Fig1/Fig1j.csv |
| 2c | Per-cell in-FA lifetime, untreated | src/fig2/plot_fig2c.py | data/Fig2/Fig2c.csv |
| 2d | Per-cell SD of per-FA lifetime | src/fig2/plot_fig2d.py | data/Fig2/Fig2d.csv |
| 2e | Per-cell in-FA lifetime, DMSO vs blebbistatin | src/fig2/plot_fig2e.py | data/Fig2/Fig2e.csv |
| 2g | Lifetime change vs time (spreading) | src/fig2/plot_fig2g.py | data/Fig2/Fig2g.csv |
| 2h | Vin-WF cardiomyocyte traces | src/fig2/plot_fig2h.py | data/Fig2/Fig2h.csv |
| 2i | Vin-cpH cardiomyocyte traces | src/fig2/plot_fig2i.py | data/Fig2/Fig2i.csv |
| 2j | Power spectra (beating, lifetime) | src/fig2/plot_fig2j.py | data/Fig2/Fig2j.csv, Fig2j_meta.csv |
| 2k | Lifetime-beating coupling r per FOV | src/fig2/plot_fig2k.py | data/Fig2/Fig2k.csv |
| 3c | E-cadherin per-FOV junction lifetime | src/fig3/plot_fig3c.py | data/Fig3/Fig3c.csv |
| 3d | E-cadherin radial lifetime profile | src/fig3/plot_fig3d.py | data/Fig3/Fig3d.csv |
| 3f | alpha-catenin per-FOV junction lifetime | src/fig3/plot_fig3f.py | data/Fig3/Fig3f.csv |
| 3g | alpha-catenin radial lifetime profile | src/fig3/plot_fig3g.py | data/Fig3/Fig3g.csv |
| 3k | Junction lifetime profile, MDCKII stretch (aCat) | src/fig3/plot_fig3k.py | data/Fig3/Fig3k_profiles.csv |
| 3l | Junctional lifetime anisotropy, MDCKII stretch (aCat) | src/fig3/plot_fig3l.py | data/Fig3/Fig3l.csv |
| 4e | HMR-1/E-cadherin lifetime, seam cells (perp. vs para.) | src/fig4/plot_fig4e.py | data/Fig4/Fig4e.csv |
| 4f | DEB-1/vinculin lifetime, pharynx vs body-wall muscle | src/fig4/plot_fig4f.py | data/Fig4/Fig4f.csv |
| 4g | LAM-2/laminin lifetime, pharynx vs body-wall muscle | src/fig4/plot_fig4g.py | data/Fig4/Fig4g.csv |
| 5c | Per-FOV in-FA lifetime, inhibitors (MEF) | src/fig5/plot_fig5c.py | data/Fig5/Fig5c.csv |
| 6b | Spermathecal lifetime contrast vs adult age | src/fig6/plot_Fig6b.py (stats: src/fig6/stats_Fig6b.py) | data/Fig6/Fig6b.csv |
| 6d | Spermathecal lifetime, control vs glp-1(ts) | src/fig6/plot_Fig6d.py | data/Fig6/Fig6d.csv |
| 6f | Spermathecal lifetime contrast vs ovulation number | src/fig6/plot_Fig6f.py | data/Fig6/Fig6f.csv |
| 6h | Spermathecal lifetime, control vs plc-1(RNAi) | src/fig6/plot_Fig6h.py | data/Fig6/Fig6h.csv |
| S1b | Relative brightness, 17 dyes | src/figS1/plot_figS1b.py | data/FigS1/FigS1b.csv |
| S1c | Mean lifetime, 17 dyes | src/figS1/plot_figS1c.py | data/FigS1/FigS1c.csv |
| S2a | Absorption/emission spectra | src/figS2/plot_figS2a.py | data/FigS2/FigS2a.csv |
| S2b | Photophysics table (no plot) | - | data/FigS2/FigS2b.csv |
| S2c | cpWHalo reversibility trace 1 | src/figS2/plot_figS2c.py | data/FigS2/FigS2c.csv |
| S2d | cpWHalo reversibility trace 2 | src/figS2/plot_figS2d.py | data/FigS2/FigS2d.csv |
| S3c | Representative pooled in-FA decay fit | src/figS3/plot_figS3c.py | data/FigS3/FigS3c.csv, FigS3c_meta.csv |
| S4g | E-cadherin per-FOV variance | src/figS4/plot_figS4g.py | data/FigS4/FigS4g.csv |
| S4h | alpha-catenin per-FOV variance | src/figS4/plot_figS4h.py | data/FigS4/FigS4h.csv |
| S4i | Radial lifetime, one image | src/figS4/plot_figS4i.py | data/FigS4/FigS4i.csv |
| S5c | Representative periphery lifetime profile, MDCKII stretch | src/figS5/plot_figS5c.py | data/FigS5/FigS5c.csv |
| S5d | Representative four-region lifetime, MDCKII stretch | src/figS5/plot_figS5d.py | data/FigS5/FigS5d.csv |
| S6d | Brood size per strain | src/figS6/plot_FigS6d.py | data/FigS6/FigS6d.csv |
| S6e | Embryonic viability per strain | src/figS6/plot_FigS6e.py | data/FigS6/FigS6e.csv |

`data/index.csv` is the data dictionary: every file, its columns and units, the
panel it maps to, and the number of rows.

### Panels not in this table
- Image and schematic panels (FLIM images, gels, cartoons, and the structure
  renders in Fig. 1b, S4b, S4e, S6b) have no source-data CSV. Structure renders
  are documented in `docs/pymol_render_notes.md`. This also covers the
  representative FLIM/microscopy image panels (Fig. 3h-j; 4a-d; 6a, 6c, 6e,
  6g, 6i, 6j; S5a, S5b; S6a) -- these are figure crops assembled directly
  for the manuscript, not files this package regenerates.

## Layout

```
data/          source-data CSVs, one per panel, grouped by figure; index.csv
src/           plotting scripts, one per panel; common/ shared style and paths
src/pipeline/  the upstream analysis: FLIM fitting and Cellpose segmentation,
               i.e. the code that produced the CSVs in data/ from the raw data
figures/       regenerated panels (SVG + PNG); generated output, not committed
               -- run make_all_figures.py to produce it
examples/      smoke-test runner, plus one small input and expected output per
               upstream analysis step
docs/          non-scriptable steps, raw-data note, PyMOL render notes
```

## The upstream analysis

The panels above regenerate from the committed CSVs alone. The code that produced
those CSVs from the raw FLIM data is in `src/pipeline/`:

| Step | Script |
|---|---|
| IRF extraction + per-pixel lifetime maps + LAS X validation | `src/pipeline/flim_pipeline.py` |
| Cellpose-SAM segmentation into label maps | `src/pipeline/cellpose_inference.py` |
| Per-ROI lifetime by NNLS on each ROI's pooled raw decay | `src/pipeline/flim_fa_fit.py` (snapshots, + heterogeneity test), `src/pipeline/flim_roi_fit.py` (time-lapse, shared labels + drift) |

These steps need the raw data (PTU files, LAS X spreadsheets, Cellpose weights),
which is deposited externally -- see `docs/RAW_DATA.md`. **None of it is needed to
reproduce a figure.** `src/pipeline/README.md` documents how the stages chain, the
two IRF-extraction conventions used in this study and which datasets used each,
and the exact per-dataset invocations.

`examples/` ships a small committed slice of real data for each of these steps
with its expected output, so the fit can be checked without a large download:

```
python examples/flim_roi_fit/check_flim_roi_fit.py    # refits a real 192x192
                                                      # photon-cube crop and
                                                      # reproduces the published
                                                      # per-FA lifetimes
python examples/segmentation/check_segmentation.py     # validates the committed
                                                      # Cellpose label mask
```

## Conventions
- Lifetime is the intensity-weighted mean fluorescence lifetime in ns; ROI
  values use a pooled-decay NNLS fit unless noted otherwise in a panel's
  script or `data/index.csv` (see `docs/NON_SCRIPTABLE_STEPS.md`).
- Vinculin constructs use paper names: Vin-WF, Vin-cpH, Vin-N-WF; C. elegans
  constructs are prefixed Ce (CeEcad-, CeVin-, CeLam-).
- Tests are two-sided and are run on independent cell-, FOV-, animal- or
  molecule-level summaries; pixels, focal adhesions, junction segments, radial
  bins and repeated time points within one cell or FOV are nested measurements,
  not biological replicates. The test used for each comparison is stated in that
  panel's script and in the figure legend. Across the package these are:
  Mann-Whitney U, Welch's t-test for unequal variances, paired t-tests for
  within-animal comparisons of two regions of the same worm, Holm-Bonferroni
  control of family-wise error across related comparisons, Kaplan-Meier survival
  with right-censoring for single-molecule dwell times (Fig. 1j), Fisher's exact
  test for the fraction of coupled cardiomyocytes, and a Monte-Carlo test against
  1,000 phase-randomized surrogates for lifetime-beating coupling (Fig. 2k,
  fixed seed 12345). Significance is p < 0.05 after adjustment where applicable.
- All two-group comparisons in Fig. 6 use two-sided Mann-Whitney U with
  Holm-Bonferroni correction, via the shared helpers in `src/fig6/_stats.py`.
  Fig. 6b's p-values are reported in its legend rather than drawn on the panel
  and are computed by `src/fig6/stats_Fig6b.py` (writes `Fig6b_stats.csv`).
  Fig. S6d and S6e deliberately keep Welch's t-test: both are ceiling-bounded
  measurements whose N2 reference piles up at the ceiling, where a rank test
  reports significance for biologically negligible shifts -- see the rationale
  in `welch_holm_vs_control()` in `src/figS6/_common.py`.

## Data availability
Source data are provided here and with the article. Raw and processed microscopy
data (FLIM PTU, images, Cellpose label masks and model weights, optical-tweezers
HDF5) are deposited at Figshare: [FIGSHARE ACCESSION - to be added on
publication]. See `docs/RAW_DATA.md`.

## Citation
[Wu et al., 2026, WHaloForce. Citation to be finalized on publication.]

The repository itself will be deposited at GitHub (accession to be added on
publication).

## License
Code is licensed under the MIT License (`LICENSE`). The source data in `data/`
are licensed under CC BY 4.0 (`LICENSE-data`).
