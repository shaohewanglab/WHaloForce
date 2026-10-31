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

Last run: 48/48 panels and 2/2 step checks passed (2026-10-31, Python 3.12
with the packages in `requirements.txt`).

## Figure panel to script and source data

| Panel | Plotted quantity | Script | Source data |
|---|---|---|---|
| 1d | Mean lifetime (int/amp weighted), cpHalo/cpWHalo/dye | src/fig1/plot_fig1d.py | data/Fig1/Fig1d.csv |
| 1f | cpHalo-554 force-ramp trace | src/fig1/plot_fig1f.py | data/Fig1/Fig1f.csv |
| 1g | cpWHalo-554 force-ramp trace | src/fig1/plot_fig1g.py | data/Fig1/Fig1g.csv |
| 1h | cpWHalo-554 force-clamp trace | src/fig1/plot_fig1h.py | data/Fig1/Fig1h.csv |
| 1i | Quenched-state dwell times (box) | src/fig1/plot_fig1i.py | data/Fig1/Fig1i.csv |
| 1j | Dwell-time survival | src/fig1/plot_fig1j.py | data/Fig1/Fig1j.csv |
| 2c | Per-FOV in-FA lifetime, untreated | src/fig2/plot_fig2c.py | data/Fig2/Fig2c.csv |
| 2d | Per-FOV SD of per-FA lifetime | src/fig2/plot_fig2d.py | data/Fig2/Fig2d.csv |
| 2e | Per-FOV in-FA lifetime, DMSO vs blebbistatin | src/fig2/plot_fig2e.py | data/Fig2/Fig2e.csv |
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
| 5d | Staged labeling, per-FOV in-FA lifetime (MEF Vcl-WF, JFX554) | src/fig5/plot_fig5d.py | data/Fig5/Fig5d.csv |
| 6b | Spermathecal lifetime contrast vs adult age | src/fig6/plot_Fig6b.py (stats: src/fig6/stats_Fig6b.py) | data/Fig6/Fig6b.csv |
| 6d | Spermathecal lifetime, control vs glp-1(ts) | src/fig6/plot_Fig6d.py | data/Fig6/Fig6d.csv |
| 6f | Spermathecal lifetime contrast vs ovulation number | src/fig6/plot_Fig6f.py | data/Fig6/Fig6f.csv |
| 6h | Spermathecal lifetime, control vs plc-1(RNAi) | src/fig6/plot_Fig6h.py | data/Fig6/Fig6h.csv |
| 6j | Spermathecal lifetime contrast, L4 vs day-3 labeling, day 4 | src/fig6/plot_Fig6j.py | data/Fig6/Fig6j.csv |
| ED1b | Relative brightness, 17 dyes | src/extdatafig1/plot_ExtDataFig1b.py | data/ExtDataFig1/ExtDataFig1b.csv |
| ED1c | Mean lifetime, 17 dyes | src/extdatafig1/plot_ExtDataFig1c.py | data/ExtDataFig1/ExtDataFig1c.csv |
| ED2a | Absorption/emission spectra | src/extdatafig2/plot_ExtDataFig2a.py | data/ExtDataFig2/ExtDataFig2a.csv |
| ED2b | Photophysics table (no plot) | - | data/ExtDataFig2/ExtDataFig2b.csv |
| ED2c | cpWHalo reversibility trace 1 | src/extdatafig2/plot_ExtDataFig2c.py | data/ExtDataFig2/ExtDataFig2c.csv |
| ED2d | cpWHalo reversibility trace 2 | src/extdatafig2/plot_ExtDataFig2d.py | data/ExtDataFig2/ExtDataFig2d.csv |
| ED3c | Representative pooled in-FA decay fit | src/extdatafig3/plot_ExtDataFig3c.py | data/ExtDataFig3/ExtDataFig3c.csv, ExtDataFig3c_meta.csv |
| ED4g | E-cadherin per-FOV variance | src/extdatafig4/plot_ExtDataFig4g.py | data/ExtDataFig4/ExtDataFig4g.csv |
| ED4h | alpha-catenin per-FOV variance | src/extdatafig4/plot_ExtDataFig4h.py | data/ExtDataFig4/ExtDataFig4h.csv |
| ED4i | Radial lifetime, one image | src/extdatafig4/plot_ExtDataFig4i.py | data/ExtDataFig4/ExtDataFig4i.csv |
| ED5c | Representative periphery lifetime profile, MDCKII stretch | src/extdatafig5/plot_ExtDataFig5c.py | data/ExtDataFig5/ExtDataFig5c.csv |
| ED5d | Representative four-region lifetime, MDCKII stretch | src/extdatafig5/plot_ExtDataFig5d.py | data/ExtDataFig5/ExtDataFig5d.csv |
| ED6d | Brood size per strain | src/extdatafig6/plot_ExtDataFig6d.py | data/ExtDataFig6/ExtDataFig6d.csv |
| ED6e | Embryonic viability per strain | src/extdatafig6/plot_ExtDataFig6e.py | data/ExtDataFig6/ExtDataFig6e.csv |
| ED7c | Staged labeling (MEF Vcl-WF, JF549) | src/extdatafig7/plot_ExtDataFig7c.py | data/ExtDataFig7/ExtDataFig7c.csv |
| ED7d | Staged labeling (MEF Vcl-WF, JF571) | src/extdatafig7/plot_ExtDataFig7d.py | data/ExtDataFig7/ExtDataFig7d.csv |
| ED7f | Staged labeling (NIH-3T3 Vin-WF, JFX554) | src/extdatafig7/plot_ExtDataFig7f.py | data/ExtDataFig7/ExtDataFig7f.csv |
| ED7g | Staged labeling (NIH-3T3 Vin-WF, JF549) | src/extdatafig7/plot_ExtDataFig7g.py | data/ExtDataFig7/ExtDataFig7g.csv |
| ED7h | Staged labeling (NIH-3T3 Vin-WF, JF571) | src/extdatafig7/plot_ExtDataFig7h.py | data/ExtDataFig7/ExtDataFig7h.csv |

"ED" is Extended Data (formerly Supplementary Fig. S1-S6, same order and panel
letters). `data/index.csv` is the data dictionary: every file, its columns and
units, the panel it maps to, and the number of rows.

### Panels not in this table
- Image and schematic panels (FLIM images, gels, cartoons, timelines, and the
  structure renders in Fig. 1b, ED4b, ED4e, ED6b) have no source-data CSV.
  Structure renders are documented in `docs/pymol_render_notes.md`. This also
  covers the representative FLIM/microscopy image panels (Fig. 3h-j; 4a-d; 5b,
  5e, 5f; 6a, 6c, 6e, 6g, 6i; ED3a; ED5a, ED5b; ED6a; ED7a, ED7e) and the
  labeling timelines (Fig. 5d left, 6i; ED7b) -- these are assembled directly
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
  Fig. 6j is a single comparison (two-sided Mann-Whitney U, no correction).
  Fig. 6f reports Pearson r with a 95% percentile CI from 50,000 whole-worm
  bootstrap resamples (seed 20261029), because 11 of the 42 worms contribute
  both spermathecae; no independent-point p value is shown.
  Extended Data Fig. 6d and 6e deliberately keep Welch's t-test: both are
  ceiling-bounded measurements whose N2 reference piles up at the ceiling, where
  a rank test reports significance for biologically negligible shifts -- see the
  rationale in `welch_holm_vs_control()` in `src/extdatafig6/_common.py`.
- FA lifetime panels (Fig. 2c-e, 5c, 5d; ED7c,d,f-h) plot one value per field of
  view (FOV); different FOVs contain different cells. The staged-labeling panels
  (Fig. 5d, ED7) share `src/common/staged.py`: two-sided Mann-Whitney U on FOVs
  pooled across experiments, Holm-Bonferroni over the two Bleb vs DMSO tests of
  a panel and, separately, over its two schedule tests.

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
