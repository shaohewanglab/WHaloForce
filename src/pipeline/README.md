# Upstream pipeline: FLIM fitting and segmentation

The scripts in `src/<fig>/` regenerate the plotted panels from the committed
source-data CSVs. The scripts in *this* folder are the steps that produced those
CSVs in the first place: they turn raw FLIM photon data into per-ROI lifetimes.

**These steps need the raw data, which is not in this repository.** PTU photon
files, LAS X fit spreadsheets and Cellpose model weights are ~200 GB in total and
live in the external deposit (see `docs/RAW_DATA.md`). Nothing here is needed to
reproduce a figure panel. The code is included so the analysis is inspectable and
re-runnable against the deposited raw data.

For a small, self-contained check of the fit core and the segmentation contract
that runs from committed inputs only, see `examples/flim_roi_fit/` and
`examples/segmentation/`.

## The five stages

```
  LAS X (GUI, not scripted)          -> tau_1..tau_3 per image + PTU export
    |                                   (docs/NON_SCRIPTABLE_STEPS.md)
    v
1 flim_pipeline.py                   -> IRF vector, per-pixel tau maps,
    |                                   intensity images, LAS X validation
    v
2 cellpose_inference.py              -> integer label maps (FAs / cells /
    |                                   junctions), QC overlays
    v
3 flim_fa_fit.py                     -> per-ROI and per-FOV lifetimes by NNLS
  or flim_roi_fit.py                    on each ROI's pooled raw decay
    |
    v
4 per-dataset aggregation and stats  -> the tidy CSVs committed under data/
    |                                   (panel-specific; see src/<fig>/)
    v
5 src/<fig>/plot_*.py                -> figures/
```

Stage 4 is deliberately not generalized: each panel aggregates differently
(per-FOV, per-junction, per-worm, radial profile, time series), and what it
produces is exactly the committed `data/FigX/FigXy.csv`. The aggregation actually used for each panel
is documented in that panel's plotting script and in `data/index.csv`.

## Entry points

| Script | What it does |
|---|---|
| `flim_pipeline.py` | Stages 1-3 of the fit: extract the IRF from a Rhodamine B reference by damped Richardson-Lucy deconvolution, fit per-pixel amplitudes with tau_1..tau_3 held fixed from LAS X, write intensity + lifetime maps, and validate per-image mean lifetimes against LAS X. |
| `cellpose_inference.py` | Run Cellpose-SAM 4.1.1 weights over the intensity images and write one uint16 label map per image, plus QC overlays. Requires the model weights. |
| `flim_fa_fit.py` | Per-ROI lifetime by exact NNLS on each ROI's **pooled** raw decay, plus a within-ROI heterogeneity test (shot-noise null, Moran's I) and a per-image pooled in-ROI row. Used for the snapshot focal-adhesion and junction datasets. |
| `flim_roi_fit.py` | The same pooled-decay NNLS without the heterogeneity test, but able to follow a **shared** label map across the frames of a time-lapse with per-frame drift correction (`--shared-labels-dir`, `--drift-csv`). Used for live/time-lapse datasets. |
| `FLIM_vis_v15.py` | Lifetime colormapping helper used by `flim_pipeline.py` for the QC images. Not called directly. |

Why the pooled decay: a per-pixel NNLS fit is biased on dim pixels (NNLS clamps
the slow-component amplitude to zero at low SNR), so averaging the per-pixel map
drags an ROI's lifetime toward its brighter sub-regions. Summing the raw decays
of an ROI's pixels and fitting once gives a far higher-SNR decay and an unbiased
`tau_int`. The pooled value is the principal ROI measurement throughout the
paper; the per-pixel map is for visualization and QC.

## IRF extraction

`extract_irf` recovers the IRF by deconvolving the reference's own decay out of
the measured Rhodamine B reference, then fits a global time-zero shift. Three
things have to be right, and each was wrong somewhere in the first pass:

**1. The reference must match the experiment's detector.** The `Channel` column
of every LAS X spreadsheet records it. The 3T3 and MEF focal-adhesion sessions
are HyD X 2; the 20260526 3T3-untreated session is HyD X 4. Detector matters more
than spectral window here: fitted Rhodamine B lifetimes differ by ~0.1 ns
between the two detectors for the same window. Match detector first, then
window; only then tune the settings below on chi-squared.

**2. The reference's decay model must match the reference.** Deconvolution needs
a model of the reference decay. A single exponential is right only for a
reference that fits well with one component. Rhodamine B aliquots stored at
-20 C in this study developed a small fast component (6-11% amplitude, 0.25-0.52
ns): LAS X fits them at chi2 ~49 with one component and ~1.9 with two. Pass the
two-component fit with `--ref-components AMP,TAU AMP,TAU`; deconvolving with a
mono-exponential kernel instead leaves a distorted IRF. The freshly prepared
2026-01-27 reference set is genuinely mono-exponential and needs only `--tau-r`.

**3. The session time-zero shift must be fitted, not assumed.** Each acquisition
session has its own laser t0 offset (`--irf-extra-shift-ns`), fitted once per
session by scanning chi-squared. Leaving it at 0 is the single largest error
source: it shows up as an oscillating residual at the rising edge of every
experimental fit. Fitted values here range from -0.08 to +0.16 ns.

`--irf-sigma-bins` damps the Richardson-Lucy iteration (plain RL is
semi-convergent on these high-SNR references and never settles). With the
reference modelled correctly, only light damping is needed: 0.1-0.2, not the 0.4
that an earlier pass used to compensate for the wrong kernel.

### Settings per session

Every session below uses the 20260529 reference set with its two-component LAS X
model. `chi2` is the median reduced chi-squared of the pooled in-ROI NNLS fit.

| Session | Panels | Detector | Reference | sigma | shift (ns) | chi2 before -> after |
|---|---|---|---|---|---|---|
| 20260526 3T3-untreated | 2c, 2d | HyD X 4 | `RhB-37C-HyDX4-554-558to650` | 0.1 | -0.06 | 18.3 -> 4.2 |
| 20260609 | 2e (LV), ED3c | HyD X 2 | `RhB-37C-HyDX2-561-565to650-1` | 0.2 | -0.08 | 44.7 -> 3.7 |
| 20260727 | 2e | HyD X 2 | `RhB-37C-HyDX2-561-565to650-1` | 0.1 | -0.02 | 46.5 -> 11.2 |
| 20260606 | 5c (Y5uM-12h) | HyD X 2 | `RhB-37C-HyDX2-561-565to650-1` | 0.1 | -0.02 | 61.5 -> 8.3 |
| 20260607 | 5c (drug-panel) | HyD X 2 | `RhB-37C-HyDX2-561-565to650-1` | 0.1 | -0.02 | 27.2 -> 3.5 |
| 20260608 | 5c (drug-panel2) | HyD X 2 | `RhB-37C-HyDX2-561-565to650-1` | 0.1 | **+0.16** | 192.2 -> 6.1 |

The MDCKII junction datasets (Fig. 3c/3d/3f/3g, ED4g-i) and the iCM datasets
(Fig. 2h-2k) were processed separately with damped RL at sigma 0.4 and are not
part of this table.

The staged-labeling sessions (Fig. 5d, Extended Data Fig. 7c,d,f-h; 20261013,
20261017, 20261024 and 20261027) use the same HyD X 2 reference
(`RhB-37C-HyDX2-561-565to650-1`, tau_R = 2.088 ns). The 20261027 decays were
0.12 ns later than the other sessions, measured on their leading edges, so that
session has an extra IRF time-zero shift of +0.10 ns. One FOV was excluded
before fitting (20261024 `3T3-185-stain-treat-DMSO-8`, Extended Data Fig. 7f).

An independent check on the 3T3-untreated choice: the 2026-01-27 HyD X 4 reference
at the same 554/558-650 window -- a different dye batch four months earlier, and
mono-exponential rather than two-component -- gives chi2 4.0 and a mean in-FA
lifetime of 3.0213 ns, against 3.0224 ns for the reference above. Two
independent references agreeing to 0.001 ns.

`--tau-r` is the reference's fitted lifetime and must be matched to `--irf-ptu`.
The reference is Rhodamine B at 2 uM in ethanol, acquired at 37 C or room
temperature with detector and spectral settings matched to the experiment. All
fitted values live in the spreadsheets under `data/RhB-reference/`.

## Which Cellpose model went with which images

Segmentation weights are as much a part of the analysis as the flags above, and
the 3T3 focal-adhesion datasets used **two** models within each dataset, chosen
per image by treatment -- as the Methods put it, "separate fine-tuned models were
used when perturbations substantially altered object morphology". The `model`
column of each deposited `cellpose_labels/<dataset>/inference_log.csv` records
what was actually used, one row per image; that log is the authority. For the six
DMSO-vs-blebbistatin datasets (`dmso-bleb20-{0727-,LV}{185,186,188}`, behind
Fig. 2c,2d,2e):

| Images | Model |
|---|---|
| DMSO | `FA_cpsam_3T3` -- Cellpose-SAM fine-tuned on hand-corrected 3T3 focal-adhesion crops |
| blebbistatin 20 uM | `FA-3T3-bleb20` -- a further fine-tune on blebbistatin-treated images, whose adhesions are larger and more elongated |

Using the DMSO model on a blebbistatin image is a silent error: it still returns
plausible adhesions, just fewer and differently shaped ones, which changes every
lifetime measured inside them. `cellpose_inference.py` therefore refuses to run
when an image matches no `--model` spec or more than one.

## Worked example: the Fig. 2e focal-adhesion lifetimes

Paths below are inside the raw-data deposit. This is the invocation for one of
the three constructs of the 2026-07-27 DMSO-vs-blebbistatin acquisition; the
other constructs and the Fig. 5c inhibitor panel differ only in the spreadsheet,
PTU folder and `--include-path-regex`.

```
# Stage 1 -- IRF + per-pixel maps + LAS X validation
python src/pipeline/flim_pipeline.py \
  --spreadsheets 20260727-3T3-bulk-185-186-188-bleb20uM-12h.xlsx \
  --ptu-folder   20260727-3T3-bulk-185-186-188-bleb20uM-12h.sptw \
  --irf-ptu      20260529-RhB.sptw/RhB-37C-HyDX2-561-565to650-1.ptu \
  --tau-r 2.088 --irf-sigma-bins 0 --rl-iters 30 --irf-rl-tol 0 \
  --output-dir   outputs/dmso-bleb20-0727-185

# Stage 2 -- Cellpose-SAM focal-adhesion segmentation (needs the weights).
# The published run used a DIFFERENT fine-tune per condition, so the two specs
# below are not interchangeable: DMSO images get the base focal-adhesion model,
# blebbistatin images get the blebbistatin fine-tune. Anchoring each token on
# '-' or a string end keeps them mutually exclusive; the script refuses to run
# if an image matches none or more than one.
python src/pipeline/cellpose_inference.py \
  --input-dirs outputs/dmso-bleb20-0727-185 \
  --model '(^|-)DMSO(-|$)=<weights>/FA_cpsam_3T3' \
          '(^|-)bleb(-|$)=<weights>/FA-3T3-bleb20' \
  --outdir cellpose_labels/dmso-bleb20-0727-185

# Stage 3 -- per-FA and per-FOV lifetimes (same IRF flags as stage 1)
python src/pipeline/flim_fa_fit.py \
  --spreadsheets 20260727-3T3-bulk-185-186-188-bleb20uM-12h.xlsx \
  --ptu-folder   20260727-3T3-bulk-185-186-188-bleb20uM-12h.sptw \
  --labels-dir   cellpose_labels/dmso-bleb20-0727-185 \
  --irf-ptu      20260529-RhB.sptw/RhB-37C-HyDX2-561-565to650-1.ptu \
  --tau-r 2.088 --irf-sigma-bins 0 --rl-iters 30 --irf-rl-tol 0 \
  --output-dir   results/dmso-bleb20-0727-185
```

Stage 3 writes `per_roi_measurements.csv` (one row per focal adhesion) and
`per_FOV_summary.csv` (one row per image, i.e. one field of view (FOV): NNLS on
the pooled decay of all in-FA pixels; different FOVs contain different cells).
This per-FOV in-FA `tau_int` is what Fig. 2c, 2e, 5c, 5d and Extended Data Fig. 7
plot; the per-FA SD within each FOV is what Fig. 2d plots. (Analysis outputs
made before 2026-10-31, including any in the raw-data deposit, name this file
`per_cell_roi.csv`; the contents are the same.)

## Reproducibility of the fit

Re-running stage 3 on macOS against the published `per_roi_measurements.csv`
(produced on Windows) reproduces it to floating-point precision:

- `area_px`, `n_pixels_used`, `n_photons_total`, `het_n_px`: exact.
- `tau_int_ns`, `tau_amp_ns`, `chi2_red`, `mean_intensity`, and the observed
  within-ROI statistics `het_var_obs` / `moran_I_obs`: agree to 1e-5, which is
  the committed CSV's own `%.5f` rounding. Per-FOV pooled `tau_inFA_nnls_ns`
  matched exactly.
- `alpha1..alpha3`: agree to <= 1e-5 relative. The three-component amplitude
  split is mildly ill-conditioned (the fast component is close to the IRF), so
  the individual amplitudes are the least stable quantity here. `tau_int` is the
  validated one.
- The columns compared against the shot-noise Monte-Carlo null
  (`het_var_null_med`, `het_var_ratio`, `het_var_pval`, `moran_I_null_med`,
  `moran_I_z`, `moran_I_pval`) are **not** reproducible ROI-by-ROI. The RNG is
  seeded (`--seed`, default 0) but drawn once per run and advanced across every
  ROI of every image, so replaying them requires re-running the identical image
  set in the identical order. An identical whole-image re-run agrees to ~1e-2 on
  `het_var_ratio`; a re-run over a subset of the ROIs (as in the example check)
  diverges freely, up to ~0.2. No figure panel plots these columns; they support
  statements about within-FA heterogeneity, which are made at the
  distribution level over hundreds of FAs.

`examples/flim_roi_fit/check_flim_roi_fit.py` encodes these tolerances as an
executable check against a committed slice of the real data.

## Environment

`environment.yml` builds an env with everything: `ptufile`, `lmfit`, `tqdm`,
`tifffile`, `scipy` for the FLIM steps, and `cellpose` 4.1.1 for segmentation.
Segmentation wants a CUDA GPU; install the CUDA torch build separately if you
have one (`--no-gpu` works but is slow). The plotting-only `requirements.txt`
does not cover these steps.
