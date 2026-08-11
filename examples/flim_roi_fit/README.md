# Example: per-ROI FLIM lifetime fit

One small input, its expected output, and a check that regenerates it.

```
python examples/flim_roi_fit/check_flim_roi_fit.py
```

## What step this is

Stage 3 of the upstream pipeline (`src/pipeline/README.md`): given a photon
histogram cube, an IRF, the LAS X lifetimes tau_1..tau_3, and a label map, measure
the fluorescence lifetime of each labelled region. Every in-ROI lifetime in the
paper comes from this step.

The measurement is a **pooled-decay** fit: the raw decays of all pixels in an ROI
are summed and fitted once by non-negative least squares against the reconvolution
model, with tau_1..tau_3 held fixed and only the amplitudes free. The
intensity-weighted mean lifetime `tau_int = sum(a_i tau_i^2) / sum(a_i tau_i)` is
the reported quantity. Fitting the pooled decay rather than averaging a per-pixel
lifetime map matters: per-pixel NNLS clamps the slow-component amplitude to zero
on dim pixels, which biases an ROI average toward its brighter sub-regions.

The check calls `src/pipeline/flim_pipeline.build_basis` and
`src/pipeline/flim_fa_fit.fit_image` directly -- the same functions the published
run used, not a re-implementation.

## Input

`input/185-bleb-1_crop.npz` (471 KB) is real data: a 192x192 pixel crop of the
photon cube of image `185-bleb-1` (Vin-WF, blebbistatin 20 uM, from the
`dmso-bleb20-0727-185` dataset behind Fig. 2e), with everything the fit needs:

| Key | Meaning |
|---|---|
| `photon_cube` | uint16 (192, 192, 132): photon counts per pixel per TCSPC time bin |
| `labels` | uint16 (192, 192): Cellpose focal-adhesion label map, 0 = background |
| `irf` | float64 (132,): the instrument response function |
| `t_ns` | float64 (132,): time axis in ns |
| `taus` | float64 (3,): the LAS X lifetimes tau_1..tau_3 for this image, held fixed |
| `dt_ns`, `n_period` | TCSPC bin width and bins per laser period |
| `t0_shift_ns`, `irf_shift_ns`, `irf_bg_cnts` | IRF time-zero shift, the per-image LAS X IRF shift, and the LAS X IRF background |
| `lasx_mean_tau_int_ns` | the LAS X whole-image mean lifetime, for reference |
| `crop_origin_yx`, `crop_size`, `source_image` | provenance of the crop |

The IRF stored here was extracted with this dataset's exact settings: the
detector-matched Rhodamine B reference `RhB-37C-HyDX2-561-565to650-1`, its
two-component LAS X decay model, `--irf-sigma-bins 0.1` and a session t0 shift of
-0.02 ns. Those are the values in `flim-refit-2026-08-07/parameters/` and in
`src/pipeline/README.md`.

## Expected output

`expected/per_roi_measurements.csv` is a **verbatim 45-row slice of the published
`per_roi_measurements.csv`** for this image, not a value regenerated for this
example. That works because the crop window was chosen so all 45 focal adhesions
inside it are fully contained, none clipped by the boundary. Each ROI's fit only
ever sums that ROI's own pixels, so cropping cannot change any of their fits.

The one column that is *not* comparable is the per-image pooled in-FA lifetime
(`per_cell_roi.csv` in a real run): that pools every in-FA pixel of the full
512x512 image, so the crop's value legitimately differs. The check reports it and
sanity-checks its range instead of comparing it.

## Running the full step instead

The complete entry point reads the PTU directly. The 512x512 `185-bleb-1.ptu`
(~14 MB), the LAS X spreadsheet and the Rhodamine B reference PTU are all in the
raw-data deposit (`docs/RAW_DATA.md`); with those in place:

```
python src/pipeline/flim_fa_fit.py \
  --spreadsheets 20260727-3T3-bulk-185-186-188-bleb20uM-12h.xlsx \
  --ptu-folder   20260727-3T3-bulk-185-186-188-bleb20uM-12h.sptw \
  --labels-dir   cellpose_labels/dmso-bleb20-0727-185 \
  --irf-ptu      20260529-RhB.sptw/RhB-37C-HyDX2-561-565to650-1.ptu \
  --tau-r 2.088 --irf-sigma-bins 0 --rl-iters 30 --irf-rl-tol 0 \
  --include-path-regex '^185-bleb-1$' \
  --output-dir   /tmp/one-image
```

## Requirements

numpy, scipy and pandas only -- `pip install -r requirements.txt` is enough. The
committed `.npz` means no PTU reader (`ptufile`) or IRF fitter (`lmfit`) is needed;
the check prints `[SKIP]` rather than failing if a dependency is missing.
