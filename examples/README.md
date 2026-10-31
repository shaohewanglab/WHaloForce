# Examples and smoke tests

```
python examples/run_examples.py
```

runs everything below: the upstream step checks, then every figure panel.

## Upstream analysis steps

One small committed input, its expected output, and a `check_*.py` per step. All
of them run from committed data only -- no raw microscopy files, no model weights.

| Folder | Step | What the check proves |
|---|---|---|
| `flim_roi_fit/` | Per-ROI FLIM lifetime fit (pooled-decay NNLS) | The shipped fit code reproduces the **published** per-FA lifetimes for a real image region, to the source CSV's own rounding |
| `segmentation/` | Cellpose-SAM label maps | The committed label mask satisfies the contract the fit depends on and is correctly registered to its intensity image |

Each folder's README explains where the input came from, why the expected output
is trustworthy, and how to run the full-scale step against the raw-data deposit.

The steps themselves live in `src/pipeline/`; see `src/pipeline/README.md` for how
they chain together and the per-dataset invocations.

### Steps with no example, and why

* **LAS X decay fitting** -- a GUI step on the instrument. It produces the
  tau_1..tau_3 the scripted fit holds fixed. See `docs/NON_SCRIPTABLE_STEPS.md`.
  The `taus` array in `flim_roi_fit/input/*.npz` is a real LAS X output.
* **Cellpose model fine-tuning** -- done in the Cellpose GUI on hand-corrected
  crops; the resulting weights are in the raw-data deposit.
* **Optical-tweezers force traces** (Fig. 1f-j, ED2c, ED2d) -- extracted from LUMICKS
  C-Trap HDF5 upstream of this package. The committed CSVs *are* the extracted
  traces, so the panel scripts already exercise everything this package owns.
  `src/fig1/QU_transition_analysis.py` documents the dwell-time analysis.
* **Structure renders, sequence alignment, super-resolution processing** -- GUI or
  separate-instrument steps; see `docs/NON_SCRIPTABLE_STEPS.md` and
  `docs/pymol_render_notes.md`.

## Figure panels

`run_examples.py` also runs every `src/<fig>/plot_*.py` and confirms it writes its
`figures/<Panel>.svg`, which validates that each committed source-data CSV still
regenerates its panel. That is the same work `make_all_figures.py` does, with
pass/fail reporting.

## Requirements

The panels and the FLIM fit check need only `pip install -r requirements.txt`.
The segmentation check additionally needs `tifffile`, and its optional
`--model` inference path needs `cellpose` 4.1.1; both are in `environment.yml`. A
check whose dependencies are missing prints `[SKIP]` and does not fail the run.
