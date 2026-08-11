# Data availability and the raw data

## What is in this repository
The committed source data (data/) are the plotted values for each figure panel,
one tidy CSV per panel. These regenerate every plot via the scripts in src/.
They are small (a few MB total).

## What is not, and where it lives
The raw microscopy data are large (on the order of 200 GB across datasets) and
are not committed. This includes:
- FLIM photon data as PTU files and Leica LIF files,
- raw and rendered TIFF image stacks and Cellpose label masks,
- optical-tweezers HDF5 (.h5) files,
- Cellpose-SAM model weights.

Raw and processed data are deposited at Figshare:
[FIGSHARE ACCESSION - to be added on publication].

Source data for reproducing all plots are provided here and with the article.

## FLIM data convention
FLIM images exported from the pipeline are stored as paired TIFFs:
- <name>_ch0.tif: intensity image (photon counts, 16-bit), used for
  segmentation and as the photon weight for lifetime measurement.
- <name>_ch1.tif: pixel-wise lifetime stored as integer = tau_ns * 1000
  (picoseconds); divide by 1000 for nanoseconds. A value of 0 means undefined
  (fit failed or below the noise floor).

## Pipeline start point
The scripted analysis starts from the LAS X exports (PTU files plus the
per-image fit spreadsheets). See NON_SCRIPTABLE_STEPS.md for the upstream
instrument steps, and src/pipeline/README.md for the scripts that consume those
exports and the exact per-dataset invocations.

## What you can run without the raw data
Everything needed for the figures: the committed source-data CSVs regenerate all
41 panels. In addition, examples/ holds a small committed slice of the real raw
data for each analysis step that this package can check without a large download:

- examples/flim_roi_fit/ -- a 192x192 crop of a real photon cube (471 KB) plus
  its label map, IRF and LAS X lifetimes. The check refits it and reproduces the
  published per-FA lifetimes for those 45 focal adhesions.
- examples/segmentation/ -- the matching intensity image and the published
  Cellpose label mask, with a check that validates the mask and its registration.
  Re-running segmentation from scratch needs the model weights (above).

Steps that have no runnable example (LAS X fitting, Cellpose fine-tuning,
optical-tweezers trace extraction, structure renders) are listed with the reason
in examples/README.md.
