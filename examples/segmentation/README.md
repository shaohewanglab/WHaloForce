# Example: Cellpose-SAM segmentation

One small intensity image, the label mask that was measured in it, and a check
that validates the mask without needing the model weights.

```
python examples/segmentation/check_segmentation.py

# optional: re-run real inference and report agreement (needs the cellpose env
# and the ~1.2 GB weights from the raw-data deposit). This is a blebbistatin
# image, so the model to use is the blebbistatin fine-tune, NOT FA_cpsam_3T3 --
# see "Which model" below.
python examples/segmentation/check_segmentation.py --model /path/to/FA-3T3-bleb20
```

## What step this is

Stage 2 of the upstream pipeline (`src/pipeline/README.md`): turn the intensity
image written by `flim_pipeline.py` into an integer label map, which
`flim_fa_fit.py` / `flim_roi_fit.py` then measure lifetimes inside. The same
script and the same procedure segment focal adhesions, whole cells and cell-cell
junctions; what is segmented is set by the model weights, not by the code.

Segmentation used Cellpose-SAM 4.1.1, either the base model or a model fine-tuned
in the Cellpose GUI on hand-corrected crops from these datasets. Fine-tuning is a
GUI step and is described in `docs/NON_SCRIPTABLE_STEPS.md`.

## Files

| File | What it is |
|---|---|
| `input/185-bleb-1_crop_ch0.tif` | uint16, 192x192. Photon-count intensity image: a crop of the real image `185-bleb-1` (Vin-WF, blebbistatin 20 uM, dataset `dmso-bleb20-0727-185`, behind Fig. 2e). Byte-identical to the photon cube in `examples/flim_roi_fit/` summed over its time bins. |
| `expected/185-bleb-1_crop_labels.tif` | uint16, 192x192. The **published** focal-adhesion label map for that region -- the actual mask the Fig. 2e lifetimes were measured in, cropped, not regenerated. 45 ROIs, all fully inside the crop. |

Label IDs are preserved from the full 512x512 image, so they run 34..102 rather
than 1..45. That is deliberate: it is what lets the expected CSV in
`examples/flim_roi_fit/` be a literal slice of the published
`per_roi_measurements.csv`. Downstream code only requires that 0 means
background and that each ROI has a distinct positive integer.

## Which model segmented this image

The 3T3 focal-adhesion datasets used **two** fine-tunes, chosen per image by
treatment, so "the Cellpose model" is not a single thing:

| Images | Model |
|---|---|
| DMSO | `FA_cpsam_3T3` |
| blebbistatin 20 uM -- **including this one** | `FA-3T3-bleb20` |

`185-bleb-1` is a blebbistatin image, so its published mask came from
`FA-3T3-bleb20`. The authority for any image is the `model` column of the
deposited `cellpose_labels/<dataset>/inference_log.csv`, which records one row per
image.

## Why the mask is the committed artifact

Only because the weights are large: ~1.2 GB per model, so they are deposited with
the raw data rather than committed here. Inference itself reproduces well.
Re-running the **correct** weights (`FA-3T3-bleb20`) on the full 512x512 image on
macOS CPU (cellpose 4.0.6, torch 2.7.1) -- a different Cellpose version and
device than the published run, which used 4.1.1 on a CUDA GPU -- recovers **131 of
131 ROIs with median per-ROI IoU 1.000** (mean 0.985, 130/131 matched at
IoU >= 0.5).

Two caveats that the `--model` path here will show:

- **Use the model that matches the treatment.** These datasets used a different
  fine-tune per condition (see the table above). Segmenting this blebbistatin
  image with the DMSO model instead returns 97 ROIs at median IoU 0.82 -- still
  plausible-looking adhesions, just fewer and differently shaped, which would
  change every lifetime measured inside them. `cellpose_inference.py` refuses to
  run when an image matches no `--model` spec or more than one, so this cannot
  happen silently there.
- **The crop is not the frame.** Running the correct model on this committed
  192x192 crop rather than the full image gives 42 ROIs vs 45, median IoU 0.875.
  Cellpose normalizes intensity per image, so a crop presents a different dynamic
  range and loses context at its edges. This is a property of the example's small
  size, not of the pipeline.

Because of both, the `--model` path reports agreement rather than asserting it.
The default check does not need the weights at all.

## What the default check verifies

No model, no GPU, no weights -- only the committed files:

1. **Contract** -- shape, `uint16` dtype, 0-means-background, distinct positive
   IDs, and every ROI at or above the `--min-size 5` floor. These are the
   assumptions `flim_fa_fit.py` and `flim_roi_fit.py` make about a label map.
2. **Coherence and registration** -- labelled pixels average 5.6x the background
   intensity, and no translation of the mask within +/-3 px fits the intensity
   image better than the committed alignment. A stale, shifted or mismatched mask
   fails this.
3. **Consistency** -- the mask is byte-identical to the `labels` array inside
   `examples/flim_roi_fit/input/185-bleb-1_crop.npz`, and the intensity image is
   exactly that example's photon cube summed over time, so the two examples
   cannot drift apart.

## Running the full step instead

With the weights in place (see `docs/RAW_DATA.md`):

```
python src/pipeline/cellpose_inference.py \
  --input-dirs outputs/dmso-bleb20-0727-185 \
  --model '(^|-)DMSO(-|$)=<weights>/FA_cpsam_3T3' \
          '(^|-)bleb(-|$)=<weights>/FA-3T3-bleb20' \
  --outdir cellpose_labels/dmso-bleb20-0727-185
```

Hence the two `REGEX=PATH` specs: images whose name carries the `DMSO` token get
the base focal-adhesion model, those carrying `bleb` get the blebbistatin
fine-tune. Anchoring each token on `-` or a string end keeps them mutually
exclusive -- a bare `bleb` would also match a `Bleb-5uM` image in a spreadsheet
that mixes doses. The script exits rather than guess if an image matches no spec
or more than one, and records the model it used per image in
`<outdir>/inference_log.csv`.

## Requirements

The default check needs `tifffile` (in `environment.yml`); it prints `[SKIP]`
rather than failing if it is missing. The optional `--model` run needs
`cellpose` (4.1.1 was used for the published segmentation) and torch.
