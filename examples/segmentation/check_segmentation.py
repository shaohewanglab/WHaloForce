#!/usr/bin/env python3
"""Check the segmentation step: one small intensity image and its label mask.

Why this check is shaped differently from the other one
-------------------------------------------------------
Segmentation is the one pipeline step that cannot be re-run from committed data:
the Cellpose-SAM 4.1.1 weights are ~1.2 GB each and live in the raw-data deposit
(docs/RAW_DATA.md). So the committed artifact is the label mask itself -- the same
mask the published lifetimes were measured in -- and this check verifies the
properties the rest of the pipeline actually depends on.

Inference does reproduce: with the correct weights, re-running the full image on
CPU rather than the original CUDA GPU recovered all 131 ROIs at median IoU 1.000.
The --model path below is still reported rather than asserted, for two reasons
documented in this folder's README: it depends on picking the fine-tune that
matches the image's treatment (this is a blebbistatin image, so FA-3T3-bleb20,
not FA_cpsam_3T3), and running on this small crop instead of the full frame
changes Cellpose's per-image intensity normalization.

By default it runs three groups of checks with no model and no GPU:

  1. Contract    -- the shape, dtype and label-ID conventions that
                    src/pipeline/flim_fa_fit.py and flim_roi_fit.py assume.
  2. Coherence   -- the mask is real and correctly registered to the image: the
                    labelled pixels are much brighter than the background, and no
                    small translation of the mask fits the image better than the
                    committed one. A shifted, stale or mismatched mask fails this.
  3. Consistency -- this mask is byte-identical to the label array shipped inside
                    examples/flim_roi_fit/input/185-bleb-1_crop.npz, so the two
                    examples cannot drift apart.

With `--model <weights>` it additionally re-runs real Cellpose inference on the
committed intensity crop and reports per-ROI agreement (IoU) against the
committed mask. That path needs the cellpose environment and the weights, and is
reported, not gated -- see the note on reproducibility above.

Usage:
    python examples/segmentation/check_segmentation.py
    python examples/segmentation/check_segmentation.py --model /path/to/FA-3T3-bleb20
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

IMAGE = "185-bleb-1_crop"
# Measured properties of the committed mask. These are the values the published
# focal-adhesion segmentation of this crop has; they are asserted so a corrupted,
# resaved or regenerated mask does not slip through unnoticed.
EXPECT_SHAPE = (192, 192)
EXPECT_N_ROIS = 45
EXPECT_MIN_AREA = 5        # cellpose_inference.py --min-size default
MIN_BRIGHTNESS_RATIO = 3.0  # mean intensity inside labels / outside
MAX_SHIFT_PROBE = 3        # px; identity must beat every shift in this window
# The weights that produced the committed mask, per the deposited
# cellpose_labels/dmso-bleb20-0727-185/inference_log.csv. 185-bleb-1 is a
# blebbistatin image, so it is the blebbistatin fine-tune, not the DMSO model.
PUBLISHED_MODEL = "FA-3T3-bleb20"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default=None,
                    help="optional path to Cellpose weights; when given, re-runs "
                         "inference on the committed crop and reports IoU. For "
                         "this image use the blebbistatin fine-tune "
                         "(FA-3T3-bleb20), not the DMSO model")
    ap.add_argument("--no-gpu", action="store_true",
                    help="force CPU inference for the optional --model run")
    args = ap.parse_args()

    try:
        import tifffile
    except ImportError:
        print("[SKIP] tifffile is not installed; it is in environment.yml "
              "(conda env create -f environment.yml).")
        return 0

    ipath = HERE / "input" / f"{IMAGE}_ch0.tif"
    lpath = HERE / "expected" / f"{IMAGE}_labels.tif"
    inten = tifffile.imread(ipath)
    labels = tifffile.imread(lpath)
    print(f"input    {ipath.relative_to(REPO)}")
    print(f"expected {lpath.relative_to(REPO)}")
    print(f"         (segmented in the published run with {PUBLISHED_MODEL})")

    failures = []

    def check(ok: bool, name: str, detail: str):
        print(f"  [{'OK ' if ok else 'FAIL'}] {name:22s} {detail}")
        if not ok:
            failures.append(name)

    # ---- 1. contract ----------------------------------------------------- #
    check(labels.shape == EXPECT_SHAPE and inten.shape == EXPECT_SHAPE,
          "shape", f"intensity {inten.shape}, labels {labels.shape}, "
                   f"expected {EXPECT_SHAPE}")
    check(labels.dtype == np.uint16, "label dtype",
          f"{labels.dtype} (uint16: what cellpose_inference.py writes)")
    ids = np.setdiff1d(np.unique(labels), [0])
    check(len(ids) == EXPECT_N_ROIS, "n_ROIs",
          f"{len(ids)} ROIs, expected {EXPECT_N_ROIS}")
    check(ids.min() >= 1, "background is 0",
          f"lowest non-zero label = {ids.min()} (0 must mean not-an-ROI)")
    areas = np.bincount(labels.ravel())[1:]
    areas = areas[areas > 0]
    check(areas.min() >= EXPECT_MIN_AREA, "min ROI area",
          f"smallest ROI {areas.min()} px >= --min-size {EXPECT_MIN_AREA} "
          f"(median {int(np.median(areas))}, max {areas.max()} px)")
    # The crop deliberately keeps the label IDs of the full 512x512 image, which
    # is what lets the expected CSV in examples/flim_roi_fit/ be a literal slice
    # of the published per_roi_measurements.csv. So IDs are NOT 1..N here.
    print(f"  [   ] label ID range        {ids.min()}..{ids.max()}, "
          f"non-contiguous by design (IDs preserved from the full image)")

    # ---- 2. coherence with the intensity image ---------------------------- #
    x = inten.astype(float)
    m = labels > 0
    ratio = x[m].mean() / max(x[~m].mean(), 1e-9)
    check(ratio >= MIN_BRIGHTNESS_RATIO, "labels on signal",
          f"mean intensity inside/outside = {ratio:.2f} "
          f"(>= {MIN_BRIGHTNESS_RATIO}); {100 * m.mean():.1f}% of px labelled")
    base = x[m].mean()
    alt_shift, alt_val = None, -np.inf
    for dy in range(-MAX_SHIFT_PROBE, MAX_SHIFT_PROBE + 1):
        for dx in range(-MAX_SHIFT_PROBE, MAX_SHIFT_PROBE + 1):
            if dy == dx == 0:
                continue
            v = x[np.roll(np.roll(m, dy, 0), dx, 1)].mean()
            if v > alt_val:
                alt_shift, alt_val = (dy, dx), v
    check(base > alt_val, "mask registration",
          f"committed mask scores {base:.1f}; best shifted alternative "
          f"{alt_shift} scores {alt_val:.1f} "
          f"(no shift within +/-{MAX_SHIFT_PROBE} px fits better)")

    # ---- 3. consistency with the FLIM example ---------------------------- #
    npz = REPO / "examples" / "flim_roi_fit" / "input" / "185-bleb-1_crop.npz"
    if npz.exists():
        z = np.load(npz, allow_pickle=False)
        check(np.array_equal(labels, z["labels"]), "matches FLIM example",
              f"identical to the labels in {npz.name}")
        check(np.array_equal(inten, z["photon_cube"].sum(axis=2).astype(np.uint16)),
              "intensity == cube sum",
              "the intensity image is the photon cube summed over time bins")
    else:
        print(f"  [   ] matches FLIM example {npz.name} not found; skipped")

    # ---- 4. optional: real inference ------------------------------------- #
    if args.model:
        try:
            from cellpose import models
        except ImportError as e:
            print(f"  [   ] live inference        SKIP ({e}); needs the cellpose env")
        else:
            print(f"  ---- live Cellpose inference with {args.model} "
                  f"(reported, not gated) ----")
            if PUBLISHED_MODEL not in Path(args.model).name:
                print(f"  [!! ] WARNING: the committed mask came from "
                      f"{PUBLISHED_MODEL}, but --model points at "
                      f"{Path(args.model).name!r}. Low IoU below most likely "
                      f"means the wrong fine-tune, not a pipeline problem.")
            model = models.CellposeModel(gpu=not args.no_gpu,
                                         pretrained_model=args.model)
            got, _, _ = model.eval(inten, channels=None, flow_threshold=0.4,
                                   cellprob_threshold=0.0, min_size=5)
            ious = []
            for i in ids:
                a = labels == i
                overlap = np.bincount(got[a].ravel())
                overlap[0] = 0
                if not overlap.any():
                    ious.append(0.0)
                    continue
                j = overlap.argmax()
                b = got == j
                ious.append((a & b).sum() / (a | b).sum())
            ious = np.array(ious)
            print(f"  [   ] n_ROIs                {int(got.max())} found vs "
                  f"{len(ids)} committed")
            print(f"  [   ] per-ROI IoU           median {np.median(ious):.3f}, "
                  f"mean {ious.mean():.3f}, "
                  f"{(ious >= 0.5).sum()}/{len(ious)} matched at IoU>=0.5")

    if failures:
        print(f"\nFAILED: {len(failures)} check(s): {failures}")
        return 1
    print("\nOK: the committed segmentation mask is valid and correctly "
          "registered.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
