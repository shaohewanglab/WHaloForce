"""Run a trained Cellpose model on every intensity image and save label maps.

This is the segmentation step of the upstream pipeline: it turns the intensity
images written by flim_pipeline.py into the integer label maps that
flim_fa_fit.py / flim_roi_fit.py measure lifetimes inside. The same script
segments focal adhesions, cells and cell-cell junctions -- what it finds is set
by the model passed to --model, not by the code.

Accepts either of two per-image naming conventions:
  * `<name>_ch0.tif`       -- LAS X export (uint16 photons)
  * `<name>_intensity.tif` -- flim_pipeline.py output (float32 photons)

REQUIRES cellpose 4.1.1 (Cellpose-SAM) and torch; a CUDA GPU is strongly
recommended but --no-gpu works. See environment.yml. Model weights are ~1.2 GB
each and are NOT in this repository -- they are in the raw-data deposit (see
docs/RAW_DATA.md), which is why this step is only smoke-checked in
examples/segmentation/ against a committed reference mask.

    # one model for every image
    python src/pipeline/cellpose_inference.py \
        --input-dirs outputs/<dataset> \
        --model <weights>/FA_cpsam_3T3 \
        --outdir cellpose_labels/<dataset>

    # a different fine-tune per treatment, as the published 3T3 runs used
    python src/pipeline/cellpose_inference.py \
        --input-dirs outputs/<dataset> \
        --model '(^|-)DMSO(-|$)=<weights>/FA_cpsam_3T3' \
                '(^|-)bleb(-|$)=<weights>/FA-3T3-bleb20' \
        --outdir cellpose_labels/<dataset>

Outputs:
  <outdir>/<name>_labels.tif   uint16 label map per source image
  <outdir>/qc/<name>.png       red-boundary overlay on grayscale intensity
  <outdir>/inference_log.csv   one row per image: n_ROIs, elapsed_s, ...

Then measure lifetimes inside those labels:

    python src/pipeline/flim_fa_fit.py --labels-dir <outdir> ...
"""
import argparse
import glob
import os
import re
import time

import numpy as np
import tifffile
import matplotlib.image as mpimg

from cellpose import models


def _find_boundaries(labels):
    """Numpy-only outer boundary mask: pixel != one of its 4-neighbours."""
    L = labels
    b = np.zeros_like(L, dtype=bool)
    b[:-1, :] |= L[:-1, :] != L[1:, :]
    b[1:, :]  |= L[:-1, :] != L[1:, :]
    b[:, :-1] |= L[:, :-1] != L[:, 1:]
    b[:, 1:]  |= L[:, :-1] != L[:, 1:]
    return b & (L > 0)

HERE = os.path.dirname(os.path.abspath(__file__))    # src/pipeline
REPO = os.path.dirname(os.path.dirname(HERE))        # repository root

# Suffixes that identify an intensity image. Order matters: try the longer
# (`_intensity.tif`) before the shorter (`_ch0.tif`) so glob assignments are
# unambiguous when both happen to coexist.
INTENSITY_SUFFIXES = ("_intensity.tif", "_ch0.tif")


def _basename_for(path):
    """Strip the intensity suffix to recover the per-image base name."""
    bn = os.path.basename(path)
    for sfx in INTENSITY_SUFFIXES:
        if bn.endswith(sfx):
            return bn[: -len(sfx)]
    return os.path.splitext(bn)[0]


def _list_intensity_images(input_dirs):
    """Yield (path, base_name) for every intensity image in input_dirs."""
    seen = set()
    for d in input_dirs:
        for sfx in INTENSITY_SUFFIXES:
            for p in sorted(glob.glob(os.path.join(d, f"*{sfx}"))):
                bn = _basename_for(p)
                key = (os.path.abspath(d), bn)
                if key in seen:
                    continue
                seen.add(key)
                yield p, bn, os.path.basename(d)


def _gray_for_qc(im):
    x = im.astype(float)
    pos = x[x > 0]
    if not pos.size:
        return np.zeros_like(x)
    lo, hi = np.percentile(pos, [1, 99.5])
    return np.clip((x - lo) / max(hi - lo, 1e-6), 0, 1)


def save_qc(name, intensity, labels, qcdir):
    gray = _gray_for_qc(intensity)
    rgb = np.dstack([gray] * 3)
    rgb[_find_boundaries(labels)] = [1, 0, 0]
    mpimg.imsave(os.path.join(qcdir, f"{name}.png"), rgb)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", nargs="+", required=True, metavar="SPEC",
                    help="Cellpose model weights. Either one path, applied to "
                         "every image, or several 'REGEX=PATH' specs selected by "
                         "matching REGEX (case-insensitive) against the image "
                         "name -- the published 3T3 runs used a per-condition "
                         "fine-tune this way, e.g. --model "
                         "'(^|-)DMSO(-|$)=W/FA_cpsam_3T3' "
                         "'(^|-)bleb(-|$)=W/FA-3T3-bleb20'. Anchor the token on "
                         "'-' or a string end as shown: a bare 'bleb' would also "
                         "match a 'Bleb-5uM' image in a spreadsheet that mixes "
                         "doses. An image that matches no spec, or more than one, "
                         "is an error rather than a silent choice.")
    ap.add_argument("--input-dirs", nargs="+", required=True,
                    help="one or more folders to scan for *_ch0.tif or "
                         "*_intensity.tif images (typically a flim_pipeline "
                         "--output-dir)")
    ap.add_argument("--outdir", required=True,
                    help="output directory for labels and QC")
    ap.add_argument("--no-gpu", action="store_true",
                    help="force CPU inference (default uses GPU if available)")
    ap.add_argument("--flow-threshold", type=float, default=0.4,
                    help="cellpose flow threshold (lower = stricter)")
    ap.add_argument("--cellprob-threshold", type=float, default=0.0,
                    help="cellpose cellprob threshold; lower => more masks")
    ap.add_argument("--min-size", type=int, default=5,
                    help="drop ROIs below this px area (default 5)")
    ap.add_argument("--limit", type=int, default=None,
                    help="process at most N images (smoke test)")
    args = ap.parse_args()

    qcdir = os.path.join(args.outdir, "qc")
    os.makedirs(qcdir, exist_ok=True)

    # Resolve the model selection strategy:
    #   * one bare path            -> that single model for every image
    #   * one or more REGEX=PATH   -> first spec whose regex matches the name
    if len(args.model) == 1 and "=" not in args.model[0]:
        model_specs = [(None, args.model[0])]
    else:
        model_specs = []
        for spec in args.model:
            if "=" not in spec:
                raise SystemExit(
                    f"--model {spec!r}: mixing a bare path with REGEX=PATH specs "
                    f"is ambiguous; give either one path or all REGEX=PATH")
            rx, _, path = spec.partition("=")
            model_specs.append((re.compile(rx, re.IGNORECASE), path))
        print("Per-condition models: "
              + ", ".join(f"{rx.pattern} -> {os.path.basename(p)}"
                          for rx, p in model_specs))

    _model_cache = {}

    def get_model_for(name):
        """Return (CellposeModel, model_path) for this image, loading lazily.

        Selecting the wrong fine-tune silently changes every downstream lifetime,
        so an image that matches no spec -- or several -- stops the run instead of
        falling back to a default or to whichever spec came first.
        """
        hits = [(rx, path) for rx, path in model_specs
                if rx is None or rx.search(name)]
        if not hits:
            raise SystemExit(
                f"no --model spec matches image {name!r}; patterns tried = "
                f"{[rx.pattern for rx, _ in model_specs]}")
        if len(hits) > 1:
            raise SystemExit(
                f"image {name!r} matches {len(hits)} --model specs "
                f"({[rx.pattern for rx, _ in hits]}); make the patterns "
                f"mutually exclusive (anchor the condition token on '-' or a "
                f"string end) so the choice is unambiguous")
        mp = hits[0][1]
        if mp not in _model_cache:
            print(f"Loading model: {mp}")
            _model_cache[mp] = models.CellposeModel(
                gpu=not args.no_gpu, pretrained_model=mp)
        return _model_cache[mp], mp

    log_rows = []
    n_done = 0
    for ipath, name, folder_label in _list_intensity_images(args.input_dirs):
        if args.limit is not None and n_done >= args.limit:
            break
        model, model_path = get_model_for(name)
        intensity = tifffile.imread(ipath)
        t0 = time.time()
        # Cellpose normalizes internally (default percentile [1, 99] on a
        # per-image basis), matching how the FA_cpsam model was trained
        # on raw 16-bit crops. Float intensities work too.
        masks, flows, styles = model.eval(
            intensity,
            channels=None,
            flow_threshold=args.flow_threshold,
            cellprob_threshold=args.cellprob_threshold,
            min_size=args.min_size,
        )
        dt = time.time() - t0
        n_fas = int(masks.max())
        tifffile.imwrite(os.path.join(args.outdir, f"{name}_labels.tif"),
                         masks.astype(np.uint16))
        save_qc(name, intensity, masks, qcdir)
        log_rows.append({
            "image": name, "n_ROIs": n_fas,
            "shape_y": intensity.shape[0], "shape_x": intensity.shape[1],
            "elapsed_s": round(dt, 2),
            "folder": folder_label,
            "model": os.path.basename(model_path),
        })
        n_done += 1
        print(f"  {name:40s}  ROIs={n_fas:5d}  {dt:5.2f}s")

    if not log_rows:
        raise SystemExit(
            f"No intensity images found in: {args.input_dirs}. Expected "
            f"files ending in {INTENSITY_SUFFIXES}.")

    # Write inference log
    import csv
    log_path = os.path.join(args.outdir, "inference_log.csv")
    with open(log_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=log_rows[0].keys())
        w.writeheader()
        w.writerows(log_rows)
    total_fas = sum(r["n_ROIs"] for r in log_rows)
    total_time = sum(r["elapsed_s"] for r in log_rows)
    print(f"\nProcessed {len(log_rows)} images, total {total_fas} ROIs, "
          f"{total_time:.1f}s on {'GPU' if not args.no_gpu else 'CPU'}.")
    print(f"Labels: {args.outdir}/")
    print(f"Next:   python src/pipeline/flim_fa_fit.py "
          f"--labels-dir {args.outdir} ...")


if __name__ == "__main__":
    main()
