#!/usr/bin/env python3
"""Per-ROI FLIM fit by NNLS on the ROI's summed decay (bias-free).

Rationale
---------
flim_pipeline.py writes a per-pixel tau_int map by per-pixel NNLS. That
per-pixel fit is biased on dim pixels (NNLS clamps the slow-component
amplitude to zero in low-SNR pixels, dragging mean tau_int_map a few tenths
of a ns below the true global value -- see the diff between mean_tau_int
and mean_tau_int_map columns in flim_pipeline's summary.csv).

For per-ROI measurements (cells, FAs, any segmented region), the correct
approach is to **pool the raw PTU decay over the ROI's pixels and run NNLS
once on that summed decay**, holding the LAS-X-fitted tau_1..tau_3 fixed.
The summed decay has far higher SNR than any single pixel, so NNLS doesn't
clamp components, and the recovered tau_int is unbiased.

Inputs
------
* The PTU folder + spreadsheet that flim_pipeline.py already ran on -- so we
  reuse the same tau_1..tau_3 per image and the same row -> PTU mapping
  logic (Region filter, t-index for time-lapse).
* The IRF vector saved by flim_pipeline.py (`outputs/.../irf_vector.npy`).
* A folder of label maps (`<labels-dir>/<image_name>_labels.tif`) where
  each ROI has a unique integer >= 1 and 0 = not-an-ROI.

The label maps don't need to be FAs -- they can be cell segmentations,
mitochondria, anything. The script is structure-agnostic.

Photon floor
------------
Per-pixel filtering (`--min-photons-per-pixel`) only removes pixels with
near-zero counts (default 5). Inside a properly-segmented ROI we trust the
mask, so the per-pixel floor only kills pixels that are essentially detector
dark (no fluorescence). Don't use 50 here -- that would bias toward the
bright sub-regions of each ROI. The relevant SNR check happens per-ROI
(`--min-photons-per-roi`, default 500) and per-area (`--min-area-px`, default 5).

Output
------
`<output-dir>/per_roi_measurements.csv` with one row per ROI:
    image_name, path, t_index, roi_id, area_px, n_pixels_used,
    n_photons_total, mean_intensity, tau_int_ns, tau_amp_ns,
    alpha1, alpha2, alpha3, alpha_bg, chi2_red, taus_used
Plus `per_roi_per_image_summary.csv` with one row per image.
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
# tifffile is imported lazily in main(), where the label maps are read, so
# that the pure-numpy fit functions above can be imported and checked with
# only numpy/scipy/pandas installed (see examples/flim_roi_fit/).
from scipy.optimize import nnls
from tqdm import tqdm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flim_pipeline as fp  # type: ignore  # noqa: E402


log = logging.getLogger("flim_roi")


def setup_logging(output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%H:%M:%S")
    sh = logging.StreamHandler(sys.stdout); sh.setFormatter(fmt); log.addHandler(sh)
    fh = logging.FileHandler(output_dir / "flim_per_roi.log"); fh.setFormatter(fmt)
    log.addHandler(fh)


def fit_per_roi(cube: np.ndarray, labels: np.ndarray,
                A: np.ndarray, keep: np.ndarray, taus: np.ndarray,
                min_photons_per_pixel: float,
                min_photons_per_roi: float,
                min_area_px: int) -> list[dict]:
    """Fit one NNLS per ROI label in `labels` (uint16 label map, 0 = bg).

    Parameters
    ----------
    cube      : (Y, X, H) photon histograms, post-IRF-bg-subtraction if any.
    labels    : (Y, X) integer label map, ROI ids >= 1, 0 = background.
    A         : (H, 4) design matrix.
    keep      : (H,) boolean keep-bins mask.
    taus      : (3,) array of fixed lifetimes (ns).
    """
    H, W, nbins = cube.shape
    if labels.shape != (H, W):
        raise ValueError(f"label shape {labels.shape} != cube ({H},{W})")
    # Per-pixel intensity and inclusion mask.
    pix_intensity = cube.sum(axis=2)
    bright_pix = pix_intensity >= float(min_photons_per_pixel)
    Af = A[keep]
    dof = max(int(keep.sum()) - 4, 1)
    rows = []
    roi_ids = np.unique(labels)
    roi_ids = roi_ids[roi_ids > 0]
    for rid in roi_ids:
        roi_full = labels == rid
        area_px = int(roi_full.sum())
        if area_px < min_area_px:
            continue
        used = roi_full & bright_pix
        n_used = int(used.sum())
        if n_used == 0:
            continue
        decay = cube[used].sum(axis=0).astype(np.float64)
        n_phot = float(decay.sum())
        if n_phot < min_photons_per_roi:
            continue
        try:
            x, _ = nnls(Af, decay[keep])
        except Exception:
            continue
        a = x[:3]
        if (a * taus).sum() <= 0:
            tau_int = float("nan")
            tau_amp = float("nan")
        else:
            tau_int = float((a * taus**2).sum() / (a * taus).sum())
            tau_amp = float((a * taus).sum() / a.sum()) if a.sum() > 0 else float("nan")
        model = Af @ x
        chi2 = float(np.sum((decay[keep] - model)**2 /
                            np.clip(model, 1.0, None)) / dof)
        rows.append({
            "roi_id": int(rid),
            "area_px": area_px,
            "n_pixels_used": n_used,
            "n_photons_total": n_phot,
            "mean_intensity": float(pix_intensity[roi_full].mean()),
            "tau_int_ns": tau_int,
            "tau_amp_ns": tau_amp,
            "alpha1": float(x[0]),
            "alpha2": float(x[1]),
            "alpha3": float(x[2]),
            "alpha_bg": float(x[3]),
            "chi2_red": chi2,
        })
    return rows


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spreadsheets", nargs="+", type=Path, required=True,
                    help="LAS X spreadsheet(s) with fitted tau_i values")
    ap.add_argument("--ptu-folder", type=Path, required=True,
                    help="folder of source PTU files")
    ap.add_argument("--labels-dir", type=Path, required=False,
                    help="folder of <image_name>_labels.tif label maps "
                         "(per-frame; use this OR --shared-labels-dir)")
    ap.add_argument("--shared-labels-dir", type=Path, required=False,
                    help="folder with ONE <Path>_labels.tif per Position "
                         "(shared across all live frames of that Position; "
                         "label IDs are preserved so the same roi_id tracks "
                         "the same cell across time). Use with --drift-csv "
                         "to shift the mask per live frame.")
    ap.add_argument("--drift-csv", type=Path, default=None,
                    help="drift_vectors.csv (from estimate_xy_drift.py). When "
                         "set together with --shared-labels-dir, the t0 mask "
                         "is shifted by -(shift_y, shift_x) per frame before "
                         "ROI summation.")
    ap.add_argument("--irf", type=Path, required=True,
                    help="path to saved irf_vector.npy (from a flim_pipeline run)")
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--min-photons-per-pixel", type=float, default=5.0,
                    help="skip pixels with <N photons (default 5; only excludes "
                         "near-zero pixels -- segmentation handles bg)")
    ap.add_argument("--min-photons-per-roi", type=float, default=500.0,
                    help="drop ROIs whose post-filter total < N photons (default 500)")
    ap.add_argument("--min-area-px", type=int, default=5,
                    help="drop ROIs smaller than this in pixels (default 5)")
    ap.add_argument("--t0-shift-ns", type=float, default=-0.0023,
                    help="IRF t0 shift used by the matching flim_pipeline run "
                         "(see its log; default uses the typical fitted value)")
    ap.add_argument("--exclude-ns", nargs="+", default=None, metavar="LO-HI",
                    help="time range(s) in ns to drop from the fit, e.g. '7.0-7.6'")
    ap.add_argument("--skip-tindex", nargs="+", type=int, default=None,
                    metavar="N",
                    help="for time-lapse runs, drop these t-index values "
                         "(typically '--skip-tindex 1' if flim_pipeline used it)")
    ap.add_argument("--limit", type=int, default=None,
                    help="process at most N images (debugging)")
    return ap.parse_args(argv)


def main(argv=None):
    import tifffile
    args = parse_args(argv)
    setup_logging(args.output_dir)
    log.info("Per-ROI FLIM fit; output -> %s", args.output_dir)

    irf = np.load(args.irf).astype(np.float64)
    n = int(irf.size)
    # If a sidecar irf_meta.json is present (written by recent flim_pipeline
    # runs), pull t0_shift from there so the user doesn't have to remember.
    import json
    meta_path = args.irf.parent / "irf_meta.json"
    if meta_path.exists():
        with open(meta_path) as f:
            meta = json.load(f)
        t0_shift_ns = float(meta.get("t0_shift_ns", args.t0_shift_ns))
        log.info("Loaded IRF from %s (n=%d, t0_shift=%+.4f ns from irf_meta.json)",
                 args.irf, n, t0_shift_ns)
    else:
        t0_shift_ns = float(args.t0_shift_ns)
        log.info("Loaded IRF from %s (n=%d, t0_shift=%+.4f ns from --t0-shift-ns)",
                 args.irf, n, t0_shift_ns)

    # Spreadsheet (Region filter + t-index added automatically by flim_pipeline).
    df = fp.load_spreadsheets(args.spreadsheets)
    time_lapse = fp.T_INDEX_COL in df.columns
    if time_lapse and args.skip_tindex:
        skip = set(int(x) for x in args.skip_tindex)
        before = len(df)
        df = df[~df[fp.T_INDEX_COL].isin(skip)].copy()
        log.info("  --skip-tindex %s: %d -> %d rows", sorted(skip), before, len(df))
    if args.limit:
        df = df.iloc[:args.limit]
    log.info("Processing %d images", len(df))

    exclude_ranges = []
    for tok in (args.exclude_ns or []):
        lo, hi = (float(v) for v in tok.split("-"))
        exclude_ranges.append((lo, hi))

    # Modes: per-frame labels (--labels-dir) OR shared per-Position labels
    # (--shared-labels-dir, optionally + --drift-csv to shift the mask per frame).
    use_shared = args.shared_labels_dir is not None
    if use_shared and args.labels_dir is not None:
        raise SystemExit("Use either --labels-dir or --shared-labels-dir, not both")
    if not use_shared and args.labels_dir is None:
        raise SystemExit("Provide --labels-dir or --shared-labels-dir")
    drift_df = None
    if args.drift_csv is not None:
        if not use_shared:
            log.warning("--drift-csv has no effect without --shared-labels-dir; ignoring")
        else:
            drift_df = pd.read_csv(args.drift_csv)
            log.info("Loaded drift vectors from %s (%d rows)",
                     args.drift_csv, len(drift_df))
    shared_label_cache: dict[str, np.ndarray] = {}
    from scipy.ndimage import shift as ndshift  # type: ignore

    def lookup_labels(position_name: str, image_name: str, t_index_local):
        """Resolve the (possibly drift-shifted) label map for one frame.

        Returns (labels: np.ndarray, label_source: str) or (None, reason).
        """
        if use_shared:
            if position_name not in shared_label_cache:
                p = args.shared_labels_dir / f"{position_name}_labels.tif"
                if not p.exists():
                    return None, f"no shared labels at {p}"
                shared_label_cache[position_name] = tifffile.imread(p)
            base = shared_label_cache[position_name]
            if drift_df is None or t_index_local is None or t_index_local == 0:
                return base, "shared (no drift)"
            d = drift_df[(drift_df["position"] == position_name) &
                         (drift_df["t_index"] == t_index_local)]
            if d.empty:
                return base, "shared (no drift row)"
            dy = float(d["shift_y_px"].iloc[0])
            dx = float(d["shift_x_px"].iloc[0])
            # Mask is in t0 coordinates; live-frame features at p are at p+shift
            # in t0 (shift = ref - mov). To get the mask in LIVE coordinates,
            # shift it by -shift. order=0 to preserve integer label IDs.
            shifted = ndshift(base, shift=(-dy, -dx), order=0,
                              mode="constant", cval=0, prefilter=False)
            return shifted.astype(base.dtype), f"shared (drift dy={dy:+.2f},dx={dx:+.2f})"
        else:
            p = args.labels_dir / f"{image_name}_labels.tif"
            if not p.exists():
                return None, f"no labels at {p}"
            return tifffile.imread(p), str(p)

    out_rows = []
    img_rows = []
    irf_info_cache = None
    for _, row in tqdm(df.iterrows(), total=len(df), desc="images"):
        name = str(row[fp.COL_PATH])
        t_index = int(row[fp.T_INDEX_COL]) if time_lapse else None
        ptu = fp.find_ptu_for_path(args.ptu_folder, name, t_index=t_index)
        if ptu is None:
            label = f"{name}_t{t_index}" if t_index else name
            log.warning("  no PTU for '%s'; skipping", label)
            continue
        image_name = ptu.stem
        labels, src = lookup_labels(name, image_name, t_index)
        if labels is None:
            log.warning("  no labels for %s (%s); skipping", image_name, src)
            continue
        if labels.max() == 0:
            img_rows.append({"image_name": image_name, "n_rois": 0,
                             "n_rois_fitted": 0})
            continue
        cube, n_period, dt_ns, _ = fp.read_ptu_cube(ptu)
        if irf_info_cache is None:
            t_ns = np.arange(min(n, n_period)) * dt_ns
            irf_info_cache = {
                "irf": irf[:min(n, n_period)],
                "t_ns": t_ns, "dt_ns": dt_ns,
                "n": min(n, n_period),
                "t0_shift_ns": t0_shift_ns,
            }
        if n_period != irf_info_cache["n"]:
            log.warning("  %s: period %d != IRF/cube cached %d; cropping",
                        image_name, n_period, irf_info_cache["n"])
        n_use = irf_info_cache["n"]
        cube = cube[..., :n_use].astype(np.float64)
        if fp.COL_IRF_BG in row:
            irf_bg = float(row.get(fp.COL_IRF_BG, 0.0) or 0.0)
            if irf_bg:
                cube = cube - irf_bg
                cube = np.clip(cube, 0, None)
        taus = np.array([float(row[c]) for c in fp.COL_TAU])
        irf_shift_ns = float(row.get(fp.COL_IRF_SHIFT, 0.0) or 0.0)
        A = fp.build_basis(irf_info_cache, taus, irf_shift_ns)[:n_use]
        keep = fp.build_keep_mask(irf_info_cache["t_ns"][:n_use], exclude_ranges)
        roi_rows = fit_per_roi(
            cube, labels, A, keep, taus,
            args.min_photons_per_pixel,
            args.min_photons_per_roi,
            args.min_area_px,
        )
        for r in roi_rows:
            r.update({"image_name": image_name,
                      "path": name,
                      "t_index": t_index if t_index else "",
                      "tau1_ns": float(taus[0]),
                      "tau2_ns": float(taus[1]),
                      "tau3_ns": float(taus[2]),
                      "lasx_mean_tau_int_ns": float(
                          row.get(fp.COL_MEAN_INT_LASX, np.nan))})
        out_rows.extend(roi_rows)
        img_rows.append({"image_name": image_name,
                         "n_rois": int(labels.max()),
                         "n_rois_fitted": len(roi_rows)})

    args.output_dir.mkdir(parents=True, exist_ok=True)
    if out_rows:
        cols_order = ["image_name", "path", "t_index", "roi_id", "area_px",
                      "n_pixels_used", "n_photons_total", "mean_intensity",
                      "tau_int_ns", "tau_amp_ns",
                      "alpha1", "alpha2", "alpha3", "alpha_bg", "chi2_red",
                      "tau1_ns", "tau2_ns", "tau3_ns", "lasx_mean_tau_int_ns"]
        out_df = pd.DataFrame(out_rows)
        out_df = out_df[[c for c in cols_order if c in out_df.columns]]
        out_csv = args.output_dir / "per_roi_measurements.csv"
        out_df.to_csv(out_csv, index=False, float_format="%.5f")
        log.info("Wrote %d ROI rows -> %s", len(out_df), out_csv)
    else:
        log.warning("No ROI rows produced.")

    img_df = pd.DataFrame(img_rows)
    img_csv = args.output_dir / "per_roi_per_image_summary.csv"
    img_df.to_csv(img_csv, index=False)
    log.info("Wrote per-image summary -> %s", img_csv)
    log.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
