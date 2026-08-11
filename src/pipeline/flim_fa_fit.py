#!/usr/bin/env python3
"""Per-FA FLIM fit by NNLS on each FA's summed decay (bias-free), plus a
within-FA heterogeneity test against a shot-noise null.

Why
---
`flim_pipeline.py` writes a per-pixel tau_int map by per-pixel NNLS. That
per-pixel fit is biased low on dim pixels (NNLS clamps the slow-component
amplitude to zero at low SNR), so an intensity-weighted average of the
per-pixel map drags the per-FA tau toward the brighter sub-regions. As in
flim_roi_fit.py (this folder), the correct per-FA lifetime is a fresh NNLS on the
FA's *pooled raw decay* (tau_1..tau_3 held fixed from the LAS X spreadsheet). The
summed decay has far higher SNR than any single pixel, so NNLS does not clamp
components and the recovered tau_int is unbiased. A deliberately permissive
per-pixel floor (default 1 photon) keeps the FA-level tau from being biased toward
bright pixels; SNR is enforced at the FA level instead (--min-photons-per-roi,
default 200; --min-area-px, default 5).

This script is for snapshot datasets: it adds the heterogeneity test and the
per-image pooled in-FA row, but drops the per-frame time-index and shared-label /
drift machinery that flim_roi_fit.py has for time-lapse data. Use flim_roi_fit.py
for live datasets. See README.md in this folder for which datasets used which.

Within-FA heterogeneity (pattern vs noise)
------------------------------------------
For each FA we also ask: is the pixel-to-pixel tau scatter larger than finite
photons alone would produce? We build a per-FA shot-noise null by treating the
FA as a single uniform lifetime (its NNLS amplitude ratios), Poisson-resampling
each pixel's *actual* photon budget from that uniform model, re-estimating
per-pixel tau, and recomputing the intensity-weighted within-FA variance. B
such draws give a null distribution; the observed variance is compared to it.

Crucially the per-pixel tau estimator used for BOTH the observed map and the
null is the same fast clamp estimator (one unconstrained lstsq via a
precomputed pseudo-inverse, then clamp negative amplitudes to zero). This
reproduces the intensity-dependent NNLS clamping bias in the synthetic data
too, so that artifact cancels in the observed-vs-null comparison and is not
mistaken for biology. The per-FA *mean* (tau_int_ns) still uses exact NNLS on
the pooled decay -- the clamp estimator is only the test statistic's estimator.

Reported per FA:
  * het_var_ratio  = observed within-FA tau variance / null-median variance
                     (>1 => real over-dispersion beyond shot noise)
  * het_var_pval   = right-tailed MC p-value of the observed variance
  * moran_I_obs    = spatial autocorrelation of the within-FA tau map, with a
  * moran_I_z/pval   bias-aware null (Moran's I computed on each MC draw too).
                     Distinguishes spatially patterned heterogeneity (gradient/
                     domains) from spatially disordered over-dispersion.

Output
------
<output-dir>/per_roi_measurements.csv  -- one row per FA
<output-dir>/per_cell_roi.csv          -- one row per image: NNLS on the pooled
                                          decay of ALL in-FA pixels (bias-free
                                          per-cell in-FA tau)

Run (needs ptufile + lmfit, same as flim_pipeline.py -- see environment.yml):

    python src/pipeline/flim_fa_fit.py \
        --spreadsheets <dataset>.xlsx \
        --ptu-folder <dataset>.sptw \
        --labels-dir cellpose_labels/<dataset> \
        --irf-ptu RhB-37C-HyDX2-561-565to650-1.ptu --tau-r 2.088 \
        --output-dir results/<dataset>

See src/pipeline/README.md for the per-dataset invocations and where the inputs
come from, and examples/flim_roi_fit/ for a small runnable check of the fit core.
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

HERE = Path(__file__).resolve().parent
# flim_pipeline.py is a sibling in this same folder; reuse its IRF extraction,
# spreadsheet loading and basis construction so both steps agree exactly.
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


# --------------------------------------------------------------------------- #
# Per-pixel tau estimator used INSIDE the heterogeneity test.
# Fast, vectorised, and applied identically to observed and null decays so its
# (intensity-dependent) bias cancels in the comparison. NOT used for the
# reported per-FA mean, which uses exact NNLS on the pooled decay.
# --------------------------------------------------------------------------- #
def _clamp_tau_int(decays_keep: np.ndarray, pinvT: np.ndarray,
                   taus: np.ndarray) -> np.ndarray:
    """tau_int for a stack of decays via lstsq-then-clamp.

    decays_keep : (..., n_keep) photon histograms over the kept bins.
    pinvT       : (n_keep, 4) = pinv(A_keep).T.
    taus        : (3,) fixed lifetimes.
    Returns tau_int with shape decays_keep.shape[:-1]; NaN where undefined.
    """
    amps = decays_keep @ pinvT                      # (..., 4) unconstrained
    a = np.clip(amps[..., :3], 0.0, None)           # clamp the 3 lifetimes
    sat = (a * taus).sum(axis=-1)                   # sum alpha_i * tau_i
    sat2 = (a * taus**2).sum(axis=-1)               # sum alpha_i * tau_i^2
    with np.errstate(divide="ignore", invalid="ignore"):
        tau = np.where(sat > 0, sat2 / sat, np.nan)
    return tau


def _wvar(tau: np.ndarray, w: np.ndarray, axis=-1):
    """Intensity-weighted variance of tau along `axis`.

    This is the same within-ROI variance definition the per-dataset aggregation
    (stage 4) uses when it summarises pixel-level scatter inside an ROI."""
    wsum = w.sum(axis=axis)
    mean = (w * tau).sum(axis=axis) / wsum
    var = (w * (tau - np.expand_dims(mean, axis))**2).sum(axis=axis) / wsum
    return var


def _moran_edges(used: np.ndarray):
    """4-connectivity edge list among the True pixels of `used` (2D bool).

    Returns (idx_i, idx_j) integer arrays indexing into the flattened list of
    used pixels in row-major order (the same order as np.where(used))."""
    order = -np.ones(used.shape, dtype=np.int64)
    ys, xs = np.where(used)
    order[ys, xs] = np.arange(ys.size)
    ii, jj = [], []
    # right neighbour
    right = used[:, :-1] & used[:, 1:]
    ry, rx = np.where(right)
    ii.append(order[ry, rx]); jj.append(order[ry, rx + 1])
    # down neighbour
    down = used[:-1, :] & used[1:, :]
    dy, dx = np.where(down)
    ii.append(order[dy, dx]); jj.append(order[dy + 1, dx])
    if not ii:
        return np.array([], int), np.array([], int)
    return np.concatenate(ii), np.concatenate(jj)


def _moran_I(tau: np.ndarray, ei: np.ndarray, ej: np.ndarray):
    """Moran's I (binary symmetric weights) for tau of shape (..., K).

    ei, ej index the K pixels; edges are undirected (counted once here). With
    S0 = 2E, I = (K / S0) * sum_{ordered ij} z_i z_j / sum z^2
            = (K / E) * sum_{edges} z_i z_j / sum z^2.
    Returns I with shape tau.shape[:-1]; NaN where undefined."""
    K = tau.shape[-1]
    E = ei.size
    if E == 0 or K < 3:
        return np.full(tau.shape[:-1], np.nan)
    z = tau - tau.mean(axis=-1, keepdims=True)
    num = (z[..., ei] * z[..., ej]).sum(axis=-1)      # sum over edges
    den = (z**2).sum(axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        I = np.where(den > 0, (K / E) * num / den, np.nan)
    return I


def heterogeneity_test(decay_k: np.ndarray, n_phot_k: np.ndarray,
                       used: np.ndarray, Af: np.ndarray, a_roi: np.ndarray,
                       taus: np.ndarray, B: int, rng) -> dict:
    """Shot-noise null for within-FA tau over-dispersion + spatial pattern.

    decay_k   : (K, n_keep) per-pixel decays over kept bins (K used pixels).
    n_phot_k  : (K,) per-pixel photon totals over kept bins (intensity weight).
    used      : (Y, X) bool mask of the K pixels (for adjacency).
    Af        : (n_keep, 4) design matrix over kept bins.
    a_roi     : (4,) NNLS amplitudes from the FA's pooled decay (uniform model).
    """
    K, n_keep = decay_k.shape
    out = dict(het_n_px=K, het_var_obs=np.nan, het_var_null_med=np.nan,
               het_var_ratio=np.nan, het_var_pval=np.nan,
               moran_I_obs=np.nan, moran_I_null_med=np.nan,
               moran_I_z=np.nan, moran_I_pval=np.nan, het_B=B)
    if K < 3:
        return out
    pinvT = np.linalg.pinv(Af).T                       # (n_keep, 4)
    w = n_phot_k.astype(np.float64)
    if w.sum() <= 0:
        return out

    # Observed within-FA tau map (clamp estimator) + its weighted variance.
    tau_obs = _clamp_tau_int(decay_k, pinvT, taus)     # (K,)
    valid = np.isfinite(tau_obs)
    if valid.sum() < 3 or w[valid].sum() <= 0:
        return out
    var_obs = float(_wvar(tau_obs[valid], w[valid]))

    # Uniform-tau model shape over kept bins, scaled to each pixel's photons.
    model = Af @ a_roi                                 # (n_keep,) expected cnts
    M = model.sum()
    if M <= 0:
        return out
    shape = np.clip(model, 0, None) / M                # sums to 1 over keep
    lam = n_phot_k[:, None] * shape[None, :]           # (K, n_keep) expected

    # Monte-Carlo null: B Poisson resamples of the uniform model.
    synth = rng.poisson(lam[None, :, :].repeat(B, axis=0)).astype(np.float64)
    tau_null = _clamp_tau_int(synth, pinvT, taus)      # (B, K)
    # Weighted variance per draw over valid pixels (mask invalid by zero weight).
    wmat = np.where(np.isfinite(tau_null), w[None, :], 0.0)
    tn = np.where(np.isfinite(tau_null), tau_null, 0.0)
    wsum = wmat.sum(axis=1)
    ok = wsum > 0
    var_null = np.full(B, np.nan)
    mean_n = np.zeros(B)
    mean_n[ok] = (wmat[ok] * tn[ok]).sum(1) / wsum[ok]
    var_null[ok] = (wmat[ok] * (tn[ok] - mean_n[ok, None])**2).sum(1) / wsum[ok]

    vn = var_null[np.isfinite(var_null)]
    if vn.size:
        med = float(np.median(vn))
        out["het_var_obs"] = var_obs
        out["het_var_null_med"] = med
        out["het_var_ratio"] = var_obs / med if med > 0 else np.nan
        out["het_var_pval"] = float((1 + (vn >= var_obs).sum()) / (vn.size + 1))

    # Moran's I: observed vs the same MC null (bias-aware).
    ei, ej = _moran_edges(used)
    if ei.size:
        I_obs = float(_moran_I(tau_obs, ei, ej))
        # null: use draws where every pixel produced a finite tau (so adjacency
        # is intact); fall back to filling NaNs with the draw mean otherwise.
        tn_full = np.where(np.isfinite(tau_null), tau_null,
                           np.nan_to_num(np.nanmean(tau_null, axis=1, keepdims=True)))
        I_null = _moran_I(tn_full, ei, ej)
        I_null = I_null[np.isfinite(I_null)]
        if np.isfinite(I_obs) and I_null.size:
            out["moran_I_obs"] = I_obs
            out["moran_I_null_med"] = float(np.median(I_null))
            sd = float(I_null.std())
            out["moran_I_z"] = (I_obs - I_null.mean()) / sd if sd > 0 else np.nan
            out["moran_I_pval"] = float((1 + (I_null >= I_obs).sum()) / (I_null.size + 1))
    return out


# --------------------------------------------------------------------------- #
# Per-image driver
# --------------------------------------------------------------------------- #
def fit_image(cube: np.ndarray, labels: np.ndarray, Af: np.ndarray,
              keep: np.ndarray, taus: np.ndarray, dof: int,
              min_pix: float, min_roi: float, min_area: int,
              B: int, rng) -> tuple[list[dict], dict | None]:
    """Per-FA NNLS + heterogeneity for one image, plus a pooled-in-FA cell row."""
    H, W, _ = cube.shape
    if labels.shape != (H, W):
        raise ValueError(f"label shape {labels.shape} != cube ({H},{W})")
    counts_keep = cube[..., keep].sum(axis=2)               # photons / pixel
    bright = counts_keep >= float(min_pix)

    def _fit_pool(used_mask):
        """Exact NNLS on the pooled decay of `used_mask` pixels."""
        decay = cube[used_mask].sum(axis=0).astype(np.float64)
        n_phot = float(decay[keep].sum())
        x, _ = nnls(Af, decay[keep])
        a = x[:3]
        sat = (a * taus).sum()
        tau_int = float((a * taus**2).sum() / sat) if sat > 0 else float("nan")
        tau_amp = float(sat / a.sum()) if a.sum() > 0 else float("nan")
        model = Af @ x
        chi2 = float(np.sum((decay[keep] - model)**2 /
                            np.clip(model, 1.0, None)) / dof)
        return x, tau_int, tau_amp, n_phot, chi2

    rows = []
    roi_ids = np.unique(labels)
    roi_ids = roi_ids[roi_ids > 0]
    for rid in roi_ids:
        roi_full = labels == rid
        area_px = int(roi_full.sum())
        if area_px < min_area:
            continue
        used = roi_full & bright
        if int(used.sum()) < 3:
            continue
        x, tau_int, tau_amp, n_phot, chi2 = _fit_pool(used)
        if n_phot < min_roi:
            continue
        decay_k = cube[used][:, keep].astype(np.float64)    # (K, n_keep)
        n_phot_k = decay_k.sum(axis=1)
        het = heterogeneity_test(decay_k, n_phot_k, used, Af, x[:4].astype(np.float64),
                                 taus, B, rng)
        row = dict(roi_id=int(rid), area_px=area_px,
                   n_pixels_used=int(used.sum()), n_photons_total=n_phot,
                   mean_intensity=float(counts_keep[roi_full].mean()),
                   tau_int_ns=tau_int, tau_amp_ns=tau_amp,
                   alpha1=float(x[0]), alpha2=float(x[1]), alpha3=float(x[2]),
                   alpha_bg=float(x[3]), chi2_red=chi2)
        row.update(het)
        rows.append(row)

    # Per-cell (per-image) pooled in-FA decay = all labelled pixels.
    cell = None
    in_fa = (labels > 0) & bright
    if int(in_fa.sum()) >= min_area:
        x, tau_int, tau_amp, n_phot, chi2 = _fit_pool(in_fa)
        if n_phot >= min_roi:
            cell = dict(tau_inFA_nnls_ns=tau_int, tau_inFA_nnls_amp_ns=tau_amp,
                        n_inFA_pixels=int(in_fa.sum()), in_FA_photons=n_phot,
                        n_FAs=int(labels.max()), chi2_red=chi2)
    return rows, cell


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spreadsheets", nargs="+", type=Path, required=True,
                    help="LAS X spreadsheet(s) with fitted tau_i")
    ap.add_argument("--ptu-folder", type=Path, required=True,
                    help="folder of source PTU files")
    ap.add_argument("--labels-dir", type=Path, required=True,
                    help="folder of <image_name>_labels.tif FA label maps")
    ap.add_argument("--output-dir", type=Path, required=True,
                    help="where to write the CSVs")
    ap.add_argument("--irf-ptu", type=Path, required=True,
                    help="Rhodamine B reference PTU for IRF extraction; must be "
                         "the same reference the matching flim_pipeline run used")
    ap.add_argument("--include-path-regex", default=None, metavar="RE",
                    help="keep only spreadsheet rows whose Path matches this "
                         "regex. The published runs used it to restrict one "
                         "spreadsheet to one construct/treatment subset (e.g. "
                         "'^185-(dmso|bleb)'); omit to fit every row.")
    ap.add_argument("--tau-r", type=float, default=fp.DEFAULT_TAU_R,
                    help="Rhodamine B reference lifetime (ns); must match --irf-ptu")
    ap.add_argument("--ref-components", nargs="+", default=None, metavar="AMP,TAU",
                    help="the reference's OWN decay model, as AMP,TAU pairs from the LAS X fit of the Rhodamine B reference (e.g. --ref-components 0.0918,0.2868 0.9082,2.1146). Overrides --tau-r for the deconvolution kernel. Use this whenever the reference needs >1 component to fit well: stored dye aliquots in this study need two, and deconvolving them with a mono-exponential kernel distorts the IRF. Omit for a genuinely mono-exponential (freshly prepared) reference.")
    # IRF extraction knobs, forwarded verbatim to flim_pipeline.extract_irf.
    # These MUST match the flim_pipeline run that produced the per-pixel maps and
    # the LAS X validation for this dataset, or the per-FA tau will sit on a
    # different IRF than the QC it was checked against. The defaults here are
    # flim_pipeline's own defaults (damped RL); the 3T3 / MEF / 3T3-distribution
    # datasets in this study were processed with plain RL instead --
    # "--irf-sigma-bins 0 --rl-iters 30 --irf-rl-tol 0". See src/pipeline/README.md.
    ap.add_argument("--rl-iters", type=int, default=400,
                    help="max Richardson-Lucy iterations (default 400)")
    ap.add_argument("--irf-sigma-bins", type=float, default=0.4,
                    help="Gaussian damping scale in time bins for the regularized "
                         "RL IRF extraction; 0 = plain RL (default 0.4)")
    ap.add_argument("--irf-rl-tol", type=float, default=1e-3,
                    help="RL convergence tolerance on the relative L1 step "
                         "(default 1e-3; 0 disables early stopping)")
    ap.add_argument("--irf-extra-shift-ns", type=float, default=0.0,
                    help="extra IRF t0 shift (ns) for this session, added to the "
                         "reference-fitted t0 (default 0)")
    ap.add_argument("--min-photons-per-pixel", type=float, default=1.0,
                    help="skip pixels below N photons (default 1; only drops "
                         "detector-dark pixels -- segmentation defines the FA)")
    ap.add_argument("--min-photons-per-roi", type=float, default=200.0,
                    help="drop FAs whose pooled decay < N photons (default 200)")
    ap.add_argument("--min-area-px", type=int, default=5,
                    help="drop FAs smaller than this (default 5)")
    ap.add_argument("--het-reps", type=int, default=200,
                    help="Monte-Carlo draws for the heterogeneity null (default 200)")
    ap.add_argument("--exclude-ns", nargs="+", default=None, metavar="LO-HI",
                    help="time range(s) in ns to drop from the fit, e.g. '7.0-7.6'")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed for the MC null")
    ap.add_argument("--limit", type=int, default=None, help="process at most N images")
    return ap.parse_args(argv)


def main(argv=None):
    import tifffile
    args = parse_args(argv)
    spreadsheets = [Path(s) for s in args.spreadsheets]
    ptu_folder = Path(args.ptu_folder)
    labels_dir = Path(args.labels_dir)
    output_dir = Path(args.output_dir)
    setup_logging(output_dir)
    log.info("Per-FA FLIM fit -> %s", output_dir)
    log.info("  spreadsheets=%s", [p.name for p in spreadsheets])
    log.info("  ptu-folder=%s  labels-dir=%s", ptu_folder, labels_dir)
    log.info("  floors: pixel>=%g, roi>=%g photons, area>=%d px; het B=%d",
             args.min_photons_per_pixel, args.min_photons_per_roi,
             args.min_area_px, args.het_reps)

    # Re-extract the IRF (fast) so we reuse the exact t0 shift flim_pipeline used.
    log.info("  IRF: rl_iters=%d sigma_bins=%.2f rl_tol=%.0e extra_shift=%+.4f ns",
             args.rl_iters, args.irf_sigma_bins, args.irf_rl_tol,
             args.irf_extra_shift_ns)
    irf_info = fp.extract_irf(args.irf_ptu, args.tau_r, output_dir / "irf_per_roi",
                              rl_iters=args.rl_iters,
                              sigma_bins=args.irf_sigma_bins,
                              rl_tol=args.irf_rl_tol,
                              extra_shift_ns=args.irf_extra_shift_ns,
                              ref_components=fp.parse_ref_components(args.ref_components))

    df = fp.load_spreadsheets(spreadsheets)
    # Optional row filter, so one LAS X spreadsheet covering several constructs
    # or treatments can be fitted one subset at a time. Out-of-scope rows have
    # no labels anyway, but filtering keeps the log clean and avoids re-reading
    # their PTUs.
    if args.include_path_regex:
        import re
        pat = re.compile(args.include_path_regex)
        n_before = len(df)
        df = df[df[fp.COL_PATH].astype(str).map(
            lambda s: pat.search(s) is not None)].reset_index(drop=True)
        log.info("--include-path-regex %r: kept %d/%d spreadsheet rows",
                 args.include_path_regex, len(df), n_before)
    if args.limit:
        df = df.iloc[:args.limit]
    log.info("Processing %d images", len(df))

    exclude_ranges = []
    for tok in (args.exclude_ns or []):
        lo, hi = (float(v) for v in tok.split("-"))
        exclude_ranges.append((lo, hi))

    rng = np.random.default_rng(args.seed)
    out_rows, cell_rows, img_rows = [], [], []
    for _, row in df.iterrows():
        name = str(row[fp.COL_PATH])
        ptu = fp.find_ptu_for_path(ptu_folder, name)
        if ptu is None:
            log.warning("  no PTU for '%s'; skipping", name)
            continue
        image_name = ptu.stem
        lp = labels_dir / f"{image_name}_labels.tif"
        if not lp.exists():
            log.warning("  no labels for %s; skipping", image_name)
            continue
        labels = tifffile.imread(lp).astype(np.int32)
        if labels.max() == 0:
            img_rows.append(dict(image_name=image_name, n_rois=0, n_rois_fitted=0))
            continue
        cube, n_period, dt_ns, _ = fp.read_ptu_cube(ptu)
        n_use = min(n_period, irf_info["n"])
        cube = cube[..., :n_use].astype(np.float64)
        irf_bg = float(row.get(fp.COL_IRF_BG, 0.0) or 0.0)
        if irf_bg:
            cube = np.clip(cube - irf_bg, 0, None)
        taus = np.array([float(row[c]) for c in fp.COL_TAU])
        irf_shift = float(row.get(fp.COL_IRF_SHIFT, 0.0) or 0.0)
        A = fp.build_basis(irf_info, taus, irf_shift)[:n_use]
        keep = fp.build_keep_mask(irf_info["t_ns"][:n_use], exclude_ranges)
        Af = A[keep]
        dof = max(int(keep.sum()) - 4, 1)

        rois, cell = fit_image(cube, labels, Af, keep, taus, dof,
                               args.min_photons_per_pixel, args.min_photons_per_roi,
                               args.min_area_px, args.het_reps, rng)
        lasx = float(row.get(fp.COL_MEAN_INT_LASX, np.nan))
        for r in rois:
            r.update(image_name=image_name, path=name,
                     tau1_ns=float(taus[0]), tau2_ns=float(taus[1]),
                     tau3_ns=float(taus[2]), lasx_mean_tau_int_ns=lasx)
        out_rows.extend(rois)
        if cell is not None:
            cell.update(image_name=image_name, path=name, lasx_mean_tau_int_ns=lasx)
            cell_rows.append(cell)
        img_rows.append(dict(image_name=image_name, n_rois=int(labels.max()),
                             n_rois_fitted=len(rois)))
        log.info("  %-22s FAs=%4d fitted=%4d", image_name, int(labels.max()), len(rois))

    output_dir.mkdir(parents=True, exist_ok=True)
    if out_rows:
        cols = ["image_name", "path", "roi_id", "area_px", "n_pixels_used",
                "n_photons_total", "mean_intensity", "tau_int_ns", "tau_amp_ns",
                "alpha1", "alpha2", "alpha3", "alpha_bg", "chi2_red",
                "het_n_px", "het_var_obs", "het_var_null_med", "het_var_ratio",
                "het_var_pval", "moran_I_obs", "moran_I_null_med", "moran_I_z",
                "moran_I_pval", "het_B",
                "tau1_ns", "tau2_ns", "tau3_ns", "lasx_mean_tau_int_ns"]
        odf = pd.DataFrame(out_rows)
        odf = odf[[c for c in cols if c in odf.columns]]
        p = output_dir / "per_roi_measurements.csv"
        odf.to_csv(p, index=False, float_format="%.5f")
        log.info("Wrote %d FA rows -> %s", len(odf), p)
        sig = (odf["het_var_pval"] < 0.05).mean()
        log.info("  median het_var_ratio=%.2f; %.0f%% of FAs p<0.05 (excess variance)",
                 odf["het_var_ratio"].median(), 100 * sig)
    else:
        log.warning("No FA rows produced.")

    if cell_rows:
        p = output_dir / "per_cell_roi.csv"
        pd.DataFrame(cell_rows).to_csv(p, index=False, float_format="%.5f")
        log.info("Wrote %d per-cell rows -> %s", len(cell_rows), p)

    pd.DataFrame(img_rows).to_csv(output_dir / "per_roi_per_image_summary.csv",
                                  index=False)
    log.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
