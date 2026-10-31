#!/usr/bin/env python3
"""Check the per-ROI FLIM fit against a committed slice of the real data.

What this exercises
-------------------
The pooled-decay NNLS fit that produces every in-ROI lifetime in the paper:
`src/pipeline/flim_pipeline.build_basis` (the reconvolution model basis) followed
by `src/pipeline/flim_fa_fit.fit_image` (exact NNLS on each ROI's summed raw
decay, plus the per-image pooled in-ROI row). Same functions, same inputs, no
re-implementation.

The example input
-----------------
`input/185-bleb-1_crop.npz` is a 192x192 pixel crop of the real photon cube of
image 185-bleb-1 (Vin-WF, blebbistatin 20 uM, dataset dmso-bleb20-0727-185, which
feeds Fig. 2e), together with the Cellpose focal-adhesion label crop, the IRF
vector, and the LAS X tau_1..tau_3 for that image. The crop window was chosen so
that all 45 focal adhesions in it are *fully contained* -- none is clipped by the
crop boundary. Because each ROI's fit only ever sums that ROI's own pixels, the
fit of the crop must reproduce the full-image fit exactly, so
`expected/per_roi_measurements.csv` is a verbatim 45-row slice of the published
`per_roi_measurements.csv` for this image -- not a value regenerated for this
example.

The full 512x512 PTU (`185-bleb-1.ptu`, ~14 MB) and the LAS X spreadsheet are in
the raw-data deposit (docs/RAW_DATA.md) if you want to run the whole entry point;
the committed .npz is the same photon data, cropped, so no PTU reader is needed
here.

Tolerances
----------
The reference CSV was produced on Windows and is rounded to %.5f, so exact
equality is not the right test. Photon counts and pixel areas are checked exactly;
the fitted lifetimes and the observed within-ROI statistics are checked to that
1e-5 rounding; the amplitudes get a relative tolerance because their split is
mildly ill-conditioned. The columns compared against the shot-noise Monte-Carlo
null are printed but not gated -- their RNG stream depends on which ROIs the run
covers, so a 45-of-131 crop cannot replay it. No figure panel plots those columns.
See src/pipeline/README.md for the measured agreement.

Usage:
    python examples/flim_roi_fit/check_flim_roi_fit.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "src" / "pipeline"))

# Columns that are a deterministic function of the photons and the basis, with
# the absolute tolerance each is checked to. 1e-5 is the reference CSV's own
# %.5f rounding; the amplitudes are O(100) with a mildly ill-conditioned split,
# so they get a relative tolerance instead (see rtol below).
EXACT_COLS = ["area_px", "n_pixels_used", "n_photons_total", "het_n_px"]
# het_var_obs and moran_I_obs are the *observed* within-ROI statistics: computed
# straight from the photons, no resampling, so they are checked like the fit. Only
# the columns compared against the Monte-Carlo null are non-deterministic.
FIT_COLS = ["tau_int_ns", "tau_amp_ns", "chi2_red", "mean_intensity",
            "het_var_obs", "moran_I_obs"]
ALPHA_COLS = ["alpha1", "alpha2", "alpha3", "alpha_bg"]
FIT_ATOL = 1e-5
# The amplitudes are checked on the np.isclose combination: alpha_bg is a small
# residual-background term that sits near zero, where a pure relative tolerance
# would be meaningless, while alpha1..alpha3 are O(100) where rounding is.
ALPHA_RTOL = 1e-3
ALPHA_ATOL = 1e-4
# Columns derived from the shot-noise Monte-Carlo null: reported, not gated. The
# published run drew from a single seeded stream advanced across every ROI of
# every image in the dataset, so a re-run over a different ROI subset -- such as
# this 45-of-131 crop -- draws different numbers by construction.
MC_COLS = ["het_var_null_med", "het_var_ratio", "het_var_pval",
           "moran_I_null_med", "moran_I_z", "moran_I_pval"]


def main() -> int:
    try:
        import pandas as pd
        import flim_pipeline as fp
        import flim_fa_fit as fa
    except ImportError as e:
        print(f"[SKIP] {e}. This check needs numpy, scipy and pandas "
              f"(pip install -r requirements.txt).")
        return 0

    npz_path = HERE / "input" / "185-bleb-1_crop.npz"
    exp_path = HERE / "expected" / "per_roi_measurements.csv"
    z = np.load(npz_path, allow_pickle=False)
    cube = z["photon_cube"].astype(np.float64)
    labels = z["labels"].astype(np.int32)
    taus = z["taus"]
    irf_bg = float(z["irf_bg_cnts"])
    if irf_bg:
        cube = np.clip(cube - irf_bg, 0, None)

    # Rebuild the model basis exactly as flim_pipeline/flim_fa_fit do.
    irf_info = {"irf": z["irf"], "t_ns": z["t_ns"], "dt_ns": float(z["dt_ns"]),
                "n": int(z["irf"].size), "t0_shift_ns": float(z["t0_shift_ns"])}
    A = fp.build_basis(irf_info, taus, float(z["irf_shift_ns"]))
    keep = fp.build_keep_mask(irf_info["t_ns"], [])
    Af, dof = A[keep], max(int(keep.sum()) - 4, 1)

    # Same floors as the shipped run for this dataset (flim_fa_fit defaults).
    # The IRF in the .npz already carries that run's reference model, damping and
    # session shift, so no IRF settings are needed here.
    rows, fov = fa.fit_image(cube, labels, Af, keep, taus, dof,
                             min_pix=1.0, min_roi=200.0, min_area=5,
                             B=200, rng=np.random.default_rng(0))
    got = pd.DataFrame(rows).sort_values("roi_id").reset_index(drop=True)
    exp = pd.read_csv(exp_path).sort_values("roi_id").reset_index(drop=True)

    print(f"input    {npz_path.relative_to(REPO)}  "
          f"({cube.shape[0]}x{cube.shape[1]} px, {cube.shape[2]} time bins)")
    print(f"expected {exp_path.relative_to(REPO)}  ({len(exp)} ROIs)")

    failures = []
    if len(got) != len(exp) or not np.array_equal(got.roi_id, exp.roi_id):
        print(f"[FAIL] fitted {len(got)} ROIs, expected {len(exp)}; "
              f"roi_id sets differ")
        return 1
    print(f"fitted   {len(got)} ROIs, roi_ids match")

    for c in EXACT_COLS:
        d = np.abs(got[c].to_numpy(float) - exp[c].to_numpy(float)).max()
        ok = d == 0
        print(f"  [{'OK ' if ok else 'FAIL'}] {c:16s} exact              "
              f"max|diff| = {d:.3g}")
        if not ok:
            failures.append(c)
    for c in FIT_COLS:
        d = np.abs(got[c].to_numpy(float) - exp[c].to_numpy(float)).max()
        ok = d <= FIT_ATOL
        print(f"  [{'OK ' if ok else 'FAIL'}] {c:16s} atol {FIT_ATOL:<9.0e} "
              f"max|diff| = {d:.3g}")
        if not ok:
            failures.append(c)
    for c in ALPHA_COLS:
        a, b = got[c].to_numpy(float), exp[c].to_numpy(float)
        ok = bool(np.allclose(a, b, rtol=ALPHA_RTOL, atol=ALPHA_ATOL))
        excess = np.abs(a - b) - (ALPHA_ATOL + ALPHA_RTOL * np.abs(b))
        print(f"  [{'OK ' if ok else 'FAIL'}] {c:16s} "
              f"rtol {ALPHA_RTOL:.0e}/atol {ALPHA_ATOL:.0e}  "
              f"worst margin = {-excess.min():.3g}")
        if not ok:
            failures.append(c)
    print("  ---- Monte-Carlo columns (reported, not gated; see docstring) ----")
    for c in MC_COLS:
        d = np.abs(got[c].to_numpy(float) - exp[c].to_numpy(float)).max()
        print(f"  [   ] {c:16s} informational      max|diff| = {d:.3g}")

    # The per-image pooled in-FA lifetime: the quantity Fig. 2c and 2e plot.
    # In the published run it is pooled over the whole 512x512 image, so it is
    # not comparable to this crop; check it is finite and physically sensible.
    if fov is None:
        print("  [FAIL] per-FOV pooled in-ROI fit returned nothing")
        failures.append("per_FOV")
    else:
        t = fov["tau_inFA_nnls_ns"]
        ok = np.isfinite(t) and 1.0 < t < 4.0
        print(f"  [{'OK ' if ok else 'FAIL'}] per-FOV pooled in-ROI tau_int = "
              f"{t:.5f} ns over {fov['n_inFA_pixels']} px "
              f"(crop only; the published per-FOV value pools the full image)")
        if not ok:
            failures.append("per_FOV")

    if failures:
        print(f"\nFAILED: {len(failures)} column(s): {failures}")
        return 1
    print("\nOK: the per-ROI FLIM fit reproduces the published values.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
