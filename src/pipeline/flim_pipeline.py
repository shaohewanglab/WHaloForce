#!/usr/bin/env python3
"""
flim_pipeline.py — Hybrid FLIM batch analysis.

Workflow
--------
LAS X performs the global 3-exponential reconvolution fit and exports, per
image, the three lifetimes tau_1..tau_3, the per-image IRF shift, and global
amplitudes/mean lifetimes (the validation reference). This pipeline:

  STAGE 1  Extract the instrument response function (IRF) from a Rhodamine B
           reference PTU by iterative (Richardson-Lucy) deconvolution, assuming
           the reference decay is a single exponential with lifetime tau_R.
           Also fit one global time-zero (t0) shift. Saves irf_vector.npy and
           irf_diagnostics.png for inspection.

  STAGE 2  For every experimental PTU, hold tau_1..tau_3 fixed (from the LAS X
           spreadsheet) and solve a per-pixel linear amplitude fit of the model
              data(t) ~ sum_i alpha_i * [IRF_shifted (*) exp(-t/tau_i)] + bg
           over a vectorised (batched) least-squares solve. Produces
           intensity-weighted mean-lifetime maps and QC images.

  STAGE 3  Aggregate per-image mean lifetimes and validate them against the LAS
           X values; warn on any image that disagrees.

Modelling notes
---------------
* The TCSPC histograms are periodic with the laser repetition period
  (number_bins_in_period bins; ~12.8 ns here). The recorded histograms contain
  a few extra bins beyond one period -- those are used only to estimate the
  (small) dark-count floor. All modelling is done on exactly one period using
  *circular* convolution, which automatically reproduces the wrapped
  incomplete-decay baseline seen before the rising edge. A per-pixel constant
  ("ones") column then absorbs residual dark counts / offset.
* The rising edge of the Rhodamine reference and of the experimental decays
  coincide (same system t0), which is consistent with the spreadsheet's
  per-image "IRF Shift ns" being 0 -- so the reference-derived IRF is already
  time-aligned to the experimental data. The fitted global t0 shift is added on
  top of the per-image spreadsheet shift when building each model basis.

Amplitude convention
--------------------
The exponential basis functions exp(-t/tau_i) are normalised to 1 at t=0, so the
fitted alpha_i are pre-exponential amplitudes (matching LAS X "Amplitude N").
With that convention:
    intensity_i  ~ alpha_i * tau_i
    tau_amp = sum(alpha_i*tau_i) / sum(alpha_i)              (amplitude weighted)
    tau_int = sum(alpha_i*tau_i^2) / sum(alpha_i*tau_i)      (intensity weighted)
Both reproduce the LAS X spreadsheet values when fed the spreadsheet amplitudes.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

# ptufile, lmfit, tifffile, matplotlib and FLIM_vis_v15 are imported lazily
# inside the functions that need them (PTU reading, IRF t0 fitting, and writing
# maps / QC images respectively). That keeps the pure-numpy parts of this module
# -- build_basis, build_keep_mask, _weighted_means -- importable with only
# numpy/scipy/pandas installed, so examples/flim_roi_fit/ can check the fit core
# against committed inputs without the full FLIM environment.

log = logging.getLogger("flim")

# --------------------------------------------------------------------------- #
# Defaults
# --------------------------------------------------------------------------- #
HERE = Path(__file__).resolve().parent
# The inputs (LAS X spreadsheet, PTU folder, Rhodamine B reference PTU) live in
# the raw-data deposit, not in this repository, so they have no in-repo default
# and must be given on the command line. See docs/RAW_DATA.md and the per-
# dataset invocations in src/pipeline/README.md.
#
# tau_R is the reference lifetime of the Rhodamine B IRF standard and depends on
# the detector and spectral window of the reference acquisition, so it must be
# matched to --irf-ptu. Values used in this study: 2.088 ns (HyD X 2, 561,
# 565-650 nm; the 3T3/MEF reference), 2.196 ns (HyD X 4, 553, 559-740 nm),
# 2.192 ns (HyD X 4, 553, 559-700 nm).
DEFAULT_TAU_R = 2.196

# Exact spreadsheet column names (note the DOUBLE space before "ns" in the
# intensity-weighted column -- match the file verbatim).
COL_PATH = "Path"
COL_TAU = ["Lifetime (τ) 1 ns", "Lifetime (τ) 2 ns", "Lifetime (τ) 3 ns"]
COL_IRF_SHIFT = "IRF Shift ns"
COL_IRF_BG = "IRF Background Cnts"
COL_MEAN_INT_LASX = "Mean τ, Intensity Weighted  ns"
COL_MEAN_AMP_LASX = "Mean τ, Amplitude Weighted ns"
COL_REGION = "Region"   # time-lapse exports only: "Overall Decay" vs "ROI N"
COL_TIME = "Time"       # time-lapse exports only: "<seconds> s" string
T_INDEX_COL = "__t_index__"  # synthetic column we set during loading


# --------------------------------------------------------------------------- #
# Small numerical helpers
# --------------------------------------------------------------------------- #
def cconv(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Circular (periodic) convolution of two equal-length vectors."""
    return np.real(np.fft.ifft(np.fft.fft(a) * np.fft.fft(b)))


def circular_shift(vec: np.ndarray, shift_bins: float) -> np.ndarray:
    """Shift `vec` by `shift_bins` (positive = later in time) with sub-bin
    accuracy using an FFT phase ramp. Periodic, so it suits the IRF on one
    laser period."""
    n = vec.size
    freqs = np.fft.fftfreq(n)
    phase = np.exp(-2j * np.pi * freqs * shift_bins)
    return np.real(np.fft.ifft(np.fft.fft(vec) * phase))


def mono_exp_kernel(t_ns: np.ndarray, tau_ns: float) -> np.ndarray:
    """Sum-normalised periodic mono-exponential decay kernel on the time grid.

    This is the measured decay produced by a delta-function IRF: a single
    exponential made periodic by the circular convolution used elsewhere.
    """
    k = np.exp(-t_ns / float(tau_ns))
    return k / k.sum()


def multi_exp_kernel(t_ns: np.ndarray, components) -> np.ndarray:
    """Sum-normalised periodic multi-exponential reference-decay kernel.

    `components` is a sequence of (amplitude, tau_ns) pairs -- the LAS X fit of
    the Rhodamine B reference itself. Amplitudes are pre-exponential (the same
    convention as the experimental fit) and are normalised here, so only their
    ratio matters.

    Why this exists: the IRF is recovered by deconvolving the reference decay out
    of the measured reference, which requires a model of that decay. A single
    exponential is only right if the reference really is mono-exponential.
    Rhodamine B aliquots stored at -20 C in this study developed a small fast
    component (~6-11% amplitude, ~0.3 ns): fitting those references with one
    component gives LAS X chi2 ~50, with two components ~1.9. Deconvolving with
    a mono-exponential kernel when the truth is bi-exponential leaves a distorted
    IRF, which then shows up as an oscillating residual at the rising edge of
    every experimental fit. Freshly prepared reference dye (the 2026-01-27 set
    here) is mono-exponential and needs only `mono_exp_kernel`.
    """
    k = np.zeros_like(np.asarray(t_ns, dtype=float))
    for amp, tau in components:
        k += float(amp) * np.exp(-t_ns / float(tau))
    return k / k.sum()


def parse_ref_components(tokens):
    """Parse ['0.0918,0.2868', '0.9082,2.1146'] -> [(0.0918, 0.2868), ...]."""
    if not tokens:
        return None
    out = []
    for tok in tokens:
        amp, _, tau = str(tok).partition(",")
        if not tau:
            raise ValueError(f"--ref-components expects AMP,TAU pairs; got {tok!r}")
        out.append((float(amp), float(tau)))
    return out


def fwhm_ns(vec: np.ndarray, dt_ns: float) -> float:
    """Full width at half maximum of a peaked vector, in ns."""
    if vec.max() <= 0:
        return float("nan")
    above = np.where(vec >= vec.max() / 2.0)[0]
    if above.size == 0:
        return float("nan")
    return (above[-1] - above[0] + 1) * dt_ns


# --------------------------------------------------------------------------- #
# PTU reading
# --------------------------------------------------------------------------- #
def read_ptu_cube(path: Path):
    """Read a PTU and return (cube[Y,X,H], n_period, dt_ns, dark_floor).

    cube is the photon histogram summed over any T frames and detector channels.
    dark_floor is estimated from histogram bins *beyond* one laser period
    (laser-off region), defaulting to 0 when no such bins exist.
    """
    import ptufile
    p = ptufile.PtuFile(str(path))
    n = int(p.number_bins_in_period)
    dt_ns = float(p.tcspc_resolution) * 1e9
    arr = p[:]  # (T, Y, X, C, H)
    full = np.asarray(arr).sum(axis=(0, 3))  # -> (Y, X, H)
    H = full.shape[-1]
    if H > n + 4:
        dark = float(np.median(full[..., n + 4:].sum(axis=(0, 1))) / (full.shape[0] * full.shape[1]))
        # dark per pixel is tiny; keep the summed-decay version separately too
    dark_decay = full[..., n:].sum(axis=(0, 1)).astype(float)
    dark_floor = float(np.median(dark_decay)) if dark_decay.size else 0.0
    cube = full[..., :n].astype(np.float64)
    p.close()
    return cube, n, dt_ns, dark_floor


def read_ptu_summed_decay(path: Path):
    """Read a PTU and return (decay[H], n_period, dt_ns) summed over all pixels
    and channels -- used for the high-SNR reference."""
    import ptufile
    p = ptufile.PtuFile(str(path))
    n = int(p.number_bins_in_period)
    dt_ns = float(p.tcspc_resolution) * 1e9
    decay = p.decode_histogram(dtype="uint64").sum(axis=0).astype(np.float64)
    p.close()
    return decay, n, dt_ns


def read_pixel_um(path: Path) -> float:
    """Pixel size (micron/px) from a PTU's ImgHdr_PixResol tag; NaN if absent.

    Used to write a per-image pixel_um sidecar so the downstream measurement
    reports ROI area in um^2 correctly even when images were acquired at
    different zoom (e.g. the MEF FN-vs-uncoat set mixes 0.190 and 0.284 um/px).
    """
    import ptufile
    p = ptufile.PtuFile(str(path))
    pix = p.tags.get("ImgHdr_PixResol")
    p.close()
    return float(pix) if pix is not None else float("nan")


# --------------------------------------------------------------------------- #
# STAGE 1 -- IRF extraction
# --------------------------------------------------------------------------- #
def extract_irf(irf_ptu: Path, tau_r: float, output_dir: Path,
                rl_iters: int = 400, sigma_bins: float = 0.4,
                rl_tol: float = 1e-3, rl_min_iters: int = 5,
                extra_shift_ns: float = 0.0, ref_components=None):
    """Extract the IRF from the Rhodamine B reference by *regularized* (damped)
    Richardson-Lucy deconvolution, stopped on convergence, and fit a global t0
    shift.

    Plain RL is semi-convergent on these high-SNR references: the reconvolution
    residual never plateaus, so a fixed iteration count over-sharpens the IRF
    toward a 1-bin delta -- and how far it collapses depends on the reference's
    SNR (the new 512x512 RhB refs collapse to ~2 bins / 194 ps by 30 iters,
    biasing every downstream mean_tau_int ~+0.06 ns high vs LAS X). To make the
    recovered width stable and SNR-independent we lightly Gaussian-smooth the
    IRF estimate each iteration (`sigma_bins`). The smoothing balances RL
    sharpening at a fixed resolution scale, so the iteration converges to a
    stable shape (||delta irf||_1 / ||irf||_1 < `rl_tol`) instead of collapsing.
    sigma_bins=0.4 gives FWHM ~388 ps, which minimizes the experimental global
    reduced chi2 across the 20250715 sessions (and matches the Leica IRF
    screenshot); it reads ~+0.09 ns vs LAS X mean_tau_int, a systematic that
    cancels in per-cell delta-tau. (sigma_bins=0.6 / FWHM ~485 ps instead
    centres mean_tau_int on LAS X but fits the decays ~3x worse.) Set
    sigma_bins=0 to recover plain RL.

    `extra_shift_ns` is added to the reference-fitted t0 shift: the IRF is
    extracted/aligned against the (2026) RhB reference, but a given experimental
    *session* may have a different laser-sync time-zero. Refining one extra time
    shift against the experimental decays (the analog of LAS X's per-image "IRF
    Shift") corrects it: 20250715 sessions want ~-0.04 ns, 20250630 wants
    ~-0.14 ns (a ~0.10 ns earlier session offset; without it the 20250630 fit
    is ~10x worse and *looks* like it wants a broader IRF).

    Returns a dict with keys: irf, t_ns, dt_ns, n, tau_r, t0_shift_ns, dark,
    chi2, fwhm_ns, peak_bin, sigma_bins, stop_iter, extra_shift_ns.
    """
    log.info("STAGE 1: extracting IRF from %s (tau_R = %.3f ns, sigma_bins = %.2f)",
             irf_ptu.name, tau_r, sigma_bins)
    decay_full, n, dt_ns = read_ptu_summed_decay(irf_ptu)
    t = np.arange(n) * dt_ns
    D = decay_full[:n]

    # Dark-count floor from the laser-off bins beyond one period.
    dark = float(np.median(decay_full[n + 4:])) if decay_full.size > n + 4 else 0.0
    Dsub = np.clip(D - dark, 1e-9, None)
    log.info("  period = %d bins, dt = %.4f ns, dark floor = %.1f cnts", n, dt_ns, dark)

    # Regularized Richardson-Lucy: Dsub = IRF (*) k  (circular), with a circular
    # Gaussian damping each step and convergence-based early stopping.
    # Reference-decay kernel. `ref_components` (the LAS X multi-exponential fit
    # of the reference) takes precedence; a bare tau_r means mono-exponential,
    # which is only valid for a reference that fits well with one component.
    if ref_components:
        k = multi_exp_kernel(t, ref_components)
        log.info("  reference model: %d components %s",
                 len(ref_components),
                 ", ".join(f"a={a:.4g}/tau={x:.4g}ns" for a, x in ref_components))
    else:
        k = mono_exp_kernel(t, tau_r)
        log.info("  reference model: mono-exponential, tau_R = %.4f ns", tau_r)
    kflip = np.roll(k[::-1], 1)  # matched filter for circular correlation
    irf = np.full(n, 1.0 / n)
    stop_iter = rl_iters
    for it in range(1, rl_iters + 1):
        est = cconv(irf, k)
        est[est <= 0] = 1e-12
        new = irf * cconv(Dsub / est, kflip)
        if sigma_bins > 0:
            new = gaussian_filter1d(new, sigma_bins, mode="wrap")
        new = np.clip(new, 0.0, None)
        new /= new.sum()
        rel = float(np.abs(new - irf).sum() / max(np.abs(irf).sum(), 1e-12))
        irf = new
        if it >= rl_min_iters and rel < rl_tol:
            stop_iter = it
            break
    if stop_iter < rl_iters:
        log.info("  regularized RL converged at iter %d (rel step < %.0e)", stop_iter, rl_tol)
    else:
        log.warning("  regularized RL hit max %d iters without reaching rel step %.0e "
                    "(IRF may be under-smoothed; raise sigma_bins or rl_iters)", rl_iters, rl_tol)

    # Refine a single global t0 shift (sub-bin) with a reconvolution fit;
    # tau_R is held fixed per spec. Expected to be small since edges align.
    def residual(params):
        s = params["shift_ns"].value / dt_ns
        amp = params["amp"].value
        bg = params["bg"].value
        model = amp * cconv(circular_shift(irf, s), k) + bg
        model = np.clip(model, 1e-9, None)
        return (D - model) / np.sqrt(model)  # Poisson-weighted

    from lmfit import Minimizer, Parameters
    params = Parameters()
    params.add("amp", value=Dsub.sum(), min=0)
    params.add("shift_ns", value=0.0, min=-1.0, max=1.0)
    params.add("bg", value=dark, min=0)
    result = Minimizer(residual, params).minimize(method="leastsq")
    t0_shift_ref = float(result.params["shift_ns"].value)
    amp_fit = float(result.params["amp"].value)
    bg_fit = float(result.params["bg"].value)

    # Fold in the empirical per-session shift (vs the experimental data).
    t0_shift = t0_shift_ref + float(extra_shift_ns)
    if extra_shift_ns:
        log.info("  t0 shift = %.4f ns (ref) %+.4f ns (session) = %.4f ns total",
                 t0_shift_ref, extra_shift_ns, t0_shift)

    fit = amp_fit * cconv(circular_shift(irf, t0_shift_ref / dt_ns), k) + bg_fit
    fit = np.clip(fit, 1e-9, None)
    chi2 = float(np.sum((D - fit) ** 2 / fit) / (n - len(result.params)))

    peak_bin = int(irf.argmax())
    fwhm = fwhm_ns(irf, dt_ns)
    log.info("  IRF peak at bin %d (%.3f ns), FWHM ~ %.0f ps, tail frac = %.3f",
             peak_bin, peak_bin * dt_ns, fwhm * 1000, float(irf[peak_bin + 5:].sum()))
    log.info("  global t0 shift = %.4f ns, reconvolution reduced chi2 = %.1f",
             t0_shift, chi2)
    if fwhm * 1000 < 250 or fwhm * 1000 > 800:
        log.warning("  IRF FWHM %.0f ps is outside the expected ~300-700 ps window "
                    "for damped RL (check tau_R / sigma_bins=%.2f).", fwhm * 1000, sigma_bins)

    # Save artifacts.
    import json
    output_dir.mkdir(parents=True, exist_ok=True)
    np.save(output_dir / "irf_vector.npy", irf.astype(np.float64))
    log.info("  saved %s", output_dir / "irf_vector.npy")
    irf_meta = {"t0_shift_ns": float(t0_shift), "dt_ns": float(dt_ns),
                "n": int(n), "tau_r_ns": float(tau_r),
                "ref_components": ([[float(a), float(x)] for a, x in ref_components]
                                   if ref_components else None),
                "fwhm_ns": float(fwhm), "dark": float(dark),
                "chi2": float(chi2), "peak_bin": int(peak_bin),
                "sigma_bins": float(sigma_bins), "stop_iter": int(stop_iter),
                "t0_shift_ref_ns": float(t0_shift_ref), "extra_shift_ns": float(extra_shift_ns)}
    with open(output_dir / "irf_meta.json", "w") as f:
        json.dump(irf_meta, f, indent=2)

    _plot_irf_diagnostics(output_dir / "irf_diagnostics.png", t, D, fit, irf,
                          dt_ns, tau_r, t0_shift, chi2, fwhm, dark)
    log.info("  saved %s", output_dir / "irf_diagnostics.png")

    return dict(irf=irf, t_ns=t, dt_ns=dt_ns, n=n, tau_r=tau_r,
                ref_components=ref_components,
                t0_shift_ns=t0_shift, dark=dark, chi2=chi2,
                fwhm_ns=fwhm, peak_bin=peak_bin,
                sigma_bins=sigma_bins, stop_iter=stop_iter,
                extra_shift_ns=float(extra_shift_ns))


def _plot_irf_diagnostics(png_path, t, D, fit, irf, dt_ns, tau_r,
                          t0_shift, chi2, fwhm, dark):
    resid = (D - fit) / np.sqrt(np.clip(fit, 1, None))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 1, figsize=(8, 9), height_ratios=[3, 2, 1.4],
                           sharex=True)

    ax[0].semilogy(t, np.clip(D, 1, None), "k.", ms=3, label="Rhodamine B decay")
    ax[0].semilogy(t, np.clip(fit, 1, None), "r-", lw=1.5,
                   label=f"reconv. fit (τ_R={tau_r:.2f} ns)")
    ax[0].axhline(max(dark, 1), color="gray", ls=":", lw=1, label=f"dark ≈ {dark:.0f}")
    ax[0].set_ylabel("counts (log)")
    ax[0].set_title(f"IRF extraction — reduced χ²={chi2:.0f}, t0 shift={t0_shift:+.3f} ns")
    ax[0].legend(loc="upper right", fontsize=8)

    ax[1].plot(t, irf / irf.max(), "b-", lw=1.5)
    ax[1].fill_between(t, irf / irf.max(), color="b", alpha=0.2)
    ax[1].set_ylabel("recovered IRF\n(norm.)")
    ax[1].set_title(f"IRF: peak {irf.argmax() * dt_ns:.3f} ns, "
                    f"FWHM ≈ {fwhm * 1000:.0f} ps", fontsize=10)

    ax[2].plot(t, resid, "g-", lw=0.8)
    ax[2].axhline(0, color="k", lw=0.6)
    ax[2].set_ylabel("weighted\nresidual")
    ax[2].set_xlabel("time (ns)")
    fig.tight_layout()
    fig.savefig(png_path, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Spreadsheets
# --------------------------------------------------------------------------- #
def load_spreadsheets(paths) -> pd.DataFrame:
    """Read and concatenate the LAS X global-fit spreadsheets.

    Time-lapse exports include extra rows for user-drawn ROIs; only the
    "Overall Decay" rows correspond to one-per-PTU global fits. If a
    `Region` column is present we drop everything else automatically. If a
    `Time` column is present we assign each row a synthetic 1-based t-index
    by gap-aware rank within its Path group (`assign_t_index`), and the
    PTU file is then `<Path>_t<idx>.ptu` instead of `<Path>.ptu`.
    """
    frames = []
    for p in paths:
        p = Path(p)
        if not p.exists():
            log.warning("spreadsheet not found, skipping: %s", p)
            continue
        df = pd.read_excel(p)
        df["__source__"] = p.name
        frames.append(df)
        log.info("  loaded %s (%d rows)", p.name, len(df))
    if not frames:
        raise SystemExit("No spreadsheets could be read.")
    out = pd.concat(frames, ignore_index=True)
    out[COL_PATH] = out[COL_PATH].astype(str).str.strip()
    if COL_REGION in out.columns:
        before = len(out)
        out = out[out[COL_REGION].astype(str).str.strip() == "Overall Decay"].copy()
        log.info("  filtered to Region=='Overall Decay': %d -> %d rows", before, len(out))
    if COL_TIME in out.columns:
        out = assign_t_index(out)
    return out


def assign_t_index(df: pd.DataFrame) -> pd.DataFrame:
    """Add a 1-based t-index per row, grouped by Path.

    LAS X's Time column has a known bug where the first reported Time is
    sometimes anomalously short (e.g. 0 then 11s instead of 0 then 60s with
    a true 60s/frame cadence). Rather than scaling row times against an
    absolute median dt, we use *row order within Path* as the default
    t-index (row 1 -> t1, row 2 -> t2, ...). When a *large* inter-row gap
    is detected (gap > 1.5 x median dt computed from interior intervals
    only, which are bug-free), we bump the t-index of subsequent rows by
    the inferred number of missing frames -- e.g. one omitted frame in the
    middle of Position 1 leaves t89 unused so the last 31 rows shift to
    t90..t120.
    """
    df = df.copy()
    t_s = df[COL_TIME].astype(str).str.replace("s", "", regex=False).str.strip()
    df["__time_s__"] = pd.to_numeric(t_s, errors="coerce")
    df[T_INDEX_COL] = 0
    for path, grp in df.groupby(COL_PATH, sort=False):
        g = grp.sort_values("__time_s__")
        times = g["__time_s__"].values.astype(float)
        n = len(times)
        # Median dt from interior gaps only (skip first and last); these
        # avoid the buggy first-frame Time and any tail anomaly.
        if n >= 4:
            interior = np.diff(times[1:-1])
            interior = interior[np.isfinite(interior) & (interior > 0)]
            median_dt = float(np.median(interior)) if interior.size else 60.0
        elif n >= 2:
            all_dts = np.diff(times)
            all_dts = all_dts[np.isfinite(all_dts) & (all_dts > 0)]
            median_dt = float(np.median(all_dts)) if all_dts.size else 60.0
        else:
            median_dt = 60.0
        # Default: rank within Path. Bump for any gap > 1.5x median dt.
        idx = np.arange(1, n + 1)
        for i in range(1, n):
            gap = times[i] - times[i - 1]
            if np.isfinite(gap) and gap > 1.5 * median_dt:
                n_missing = int(round(gap / median_dt)) - 1
                if n_missing > 0:
                    idx[i:] += n_missing
        df.loc[g.index, T_INDEX_COL] = idx
    return df


def find_ptu_for_path(ptu_folder: Path, path_value: str,
                      t_index: int | None = None) -> Path | None:
    """Resolve a spreadsheet 'Path' (+ optional t-index) to a PTU file.

    Tries the literal Path first. As a fallback, strips a trailing
    construct suffix like '-D48' / '-D90' / '-{letter}{digits}' from
    the position name -- useful when the live spreadsheet's Path includes
    the construct ("Position 1-D48") but the t0 PTU folder was named
    without it ("Position 1.ptu"). This lets one patched Path column
    cover both per-frame label lookup AND per-PTU lookup.
    """
    import re as _re
    stem = Path(path_value).name
    if stem.lower().endswith(".ptu"):
        stem = stem[:-4]
    if t_index is not None and int(t_index) > 0:
        stem = f"{stem}_t{int(t_index)}"
    cand = ptu_folder / f"{stem}.ptu"
    if cand.exists():
        return cand
    # Fallback: strip a trailing construct tag like "-D48" before the
    # optional "_t<N>" suffix.
    bare = _re.sub(r"-[A-Za-z]\d+(_t\d+)?$", r"\1", stem)
    if bare != stem:
        cand = ptu_folder / f"{bare}.ptu"
        if cand.exists():
            return cand
    for f in ptu_folder.glob("*.ptu"):
        if f.stem.lower() == stem.lower():
            return f
    return None


# --------------------------------------------------------------------------- #
# STAGE 2 -- per-pixel fitting
# --------------------------------------------------------------------------- #
def build_basis(irf_info: dict, taus, irf_shift_ns: float) -> np.ndarray:
    """Build the (n_time x 4) design matrix for one image.

    Columns: [conv(IRF_shifted, exp(-t/tau_i)) for i=1..3], plus a constant
    background column of ones. The shift applied is the global t0 plus this
    image's per-image IRF shift.
    """
    irf = irf_info["irf"]
    t = irf_info["t_ns"]
    dt = irf_info["dt_ns"]
    n = irf_info["n"]
    total_shift_bins = (irf_info["t0_shift_ns"] + irf_shift_ns) / dt
    irf_s = circular_shift(irf, total_shift_bins)
    cols = []
    for tau in taus:
        decay = np.exp(-t / float(tau))      # normalised to 1 at t=0
        cols.append(cconv(irf_s, decay))
    cols.append(np.ones(n))
    return np.column_stack(cols)


def _weighted_means(alphas, taus):
    """Return (tau_amp, tau_int) from amplitudes (axis 0 = component)."""
    sum_a = alphas.sum(axis=0)
    sum_at = (alphas * taus[..., None] if alphas.ndim > 1 else alphas * taus).sum(axis=0)
    sum_at2 = (alphas * (taus ** 2)[..., None] if alphas.ndim > 1
               else alphas * taus ** 2).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        tau_amp = np.where(sum_a > 0, sum_at / sum_a, 0.0)
        tau_int = np.where(sum_at > 0, sum_at2 / sum_at, 0.0)
    return tau_amp, tau_int


def build_keep_mask(t_ns: np.ndarray, exclude_ranges) -> np.ndarray:
    """Boolean mask over time bins, False inside any (lo, hi) ns range to drop
    those bins from the fit (e.g. to ignore a known detector artifact)."""
    keep = np.ones(t_ns.size, dtype=bool)
    for lo, hi in (exclude_ranges or []):
        keep &= ~((t_ns >= lo) & (t_ns <= hi))
    return keep


def fit_image_per_pixel(ptu_path: Path, irf_info: dict, row: pd.Series,
                        solver: str = "nnls", min_photons: float = 20.0,
                        exclude_ranges=None) -> dict:
    """Per-pixel amplitude fit for one experimental image.

    Only pixels whose total photon count exceeds `min_photons` are fitted; the
    rest are left at zero (this both speeds up the solve and avoids fitting pure
    noise, which otherwise biases the amplitude split toward the fast component).

    For the validation summary we ALSO fit the whole-image summed decay (NNLS),
    which is the direct analog of the LAS X global fit and is the robust value to
    compare against the spreadsheet.

    Returns a dict of 2D maps plus aggregate values.
    """
    cube, n, dt_ns, _dark = read_ptu_cube(ptu_path)
    if n != irf_info["n"]:
        log.warning("  %s: period %d != IRF period %d; cropping to min",
                    ptu_path.name, n, irf_info["n"])
        n = min(n, irf_info["n"])
        cube = cube[..., :n]
    H, W = cube.shape[:2]
    taus = np.array([float(row[c]) for c in COL_TAU])
    irf_shift = float(row.get(COL_IRF_SHIFT, 0.0) or 0.0)
    irf_bg = float(row.get(COL_IRF_BG, 0.0) or 0.0)

    A = build_basis(irf_info, taus, irf_shift)[:n]           # (n, 4)
    Y = cube.reshape(-1, n).astype(np.float64)               # (n_pix, n)
    if irf_bg:
        Y = Y - irf_bg                                       # fixed IRF offset
    intensity = Y.clip(min=0).sum(axis=1)                    # total photons/pixel
    npix = Y.shape[0]
    mask = intensity > min_photons                           # bright-pixel mask

    # Optionally drop time bins (e.g. a known artifact range) from the fit.
    keep = build_keep_mask(irf_info["t_ns"][:n], exclude_ranges)
    Af = A[keep]
    dof = max(int(keep.sum()) - 4, 1)

    # ---- global summed-decay fit (LAS X analog, robust validation value) ----
    # Sum only over bright pixels: dim/background pixels each carry a non-
    # fluorescence baseline (dark current, scatter) that the 3-exp model
    # cannot capture, so including them biases NNLS amplitudes into tau_1.
    from scipy.optimize import nnls
    Dsum = Y[mask].sum(axis=0) if mask.any() else Y.sum(axis=0)
    xg, _ = nnls(Af, Dsum[keep])
    ag = xg[:3]
    gmodel = Af @ xg
    g_chi2 = float(np.sum((Dsum[keep] - gmodel) ** 2 / np.clip(gmodel, 1.0, None)) / dof)
    g_amp = float((ag * taus).sum() / ag.sum()) if ag.sum() > 0 else float("nan")
    g_int = float((ag * taus ** 2).sum() / (ag * taus).sum()) if (ag * taus).sum() > 0 else float("nan")

    # ---- per-pixel fit over bright pixels only ----
    Ym = Y[mask][:, keep].T                                  # (n_keep, n_bright)
    X = np.zeros((4, npix))
    if Ym.shape[1]:
        if solver == "nnls":
            Xm = np.zeros((4, Ym.shape[1]))
            for j in range(Ym.shape[1]):
                Xm[:, j], _ = nnls(Af, Ym[:, j])
        else:  # lstsq + clip
            Xm, *_ = np.linalg.lstsq(Af, Ym, rcond=None)
            Xm[:3] = np.clip(Xm[:3], 0.0, None)
        X[:, mask] = Xm

    alphas = X[:3]                                           # (3, n_pix)
    tau_amp, tau_int = _weighted_means(alphas, taus)

    # Per-pixel reduced chi2 (Poisson-weighted) over kept bins, where fitted.
    model = Af @ X
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2_pix = np.sum((Y[:, keep].T - model) ** 2 / np.clip(model, 1.0, None), axis=0) / dof
    chi2_pix[~mask] = 0.0

    def m(v):
        return v.reshape(H, W).astype(np.float32)

    maps = dict(
        alpha1=m(alphas[0]), alpha2=m(alphas[1]), alpha3=m(alphas[2]),
        tau_int=m(tau_int), tau_amp=m(tau_amp),
        intensity=m(intensity), chi2=m(chi2_pix),
    )

    # Map-based aggregate (intensity-weighted over fitted pixels) -- reported for
    # transparency; the global fit value above is the primary validation number.
    Atot = alphas[:, mask].sum(axis=1)
    map_int = float((Atot * taus ** 2).sum() / (Atot * taus).sum()) if (Atot * taus).sum() > 0 else float("nan")
    map_amp = float((Atot * taus).sum() / Atot.sum()) if Atot.sum() > 0 else float("nan")

    agg = dict(
        n_photons_total=float(intensity.sum()),
        n_pixels_fit=int(mask.sum()),
        mean_tau_int=g_int,            # from global summed-decay fit (LAS X analog)
        mean_tau_amp=g_amp,
        mean_tau_int_map=map_int,      # from per-pixel maps (intensity weighted)
        mean_tau_amp_map=map_amp,
        global_chi2=g_chi2,
        mean_chi2=float(np.nanmedian(chi2_pix[mask])) if mask.any() else float("nan"),
    )
    return dict(maps=maps, agg=agg, taus=taus)


def save_outputs(image_name: str, result: dict, output_dir: Path,
                 flim_min: float, flim_max: float):
    """Write the per-image TIFFs and QC panel.

    Both raw maps are written as uint16:
      intensity: raw photon counts, clipped to [0, 65535].
      tau_int:   tau * 1000, i.e. ps, rounded and clipped to [0, 65535].
                 0 means no valid lifetime (matches the legacy LAS X export).
    """
    maps = result["maps"]
    img_dir = output_dir
    img_dir.mkdir(parents=True, exist_ok=True)

    intensity_u16 = np.clip(np.rint(maps["intensity"]), 0, 65535).astype(np.uint16)
    tau_int_u16 = np.clip(np.rint(maps["tau_int"] * 1000.0), 0, 65535).astype(np.uint16)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import tifffile
    # Colormapping helper; a sibling of this file, so make sure it is importable
    # even when flim_pipeline was imported from elsewhere.
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    import FLIM_vis_v15 as fvis
    tifffile.imwrite(str(img_dir / f"{image_name}_tau_int.tif"), tau_int_u16)
    tifffile.imwrite(str(img_dir / f"{image_name}_intensity.tif"), intensity_u16)

    # Colormapped TIFF via FLIM_vis_v15 (hue=lifetime, brightness=photons).
    intensity = maps["intensity"]
    tau_int = maps["tau_int"]
    pos = intensity[intensity > 0]
    pc_min = float(np.percentile(pos, 1)) if pos.size else 0.0
    pc_max = float(np.percentile(pos, 99)) if pos.size else 1.0
    rgb = fvis.apply_flim_coloring(
        photon=intensity, flim_ns=tau_int,
        pc_min=pc_min, pc_max=pc_max,
        flim_ns_min=flim_min, flim_ns_max=flim_max,
        invalid_mask=(tau_int <= 0),
    )
    fvis.write_rgb_tiff(img_dir / f"{image_name}_tau_int_colormap.tif", rgb)

    # 4-panel QC PNG.
    fig, ax = plt.subplots(2, 2, figsize=(10, 9))
    im0 = ax[0, 0].imshow(intensity, cmap="gray")
    ax[0, 0].set_title("intensity (photons)")
    fig.colorbar(im0, ax=ax[0, 0], fraction=0.046)
    im1 = ax[0, 1].imshow(tau_int, cmap=fvis.CUSTOM_CMAP, vmin=flim_min, vmax=flim_max)
    ax[0, 1].set_title("τ_int (ns)")
    fig.colorbar(im1, ax=ax[0, 1], fraction=0.046)
    im2 = ax[1, 0].imshow(maps["tau_amp"], cmap=fvis.CUSTOM_CMAP,
                          vmin=flim_min, vmax=flim_max)
    ax[1, 0].set_title("τ_amp (ns)")
    fig.colorbar(im2, ax=ax[1, 0], fraction=0.046)
    chi2 = maps["chi2"]
    vmax = float(np.percentile(chi2[np.isfinite(chi2)], 99)) if np.any(np.isfinite(chi2)) else 1.0
    im3 = ax[1, 1].imshow(chi2, cmap="magma", vmin=0, vmax=max(vmax, 1e-3))
    ax[1, 1].set_title("reduced χ²")
    fig.colorbar(im3, ax=ax[1, 1], fraction=0.046)
    for a in ax.ravel():
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(image_name)
    fig.tight_layout()
    fig.savefig(img_dir / f"{image_name}_qc.png", dpi=120)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# STAGE 3 -- validation
# --------------------------------------------------------------------------- #
def validate_against_lasx(summary_rows, output_dir: Path,
                          abs_tol: float = 0.1, rel_tol: float = 0.05):
    """Write summary.csv and validate against LAS X.

    Intensity-weighted mean lifetime (tau_int) is the primary validation target
    and is checked per-image with the requested tolerances. Amplitude-weighted
    mean (tau_amp) is reported as a cross-check only: with the fast component
    (tau_1 ~ 0.2 ns) nearly degenerate with the IRF, the amplitude split is
    ill-posed for fixed-lifetime unmixing and is NOT expected to reproduce LAS X,
    so it is summarised in aggregate rather than warned per-image.
    """
    df = pd.DataFrame(summary_rows)
    csv_path = output_dir / "summary.csv"
    df.to_csv(csv_path, index=False)
    log.info("STAGE 3: wrote %s (%d images)", csv_path, len(df))

    # Primary: tau_int agreement.
    n_warn = 0
    for _, r in df.iterrows():
        py, lasx = r["mean_tau_int"], r["mean_tau_int_LASX"]
        if not (np.isfinite(py) and np.isfinite(lasx)) or lasx == 0:
            continue
        diff = abs(py - lasx)
        rel = diff / abs(lasx)
        if diff > abs_tol and rel > rel_tol:
            n_warn += 1
            log.warning("  τ_int mismatch %s: Python=%.3f LAS X=%.3f (Δ=%+.3f ns, %.1f%%) "
                        "[global χ²=%.1f]", r["image_name"], py, lasx,
                        py - lasx, rel * 100, r.get("global_chi2", float("nan")))

    # Summary stats per group.
    for grp in sorted({str(n).split("-")[0] for n in df["image_name"]}):
        sub = df[df["image_name"].astype(str).str.startswith(grp + "-")]
        if sub.empty:
            continue
        di = (sub["mean_tau_int"] - sub["mean_tau_int_LASX"])
        da = (sub["mean_tau_amp"] - sub["mean_tau_amp_LASX"])
        log.info("  [%s] τ_int Δ = %+.3f±%.3f ns (median |Δ| %.3f); "
                 "τ_amp Δ = %+.3f ns (cross-check only); median global χ² = %.2f",
                 grp, di.mean(), di.std(), di.abs().median(),
                 da.mean(), sub["global_chi2"].median())
    log.info("  τ_int: %d/%d images exceed tolerance (%.2f ns & %.0f%%).",
             n_warn, len(df), abs_tol, rel_tol * 100)
    log.info("  Note: τ_amp is an ill-posed cross-check (fast component ≈ IRF) and "
             "is expected to read low vs LAS X; τ_int is the validated quantity.")
    return df


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def setup_logging(output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%H:%M:%S")
    sh = logging.StreamHandler(sys.stdout); sh.setFormatter(fmt); log.addHandler(sh)
    fh = logging.FileHandler(output_dir / "flim_pipeline.log"); fh.setFormatter(fmt)
    log.addHandler(fh)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--spreadsheets", nargs="+", type=Path, required=True,
                   help="LAS X global-fit spreadsheet(s) (.xlsx) for this dataset")
    p.add_argument("--ptu-folder", type=Path, required=True,
                   help="folder of the experimental PTU files (a LAS X .sptw)")
    p.add_argument("--irf-ptu", type=Path, required=True,
                   help="Rhodamine B reference PTU the IRF is extracted from")
    p.add_argument("--tau-r", type=float, default=DEFAULT_TAU_R,
                   help="reference lifetime of the Rhodamine B standard (ns); "
                        "must match --irf-ptu (default %(default)s)")
    p.add_argument("--ref-components", nargs="+", default=None, metavar="AMP,TAU",
                   help="the reference's OWN decay model, as AMP,TAU pairs from the LAS X fit of the Rhodamine B reference (e.g. --ref-components 0.0918,0.2868 0.9082,2.1146). Overrides --tau-r for the deconvolution kernel. Use this whenever the reference needs >1 component to fit well: stored dye aliquots in this study need two, and deconvolving them with a mono-exponential kernel distorts the IRF. Omit for a genuinely mono-exponential (freshly prepared) reference.")
    p.add_argument("--output-dir", type=Path, required=True,
                   help="where to write the maps, IRF vector and summary.csv")
    p.add_argument("--rl-iters", type=int, default=400,
                   help="Max regularized Richardson-Lucy iterations (a cap; the "
                        "extraction normally stops earlier on convergence).")
    p.add_argument("--irf-sigma-bins", type=float, default=0.4,
                   help="Gaussian damping scale (in time bins) for the regularized "
                        "RL IRF extraction. 0.4 -> FWHM ~388 ps, the chi2-optimal "
                        "width for the 20250715 sessions (and the screenshot match). "
                        "0.6 -> ~485 ps centres mean_tau_int on LAS X but fits ~3x "
                        "worse. Set 0 for plain (un-regularized) RL, which "
                        "over-sharpens high-SNR references and biases tau_int high.")
    p.add_argument("--irf-rl-tol", type=float, default=1e-3,
                   help="Convergence tolerance on the relative L1 step of the IRF "
                        "estimate; RL stops once below this (default 1e-3).")
    p.add_argument("--irf-extra-shift-ns", type=float, default=0.0,
                   help="Extra IRF time shift (ns) added to the reference-fitted "
                        "t0, to align the IRF to this acquisition session's laser "
                        "time-zero (analog of LAS X 'IRF Shift'). Empirically "
                        "-0.04 for 20250715 sessions, -0.14 for 20250630.")
    p.add_argument("--flim-min", type=float, default=1.5, help="colormap τ min (ns)")
    p.add_argument("--flim-max", type=float, default=3.0, help="colormap τ max (ns)")
    p.add_argument("--solver", choices=["nnls", "lstsq"], default="nnls",
                   help="Per-pixel solver: 'nnls' (non-negative, default) or "
                        "'lstsq' (batched least-squares with negatives clipped).")
    p.add_argument("--min-photons", type=float, default=50.0,
                   help="Only fit pixels with total photons above this (default 50); "
                        "dimmer pixels are left at zero. The same filter is "
                        "applied to the global summed-decay validation fit, "
                        "so dark/background pixels don't drag the model.")
    p.add_argument("--exclude-ns", nargs="+", default=None, metavar="LO-HI",
                   help="Time range(s) in ns to drop from the fit, e.g. '7.0-7.6'. "
                        "Useful to ignore a known detector artifact. Repeatable.")
    p.add_argument("--stage", choices=["1", "all"], default="all",
                   help="'1' = IRF extraction only then stop; 'all' = full pipeline.")
    p.add_argument("--limit", type=int, default=None,
                   help="Process only the first N images (debugging).")
    p.add_argument("--skip-tindex", nargs="+", type=int, default=None,
                   metavar="N",
                   help="Drop these t-index values from a time-lapse run "
                        "(e.g. '--skip-tindex 1' for live data where the "
                        "first frame has focus drift after drug addition). "
                        "Only relevant when the spreadsheet has a Time column.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    setup_logging(args.output_dir)
    log.info("FLIM pipeline starting; output -> %s", args.output_dir)

    # Parse optional artifact-exclusion time ranges ("lo-hi" in ns).
    exclude_ranges = []
    for tok in (args.exclude_ns or []):
        lo, hi = (float(v) for v in tok.split("-"))
        exclude_ranges.append((lo, hi))
    if exclude_ranges:
        log.info("Excluding time ranges from fit (ns): %s", exclude_ranges)

    # STAGE 1
    irf_info = extract_irf(args.irf_ptu, args.tau_r, args.output_dir,
                           rl_iters=args.rl_iters, sigma_bins=args.irf_sigma_bins,
                           rl_tol=args.irf_rl_tol, extra_shift_ns=args.irf_extra_shift_ns,
                           ref_components=parse_ref_components(args.ref_components))

    if args.stage == "1":
        log.info("Stage 1 complete. Inspect outputs/irf_diagnostics.png, then "
                 "re-run with --stage all to process the experimental images.")
        return 0

    # STAGE 2
    df = load_spreadsheets(args.spreadsheets)
    time_lapse = T_INDEX_COL in df.columns
    if time_lapse and args.skip_tindex:
        skip = set(int(x) for x in args.skip_tindex)
        before = len(df)
        df = df[~df[T_INDEX_COL].isin(skip)].copy()
        log.info("  --skip-tindex %s: %d -> %d rows", sorted(skip), before, len(df))
    if args.limit:
        df = df.iloc[:args.limit]
    log.info("STAGE 2: per-pixel fitting for %d images", len(df))

    summary_rows = []
    pixel_rows = []
    from tqdm import tqdm
    for _, row in tqdm(df.iterrows(), total=len(df), desc="images"):
        name = str(row[COL_PATH])
        t_index = int(row[T_INDEX_COL]) if time_lapse else None
        ptu = find_ptu_for_path(args.ptu_folder, name, t_index=t_index)
        if ptu is None:
            label = f"{name}_t{t_index}" if t_index else name
            log.warning("  no PTU found for '%s'; skipping", label)
            continue
        pixel_rows.append(dict(image_name=ptu.stem, pixel_um=read_pixel_um(ptu)))
        try:
            res = fit_image_per_pixel(ptu, irf_info, row,
                                      solver=args.solver, min_photons=args.min_photons,
                                      exclude_ranges=exclude_ranges)
        except Exception as e:  # noqa: BLE001
            log.exception("  failed on %s: %s", ptu.stem, e)
            continue
        image_name = ptu.stem      # e.g. "Position 1" or "Position 1_t42"
        save_outputs(image_name, res, args.output_dir, args.flim_min, args.flim_max)
        agg = res["agg"]
        lasx_int = float(row.get(COL_MEAN_INT_LASX, np.nan))
        lasx_amp = float(row.get(COL_MEAN_AMP_LASX, np.nan))
        summary_rows.append(dict(
            image_name=image_name,
            path=name,
            t_index=t_index if t_index else "",
            n_photons_total=agg["n_photons_total"],
            n_pixels_fit=agg["n_pixels_fit"],
            mean_tau_int=agg["mean_tau_int"],
            mean_tau_int_LASX=lasx_int,
            diff_tau_int=agg["mean_tau_int"] - lasx_int,
            mean_tau_amp=agg["mean_tau_amp"],
            mean_tau_amp_LASX=lasx_amp,
            diff_tau_amp=agg["mean_tau_amp"] - lasx_amp,
            mean_tau_int_map=agg["mean_tau_int_map"],
            mean_tau_amp_map=agg["mean_tau_amp_map"],
            global_chi2=agg["global_chi2"],
            mean_chi2=agg["mean_chi2"],
        ))

    # Per-image pixel-size sidecar (for correct ROI area across mixed zoom).
    if pixel_rows:
        px_path = args.output_dir / "pixel_um.csv"
        pd.DataFrame(pixel_rows).to_csv(px_path, index=False)
        log.info("Wrote per-image pixel sizes -> %s", px_path)

    # STAGE 3
    if summary_rows:
        validate_against_lasx(summary_rows, args.output_dir)
    log.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
