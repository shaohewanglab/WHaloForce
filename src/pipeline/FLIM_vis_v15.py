#!/usr/bin/env python3
"""
FLIM_vis.py (v15)

Create RGB TIFFs where hue encodes FLIM lifetime (ns) using a custom
Blue→Cyan→Green→Yellow→Red colormap, and intensity encodes photon counts.

Supported modes:
1) Export-root mode:
   - Point the script at the *root export folder* produced by the LIF exporter
     (e.g., containing "01_FastFLIM_IRF_corrected/", "02_Phasor_phase_lifetime/", ...,
     and "06_Intensity_raw/").
   - The script will:
       * find all "lifetime" folders (and the FastFLIM folder),
       * read each lifetime TIFF,
       * pair it with the matching raw intensity TIFF from the intensity folder
         (same relative path under "scene_.../"),
       * write an RGB TIFF into a *parallel* output folder next to the lifetime folder,
         e.g. "02_Phasor_phase_lifetime_FLIM_RGB/...",
       * write one legend PNG per output folder.

2) Folder-pairing mode (legacy):
   - Input is a folder containing TIFF files named with `ch0/ch1/ch2...`.
   - Pairs are inferred as (ch0=photon, ch1=flim), (ch2=photon, ch3=flim), ...

3) Two-file mode (legacy):
   - Input is two files: <photon> <flim>.

Notes:
- Lifetime TIFFs are assumed to be 16-bit grayscale where
  gray_value = lifetime_ns * flim_unit_scale (default flim_unit_scale=1000.0).
- Pixels with lifetime gray value == 0 are treated as invalid lifetime
  and forced to black (even if photon counts are nonzero).

v9 changes:
- In export-root mode, output folder names include the FLIM lifetime
  visualization range, e.g. "..._FLIM_RGB_1.0to2.5" for "-l 1.0 -L 2.5".
  This allows keeping multiple color scales side-by-side.

v10 changes:
- In export-root mode, also process the raw intensity folder (default
  "06_Intensity_raw") and export per-slice 8-bit grayscale TIFFs in a
  parallel folder. Scaling is percentile-based with a single saturation
  parameter: --intensity-sat S uses low/high percentiles (S, 100-S).
  Default S=1 → 1–99%.

v11 changes:
- Fix argparse help-string formatting error by escaping percent signs ("%%").


v12 changes:
- Add median-centered lifetime scaling via --flim-span. When provided, each image is scaled to
  [median - span/2, median + span/2] using the median lifetime of that image (computed over valid pixels).
- RGB output folder tags now encode either fixed bounds (_<min>to<max>) or span-based scaling (_medspan<span>).

v13 changes:
- Write a CSV summary listing per-image FLIM scaling values (median, vmin, vmax, actual span, valid-pixel count).
- Ensure median lifetime is computed over valid pixels only (lifetime gray > 0), excluding background zeros.

v14 changes:
- In --flim-span mode, compute the median lifetime after excluding dim pixels using an intensity threshold.
  This helps prevent low-signal/noisy pixels from biasing the per-image median used for the color scale.
- Add CLI options to control the intensity gate used for median computation, and report the gate statistics
  (threshold + pixel counts) in the summary CSV.



v15 changes:
- In export-root mode, do not overwrite 8-bit intensity outputs by default; skip conversion for files that already exist (use --overwrite-intensity-8bit to force).
- Add --only-prefix option to process only selected lifetime folders in export-root mode based on the leading numeric folder prefix (e.g. 02, 03, 04).
"""

import argparse
import csv
import re
from pathlib import Path
import sys

import numpy as np
import tifffile
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


# ------------------------ Colormap -----------------------------------------

def create_custom_colormap() -> LinearSegmentedColormap:
    """Blue → Cyan → Green → Yellow → Red continuous colormap."""
    cdict = {
        "red": [
            (0.0, 0.0, 0.0),
            (0.25, 0.0, 0.0),
            (0.5, 0.0, 0.0),
            (0.75, 1.0, 1.0),
            (1.0, 1.0, 1.0),
        ],
        "green": [
            (0.0, 0.0, 0.0),
            (0.25, 1.0, 1.0),
            (0.5, 1.0, 1.0),
            (0.75, 1.0, 1.0),
            (1.0, 0.0, 0.0),
        ],
        "blue": [
            (0.0, 1.0, 1.0),
            (0.25, 1.0, 1.0),
            (0.5, 0.0, 0.0),
            (0.75, 0.0, 0.0),
            (1.0, 0.0, 0.0),
        ],
    }
    return LinearSegmentedColormap("custom_colormap", cdict)

CUSTOM_CMAP = create_custom_colormap()

# ------------------------ Utilities ----------------------------------------

def _fmt_float_for_tag(x: float) -> str:
    """Format a float for folder-name tags.

    Rules:
    - Integers are shown with one decimal (1.0, 3.0) to match common FLIM range notation.
    - Non-integers are shown up to 6 decimals with trailing zeros removed.
    """
    if x is None:
        return "none"
    try:
        xf = float(x)
    except Exception:
        return str(x)

    if np.isfinite(xf) and abs(xf - round(xf)) < 1e-9:
        return f"{xf:.1f}"
    s = f"{xf:.6f}".rstrip("0").rstrip(".")
    return s if s else "0.0"


def flim_range_tag(flim_ns_min: float, flim_ns_max: float) -> str:
    """Return a folder-name tag like '_1.0to2.5' from FLIM bounds."""
    return f"_{_fmt_float_for_tag(flim_ns_min)}to{_fmt_float_for_tag(flim_ns_max)}"



def flim_span_tag(span: float) -> str:
    """Return a folder-name tag like '_medspan1.5' from a span in ns."""
    return f"_medspan{_fmt_float_for_tag(span)}"


def flim_vis_tag(args: argparse.Namespace) -> str:
    """Return the folder-name tag used to distinguish visualization scales.

    - If --flim-span is provided: use '_medspan<span>'
    - Otherwise: use '_<min>to<max>'
    """
    if getattr(args, "flim_span", None) is not None:
        return flim_span_tag(float(args.flim_span))
    return flim_range_tag(float(args.flim_ns_min), float(args.flim_ns_max))


def intensity_sat_tag(sat: float) -> str:
    """Return a folder-name tag like '_p1to99' from saturation percentile.

    We interpret `sat` as symmetric clipping in percent:
      low = sat
      high = 100 - sat
    """
    satf = float(sat)
    lo = _fmt_float_for_tag(satf)
    hi = _fmt_float_for_tag(100.0 - satf)
    return f"_p{lo}to{hi}"


def ensure_2d(arr: np.ndarray) -> np.ndarray:
    arr = np.squeeze(arr)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D image after squeeze, got {arr.shape}")
    return arr


def normalize(v: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    v = np.clip(v, vmin, vmax)
    denom = (vmax - vmin) if vmax > vmin else 1.0
    return (v - vmin) / denom


def pct(v: np.ndarray, p: float) -> float:
    return float(np.percentile(v, p))



def compute_used_ranges(
    photon: np.ndarray,
    flim_gray: np.ndarray,
    flim_unit_scale: float,
    pc_min: float | None,
    pc_max: float | None,
    pc_pmin: float | None,
    pc_pmax: float | None,
    flim_ns_min: float | None,
    flim_ns_max: float | None,
    flim_span: float | None = None,
    median_intensity_gate: bool = True,
    median_intensity_min: float | None = None,
    median_intensity_pmin: float | None = 10.0,
) -> tuple[float, float, float, float, np.ndarray, float | None, dict[str, object]]:
    """Return (pc_min, pc_max, flim_ns_min, flim_ns_max, flim_ns_array, flim_median, flim_median_meta).

    Photon scaling:
      - If pc_min/pc_max are not provided, robust percentiles are used (pc_pmin/pc_pmax).

    Lifetime scaling:
      - If flim_span is provided, each image is scaled to:
            [median - flim_span/2, median + flim_span/2]
        where median is computed over *valid lifetime pixels* (flim_gray > 0).
        In addition, if `median_intensity_gate` is True, the median is computed
        after excluding dim pixels using either an absolute intensity threshold
        (`median_intensity_min`) or a percentile threshold (`median_intensity_pmin`)
        computed over the intensities of valid-lifetime pixels.
      - Otherwise, fixed bounds are used (defaults 1.5–3.0 ns unless overridden).
    """
    photon = ensure_2d(np.asarray(photon))
    flim_gray = ensure_2d(np.asarray(flim_gray))
    flim_ns = flim_gray.astype(np.float32) / float(flim_unit_scale)

    # photon range (robust by default)
    if pc_min is None or pc_max is None:
        pc_vals = photon[photon > 0] if np.any(photon > 0) else photon
        if pc_min is None:
            pc_min = pct(pc_vals, pc_pmin) if pc_pmin is not None else float(np.min(pc_vals))
        if pc_max is None:
            pc_max = pct(pc_vals, pc_pmax) if pc_pmax is not None else float(np.max(pc_vals))

    flim_median: float | None = None
    flim_median_meta: dict[str, object] = {}

    # FLIM range: span-based (median-centered) if requested
    if flim_span is not None:
        span = float(flim_span)
        if span <= 0:
            raise ValueError("--flim-span must be > 0")
        valid = np.isfinite(flim_ns) & (flim_gray > 0)

        candidate_pixels = int(np.count_nonzero(valid))
        used_mask = valid
        used_pixels = candidate_pixels
        intensity_threshold: float | None = None
        fallback_used = False

        # Intensity-gated median (helps avoid dim/noisy pixels biasing the per-image center).
        if median_intensity_gate and candidate_pixels > 0:
            # Determine threshold.
            if median_intensity_min is not None:
                intensity_threshold = float(median_intensity_min)
            else:
                p = None if median_intensity_pmin is None else float(median_intensity_pmin)
                if p is not None:
                    # Clamp percentile to [0, 100]
                    p = float(np.clip(p, 0.0, 100.0))
                    vals = np.asarray(photon, dtype=np.float32)[valid]
                    vals = vals[np.isfinite(vals)]
                    # Prefer strictly positive values if available (avoids percentile=0 when zeros present).
                    vals_pos = vals[vals > 0]
                    if vals_pos.size > 0:
                        vals = vals_pos
                    if vals.size > 0:
                        intensity_threshold = float(np.nanpercentile(vals, p))

            # Apply the gate to compute the median.
            if intensity_threshold is not None and np.isfinite(intensity_threshold):
                used_mask = valid & (np.asarray(photon, dtype=np.float32) >= intensity_threshold)
                used_pixels = int(np.count_nonzero(used_mask))
                if used_pixels == 0:
                    # If gating eliminates all pixels, fall back to the ungated valid mask.
                    fallback_used = True
                    used_mask = valid
                    used_pixels = candidate_pixels
                    intensity_threshold = None

        if used_pixels > 0:
            flim_median = float(np.median(flim_ns[used_mask]))
            flim_ns_min_used = flim_median - span / 2.0
            flim_ns_max_used = flim_median + span / 2.0
            # Lifetimes cannot be negative; clip lower bound to 0.
            if flim_ns_min_used < 0.0:
                flim_ns_min_used = 0.0
            # Guard against degenerate range after clipping.
            if flim_ns_max_used <= flim_ns_min_used:
                flim_ns_max_used = flim_ns_min_used + 1e-6

            flim_median_meta = {
                "median_intensity_gate": bool(median_intensity_gate),
                "median_intensity_min": float(median_intensity_min) if median_intensity_min is not None else "",
                "median_intensity_pmin": float(median_intensity_pmin) if (median_intensity_min is None and median_intensity_pmin is not None) else "",
                "median_intensity_threshold": float(intensity_threshold) if intensity_threshold is not None else "",
                "median_candidate_pixels": candidate_pixels,
                "median_pixels_used": used_pixels,
                "median_used_fraction": (used_pixels / candidate_pixels) if candidate_pixels else "",
                "median_gate_fallback": bool(fallback_used),
            }

            return float(pc_min), float(pc_max), float(flim_ns_min_used), float(flim_ns_max_used), flim_ns, flim_median, flim_median_meta

    # FLIM range: fixed defaults unless overridden
    if flim_ns_min is None:
        flim_ns_min = 1.5
    if flim_ns_max is None:
        flim_ns_max = 3.0

    return float(pc_min), float(pc_max), float(flim_ns_min), float(flim_ns_max), flim_ns, flim_median, flim_median_meta


def compute_intensity_8bit(
    intensity: np.ndarray,
    sat: float = 1.0,
) -> tuple[np.ndarray, float, float]:
    """Scale an intensity image to uint8 using symmetric percentile saturation.

    Parameters
    ----------
    intensity:
        2D (or squeezeable) intensity image.
    sat:
        Saturation level in percent. We map percentiles (sat, 100-sat) to (0,255).
        Default sat=1.0 corresponds to 1–99% scaling.

    Returns
    -------
    img_u8:
        uint8 image scaled to [0,255].
    vmin, vmax:
        The intensity values mapped to 0 and 255, respectively (per-image).
    """
    img = ensure_2d(np.asarray(intensity))
    im = np.asarray(img, dtype=np.float32)

    finite = np.isfinite(im)
    vals = im[finite & (im > 0)]
    if vals.size == 0:
        vals = im[finite]

    if vals.size == 0:
        return np.zeros_like(img, dtype=np.uint8), 0.0, 0.0

    satf = float(sat)
    p_lo = satf
    p_hi = 100.0 - satf

    # Robust percentile scaling. Use nanpercentile for safety.
    vmin = float(np.nanpercentile(vals, p_lo))
    vmax = float(np.nanpercentile(vals, p_hi))

    if not np.isfinite(vmin):
        vmin = 0.0
    if not np.isfinite(vmax):
        vmax = vmin

    denom = vmax - vmin
    if denom <= 0:
        scaled = np.zeros_like(im, dtype=np.float32)
    else:
        scaled = (im - vmin) / denom

    scaled = np.clip(scaled, 0.0, 1.0)
    scaled[~finite] = 0.0
    img_u8 = np.round(scaled * 255.0).astype(np.uint8)
    return img_u8, vmin, vmax


# ------------------------ Core logic ---------------------------------------

def apply_flim_coloring(
    photon: np.ndarray,
    flim_ns: np.ndarray,
    pc_min: float,
    pc_max: float,
    flim_ns_min: float,
    flim_ns_max: float,
    invalid_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Return uint8 (Y,X,3) RGB where hue encodes lifetime and intensity encodes photons."""
    photon = ensure_2d(np.asarray(photon))
    flim_ns = ensure_2d(np.asarray(flim_ns))
    if photon.shape != flim_ns.shape:
        raise ValueError(f"Photon and FLIM arrays must have same shape, got {photon.shape} vs {flim_ns.shape}")

    pc_norm = normalize(photon.astype(np.float32), pc_min, pc_max)
    flim_norm = normalize(flim_ns.astype(np.float32), flim_ns_min, flim_ns_max)

    if invalid_mask is not None:
        invalid_mask = ensure_2d(np.asarray(invalid_mask).astype(bool))
        if invalid_mask.shape != photon.shape:
            raise ValueError(f"invalid_mask shape {invalid_mask.shape} does not match image shape {photon.shape}")
        # Force invalid-lifetime pixels to black regardless of photon counts
        pc_norm = pc_norm.copy()
        pc_norm[invalid_mask] = 0.0

    rgb = CUSTOM_CMAP(flim_norm)[..., :3]  # float in [0,1]
    out = (rgb * pc_norm[..., None])
    return np.clip(out * 255.0, 0, 255).astype(np.uint8)


# ------------------------ I/O ----------------------------------------------

def write_rgb_tiff(
    out_path: Path,
    rgb: np.ndarray,
    compress: str = "deflate",
    compression_level: int | None = 6,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        str(out_path),
        rgb,
        photometric="rgb",
        compression=compress,
        compressionargs={"level": compression_level} if compression_level is not None else None,
    )


def write_gray_u8_tiff(
    out_path: Path,
    img_u8: np.ndarray,
    compress: str = "deflate",
    compression_level: int | None = 6,
) -> None:
    """Write an 8-bit grayscale TIFF."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        str(out_path),
        np.asarray(img_u8, dtype=np.uint8),
        photometric="minisblack",
        compression=compress,
        compressionargs={"level": compression_level} if compression_level is not None else None,
    )


def save_legend_png(png_path: Path, vmin: float, vmax: float, title: str = "FLIM lifetime (ns)") -> None:
    png_path.parent.mkdir(parents=True, exist_ok=True)
    grad = np.linspace(0, 1, 512)[None, :]
    fig, ax = plt.subplots(figsize=(8, 1.0), dpi=150)
    ax.imshow(grad, aspect="auto", cmap=CUSTOM_CMAP, extent=[vmin, vmax, 0, 1])
    ax.set_yticks([])
    ax.set_xlabel(title)
    ax.set_xlim(vmin, vmax)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    fig.savefig(str(png_path), bbox_inches="tight")
    plt.close(fig)


# ------------------------ Legacy folder pairing logic ----------------------

CH_REGEX = re.compile(r"^(?P<prefix>.*?)(?:[_-]?ch)(?P<idx>\d+)(?P<suffix>[^\d].*)?$", re.IGNORECASE)


def find_pairs_in_folder(folder: Path) -> list[tuple[Path, Path, str, int]]:
    """Return list of (photon_path, flim_path, prefix, pair_index).

    Pair index is floor(ch/2). Files are paired as (ch0,ch1), (ch2,ch3), ...
    The output prefix removes the `ch<idx>` token from the filename stem.
    """
    files = sorted([p for p in folder.iterdir() if p.suffix.lower() in {".tif", ".tiff"}])
    buckets: dict[tuple[str, int], dict[str, Path]] = {}
    for p in files:
        m = CH_REGEX.match(p.stem)
        if not m:
            continue
        prefix = (m.group("prefix") or "") + (m.group("suffix") or "")
        ch = int(m.group("idx"))
        pair_idx = ch // 2
        key = (prefix, pair_idx)
        role = "photon" if (ch % 2 == 0) else "flim"
        buckets.setdefault(key, {})[role] = p

    pairs: list[tuple[Path, Path, str, int]] = []
    for (prefix, pair_idx), d in sorted(buckets.items()):
        if "photon" in d and "flim" in d:
            pairs.append((d["photon"], d["flim"], prefix, pair_idx))
    return pairs


# ------------------------ Export-root mode helpers -------------------------

def _has_tifs_recursively(folder: Path) -> bool:
    if not folder.is_dir():
        return False
    for p in folder.rglob("*"):
        if p.is_file() and p.suffix.lower() in (".tif", ".tiff"):
            return True
    return False


def looks_like_rgb_output_folder(folder: Path, rgb_suffix: str | None) -> bool:
    """Heuristic to detect this script's RGB output folders.

    We avoid skipping an arbitrary folder just because its name contains
    the suffix; instead we also check for the per-folder legend.png marker.
    """
    if not folder.is_dir():
        return False
    name = folder.name.lower()

    has_suffix_token = False
    if rgb_suffix:
        s = rgb_suffix.lower()
        has_suffix_token = name.endswith(s) or f"{s}_" in name

    if (has_suffix_token or ("flim_rgb" in name)) and (folder / "legend.png").is_file():
        return True
    return False


def is_lifetime_folder_name(
    name: str,
    intensity_folder: str,
    rgb_suffix: str | None = None,
) -> bool:
    """Return True if folder name looks like a lifetime/FastFLIM map folder.

    In export-root mode we want to *exclude* folders produced by this script
    (RGB output folders), so we also filter out folders that:
      * end with the configured rgb_suffix (default: _FLIM_RGB), or
      * contain the token 'flim_rgb' (covers the default suffix even if renamed).
    """
    lower = name.lower()

    if lower == intensity_folder.lower():
        return False

    # Skip this script's output folders (avoid re-processing on re-runs).
    # Output folders are named like: <lifetime_folder><rgb_suffix>_<min>to<max>
    if rgb_suffix:
        s = rgb_suffix.lower()
        if lower.endswith(s) or f"{s}_" in lower:
            return False
    if "flim_rgb" in lower:
        return False

    if "phasor_plot" in lower or "phasorplot" in lower:
        return False
    if "legend" in lower:
        return False

    # Lifetime maps: most contain "lifetime"; FastFLIM folder doesn't.
    return ("lifetime" in lower) or ("fastflim" in lower)


def detect_export_root(folder: Path, intensity_folder: str, rgb_suffix: str | None = None) -> bool:
    """Heuristic: export root contains the intensity folder and at least one lifetime folder."""
    if not folder.is_dir():
        return False
    intensity_path = folder / intensity_folder
    if not intensity_path.is_dir():
        return False
    # Must also have at least one candidate lifetime folder
    for child in folder.iterdir():
        if not child.is_dir():
            continue
        # Robustly skip this script's RGB output folders even if the suffix changed
        # between runs (legend.png is our marker file).
        if (child / "legend.png").is_file():
            continue
        if is_lifetime_folder_name(child.name, intensity_folder, rgb_suffix) and _has_tifs_recursively(child):
            return True
    return False


def find_lifetime_folders(
    root: Path,
    intensity_folder: str,
    rgb_suffix: str | None = None,
    only_prefixes: list[str] | None = None,
) -> list[Path]:
    """Return candidate lifetime-map folders under an export root.

    Parameters
    ----------
    only_prefixes:
        Optional list of numeric prefixes (e.g. ["02", "04"]) to restrict which
        lifetime folders are processed. Prefixes are matched against the leading
        number before the first underscore in the folder name (e.g. "02_...").
    """
    out: list[Path] = []
    prefix_set = set(only_prefixes) if only_prefixes else None

    for child in sorted(root.iterdir()):
        if not child.is_dir():
            continue
        # Skip RGB output folders (legend.png marker) regardless of naming.
        if (child / "legend.png").is_file():
            continue
        if not is_lifetime_folder_name(child.name, intensity_folder, rgb_suffix):
            continue

        if prefix_set is not None:
            m = re.match(r"^(\d+)_", child.name)
            if not m:
                continue
            pref = m.group(1).zfill(2)
            if pref not in prefix_set:
                continue

        if not _has_tifs_recursively(child):
            continue
        out.append(child)

    return out
def process_intensity_root_8bit(root: Path, args: argparse.Namespace) -> None:
    """Convert raw intensity TIFFs to per-image 8-bit grayscale TIFFs.

    This is only used in export-root mode. Files are written into a parallel
    folder next to the raw intensity folder, named like:

        <intensity_folder><intensity8_suffix>_p1to99/

    where the tag reflects the chosen saturation percentiles.

    By default, this step is *idempotent*: if the 8-bit intensity output files
    already exist, they are not rewritten (use --overwrite-intensity-8bit to force).
    """
    if not args.export_intensity_8bit:
        return

    intensity_root = root / args.intensity_folder
    if not intensity_root.is_dir():
        if args.verbose:
            print(f"[WARN] Intensity folder not found (skipping intensity export): {intensity_root}")
        return

    sat = float(args.intensity_sat)
    if sat < 0 or sat >= 50:
        raise SystemExit("--intensity-sat must be in [0, 50). Example: 1 → 1–99% scaling")

    outfolder = root / f"{args.intensity_folder}{args.intensity8_suffix}{intensity_sat_tag(sat)}"

    tif_paths = [p for p in intensity_root.rglob("*") if p.is_file() and p.suffix.lower() in {".tif", ".tiff"}]
    tif_paths.sort()

    if not tif_paths:
        if args.verbose:
            print(f"[WARN] No TIFFs found under intensity folder: {intensity_root}")
        return

    # If all expected output files already exist, skip the entire intensity export step
    # (useful when re-running the script to regenerate lifetime RGB with different scales).
    if outfolder.is_dir() and not args.overwrite_intensity_8bit:
        any_missing = False
        for in_path in tif_paths:
            rel = in_path.relative_to(intensity_root)
            out_path = outfolder / rel
            if out_path.suffix.lower() not in (".tif", ".tiff"):
                out_path = out_path.with_suffix(".tif")
            if not out_path.exists():
                any_missing = True
                break

        if not any_missing:
            print(f"\n[{args.intensity_folder}] 8-bit intensity outputs already exist → {outfolder.name} (skipping)")
            return

    outfolder.mkdir(parents=True, exist_ok=True)

    print(f"\n[{args.intensity_folder}] {len(tif_paths)} raw intensity TIFFs → {outfolder.name}")

    n_written = 0
    n_skipped = 0
    for in_path in tif_paths:
        rel = in_path.relative_to(intensity_root)
        out_path = outfolder / rel
        if out_path.suffix.lower() not in (".tif", ".tiff"):
            out_path = out_path.with_suffix(".tif")

        if out_path.exists() and not args.overwrite_intensity_8bit:
            n_skipped += 1
            if args.verbose:
                print(f"  [SKIP] Exists: {out_path}")
            continue

        img = tifffile.imread(str(in_path))
        img_u8, vmin, vmax = compute_intensity_8bit(img, sat=sat)

        write_gray_u8_tiff(out_path, img_u8, compress=args.compress, compression_level=args.compression_level)
        n_written += 1
        if args.verbose:
            print(f"  Saved: {out_path}  (vmin={vmin:.3g}, vmax={vmax:.3g})")

    print(f"  Done. Wrote {n_written} 8-bit intensity TIFF(s), skipped {n_skipped}.")

def process_export_root(
    root: Path,
    lifetime_folders: list[Path],
    args: argparse.Namespace,
) -> None:
    intensity_root = root / args.intensity_folder
    if not intensity_root.is_dir():
        raise SystemExit(f"Intensity folder not found: {intensity_root}")

    # Also export 8-bit grayscale intensity images (optional, enabled by default).
    process_intensity_root_8bit(root, args)

    if not lifetime_folders:
        if getattr(args, "only_prefix", None):
            raise SystemExit(
                f"No lifetime folders found under: {root} matching --only-prefix {args.only_prefix}. "
                "Tip: omit --only-prefix to process all lifetime folders."
            )
        raise SystemExit(f"No lifetime folders found under: {root}")

    # Tag output folders with either fixed bounds or span-based scaling so multiple
    # visualization settings can be kept side-by-side.
    scale_tag = flim_vis_tag(args)

    use_span = args.flim_span is not None
    span = float(args.flim_span) if use_span else None

    total = 0

    # Summary rows for CSV (per-image scaling used).
    summary_rows: list[dict[str, object]] = []
    for lf in lifetime_folders:
        outfolder = lf.parent / f"{lf.name}{args.suffix}{scale_tag}"
        outfolder.mkdir(parents=True, exist_ok=True)

        tif_paths = [p for p in lf.rglob("*") if p.is_file() and p.suffix.lower() in {".tif", ".tiff"}]
        tif_paths.sort()

        if not tif_paths:
            continue

        print(f"\n[{lf.name}] {len(tif_paths)} lifetime TIFFs → {outfolder.name}")

        # Collect per-image medians for a representative legend in span mode.
        medians: list[float] = []

        for flim_path in tif_paths:
            rel = flim_path.relative_to(lf)
            photon_path = intensity_root / rel

            flim_gray = tifffile.imread(str(flim_path))
            flim_gray = ensure_2d(np.asarray(flim_gray))

            if photon_path.exists():
                photon = tifffile.imread(str(photon_path))
                photon = ensure_2d(np.asarray(photon))
            else:
                photon = np.ones_like(flim_gray, dtype=np.float32)
                if args.verbose:
                    print(f"  [WARN] Missing intensity: {photon_path} → using constant brightness")

            pc_min, pc_max, flim_ns_min, flim_ns_max, flim_ns, flim_median, median_meta = compute_used_ranges(
                photon=photon,
                flim_gray=flim_gray,
                flim_unit_scale=args.flim_unit_scale,
                pc_min=args.pc_min,
                pc_max=args.pc_max,
                pc_pmin=args.pc_pmin,
                pc_pmax=args.pc_pmax,
                flim_ns_min=args.flim_ns_min,
                flim_ns_max=args.flim_ns_max,
                flim_span=args.flim_span,
                median_intensity_gate=(not args.no_median_intensity_gate),
                median_intensity_min=args.median_intensity_min,
                median_intensity_pmin=args.median_intensity_pmin,
            )
            if flim_median is not None:
                medians.append(float(flim_median))

            # Valid lifetime pixels exclude zero-coded background.
            valid_mask = np.isfinite(flim_ns) & (flim_gray > 0)
            n_valid = int(np.count_nonzero(valid_mask))
            n_total = int(flim_gray.size)

            # Record the actual lifetime scaling used for this image.
            summary_rows.append(
                {
                    "lifetime_folder": lf.name,
                    "relative_path": str(rel.as_posix()),
                    "rgb_output_folder": outfolder.name,
                    "rgb_output_path": str((outfolder / rel).as_posix()),
                    "intensity_found": bool(photon_path.exists()),
                    "mode": "median_span" if use_span else "fixed",
                    "requested_span_ns": float(span) if use_span and span is not None else "",
                    "median_ns": float(flim_median) if flim_median is not None else "",
                    "median_intensity_gate": median_meta.get("median_intensity_gate", ""),
                    "median_intensity_min": median_meta.get("median_intensity_min", ""),
                    "median_intensity_pmin": median_meta.get("median_intensity_pmin", ""),
                    "median_intensity_threshold": median_meta.get("median_intensity_threshold", ""),
                    "median_candidate_pixels": median_meta.get("median_candidate_pixels", n_valid),
                    "median_pixels_used": median_meta.get("median_pixels_used", n_valid),
                    "median_used_fraction": median_meta.get("median_used_fraction", ""),
                    "median_gate_fallback": median_meta.get("median_gate_fallback", ""),
                    "vmin_ns": float(flim_ns_min),
                    "vmax_ns": float(flim_ns_max),
                    "actual_span_ns": float(flim_ns_max - flim_ns_min),
                    "valid_pixels": n_valid,
                    "total_pixels": n_total,
                    "valid_fraction": (n_valid / n_total) if n_total else "",
                }
            )

            rgb = apply_flim_coloring(
                photon=photon,
                flim_ns=flim_ns,
                pc_min=pc_min,
                pc_max=pc_max,
                flim_ns_min=flim_ns_min,
                flim_ns_max=flim_ns_max,
                invalid_mask=(flim_gray <= 0),
            )

            out_path = outfolder / rel
            if out_path.suffix.lower() not in (".tif", ".tiff"):
                out_path = out_path.with_suffix(".tif")

            if out_path.exists() and not args.overwrite:
                if args.verbose:
                    print(f"  [SKIP] Exists: {out_path}")
                continue

            write_rgb_tiff(out_path, rgb, compress=args.compress, compression_level=args.compression_level)
            total += 1
            if args.verbose:
                if use_span and flim_median is not None:
                    print(f"  Saved: {out_path}  (median={flim_median:.3g} ns, range={flim_ns_min:.3g}–{flim_ns_max:.3g} ns)")
                else:
                    print(f"  Saved: {out_path}")

        # One legend per output folder:
        if use_span and span is not None:
            if medians:
                center = float(np.median(np.asarray(medians, dtype=np.float32)))
            else:
                center = float((args.flim_ns_min + args.flim_ns_max) / 2.0)
            legend_vmin = max(0.0, center - span / 2.0)
            legend_vmax = center + span / 2.0
            legend_title = f"FLIM lifetime (ns) (median-centered; span={span:g} ns)"
        else:
            legend_vmin = float(args.flim_ns_min)
            legend_vmax = float(args.flim_ns_max)
            legend_title = "FLIM lifetime (ns)"

        legend_png = outfolder / "legend.png"
        save_legend_png(legend_png, vmin=legend_vmin, vmax=legend_vmax, title=legend_title)
        print(f"  Legend: {legend_png}")

    print(f"\nDone. Wrote {total} RGB TIFF(s).")

    # Write one consolidated summary CSV per run (in export-root mode).
    if summary_rows:
        csv_name = f"FLIM_RGB_scale_summary{scale_tag}.csv"
        csv_path = root / csv_name
        fieldnames = list(summary_rows[0].keys())
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(summary_rows)
        print(f"Summary CSV: {csv_path}")


# ------------------------ CLI ----------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="FLIM to RGB visualization (hue=lifetime, intensity=photon count). "
                    "Supports export-root mode for LIF-exported folder structures."
    )

    # Positional inputs: either FOLDER or two FILES
    p.add_argument("inputs", nargs="+", help="Either an export root folder, a legacy folder, or two files: <photon> <flim>")

    # Scaling options (with shorthands)
    p.add_argument(
        "-s", "--flim-unit-scale",
        type=float,
        default=1000.0,
        help="Gray values per ns (default: 1000 → 2400 gray = 2.400 ns)",
    )

    p.add_argument("-m", "--pc-min", type=float, default=None, help="Photon min (overrides percentile)")
    p.add_argument("-M", "--pc-max", type=float, default=None, help="Photon max (overrides percentile)")
    p.add_argument("-q", "--pc-pmin", type=float, default=1.0, help="Photon lower percentile (default: 1)")
    p.add_argument("-Q", "--pc-pmax", type=float, default=99.0, help="Photon upper percentile (default: 99)")

    # Fixed FLIM range with explicit overrides
    p.add_argument("-l", "--flim-ns-min", type=float, default=1.5, help="FLIM lower bound in ns (default: 1.5)")
    p.add_argument("-L", "--flim-ns-max", type=float, default=3.0, help="FLIM upper bound in ns (default: 3.0)")


    p.add_argument(
        "--flim-span",
        type=float,
        default=None,
        help="Alternative to -l/-L: scale EACH image to a range centered on its median lifetime. "
             "Range = [median - span/2, median + span/2] in ns. "
             "Example: --flim-span 1.5 and median 1.3 → 0.55–2.05 ns.",
    )

    # Span-median computation: optional intensity gating to exclude dim/noisy pixels
    p.add_argument(
        "--median-intensity-min",
        type=float,
        default=None,
        help="In --flim-span mode only: absolute photon-count threshold used to exclude dim pixels when computing the per-image median lifetime. "
             "If provided, overrides --median-intensity-pmin.",
    )
    p.add_argument(
        "--median-intensity-pmin",
        type=float,
        default=10.0,
        help="In --flim-span mode only: percentile (0–100) of photon counts among valid-lifetime pixels (lifetime>0) used as an intensity threshold for median computation. "
             "Pixels dimmer than this threshold are excluded from the median. Default: 10.",
    )
    p.add_argument(
        "--no-median-intensity-gate",
        action="store_true",
        help="In --flim-span mode only: disable intensity gating for median computation (use all lifetime>0 pixels).",
    )


    # Output controls (legacy two-file and legacy folder mode)
    p.add_argument("-o", "--output", type=Path, default=None, help="Output RGB TIFF path (only in two-file mode). If omitted, inferred.")
    p.add_argument("-O", "--outdir", type=Path, default=None, help="Output directory (legacy folder mode). If omitted, inferred.")
    p.add_argument("-c", "--compress", type=str, default="deflate", help="TIFF compression (default: deflate)")
    p.add_argument("-C", "--compression-level", type=int, default=6, help="Compression level if applicable (default: 6)")

    # Export-root mode controls
    p.add_argument(
        "--intensity-folder",
        type=str,
        default="06_Intensity_raw",
        help='Folder name (under export root) that contains raw intensity TIFFs (default: "06_Intensity_raw")',
    )
    p.add_argument(
        "--suffix",
        type=str,
        default="_FLIM_RGB",
        help='Suffix appended to each lifetime folder name to create the parallel RGB folder (default: "_FLIM_RGB")',
    )
    p.add_argument(
        "--only-prefix",
        action="append",
        default=None,
        help="In export-root mode: only process lifetime folders whose name starts with the given numeric prefix "
             "(e.g. 02 for '02_Phasor_phase_lifetime'). You can provide this option multiple times, or pass a comma-separated list.",
    )


    # Raw intensity → 8-bit export (export-root mode)
    g_int = p.add_mutually_exclusive_group()
    g_int.add_argument(
        "--intensity-8bit",
        dest="export_intensity_8bit",
        action="store_true",
        help="Also export 8-bit intensity TIFFs (default).",
    )
    g_int.add_argument(
        "--no-intensity-8bit",
        dest="export_intensity_8bit",
        action="store_false",
        help="Do not export 8-bit intensity TIFFs.",
    )
    p.set_defaults(export_intensity_8bit=True)
    p.add_argument(
        "--intensity-sat",
        type=float,
        default=1.0,
        help="Intensity saturation percentile S used as low/high percentiles (S, 100-S) for 8-bit scaling. "
             "Default 1.0 → 1–99%%.",
    )
    p.add_argument(
        "--intensity8-suffix",
        type=str,
        default="_8bit",
        help="Suffix appended to the raw intensity folder name for 8-bit outputs (default: '_8bit').",
    )
    p.add_argument(
        "--overwrite-intensity-8bit",
        action="store_true",
        help="Overwrite existing 8-bit intensity TIFFs if present. Default: skip files that already exist.",
    )

    # Overwrite behavior: default is overwrite (use --no-overwrite to skip).
    g_over = p.add_mutually_exclusive_group()
    g_over.add_argument(
        "--overwrite",
        dest="overwrite",
        action="store_true",
        help="Overwrite existing RGB TIFFs (default).",
    )
    g_over.add_argument(
        "--no-overwrite",
        dest="overwrite",
        action="store_false",
        help="Do not overwrite existing RGB TIFFs (skip if present).",
    )
    p.set_defaults(overwrite=True)
    p.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose logging (prints every file saved/skipped).",
    )

    args = p.parse_args()

    # Normalize --only-prefix into a sorted list of two-digit strings (e.g. ["01", "02"]).
    # This option is only used in export-root mode.
    if args.only_prefix is not None:
        prefixes: list[str] = []
        for token in args.only_prefix:
            if token is None:
                continue
            for part in str(token).split(","):
                part = part.strip()
                if not part:
                    continue
                m = re.match(r"^(\d+)", part)
                if not m:
                    raise SystemExit(f"--only-prefix must start with digits (e.g. 02), got: {part!r}")
                prefixes.append(m.group(1).zfill(2))
        args.only_prefix = sorted(set(prefixes)) if prefixes else None

    return args


# ------------------------ Legacy-mode path inference -----------------------

def infer_outdir_for_folder(folder: Path, explicit_outdir: Path | None, args: argparse.Namespace) -> Path:
    """Infer output directory for legacy folder-pairing mode.

    If not explicitly provided, name it after the input folder plus the
    configured suffix and FLIM range tag.
    """
    if explicit_outdir is not None:
        return explicit_outdir
    return folder.parent / f"{folder.name}{args.suffix}{flim_vis_tag(args)}"


def infer_paths_for_pair(photon: Path, explicit_output: Path | None, args: argparse.Namespace) -> tuple[Path, Path]:
    """Return (outdir, outfile_path)."""
    if explicit_output is not None:
        outdir = explicit_output.parent
        outfile = explicit_output
    else:
        tag = flim_vis_tag(args)
        outdir = photon.parent / f"{photon.stem}{args.suffix}{tag}"
        outfile = outdir / f"{photon.stem}{args.suffix}{tag}.tif"
    return outdir, outfile


def process_pair(
    photon_path: Path,
    flim_path: Path,
    outdir: Path,
    outfile: Path,
    args: argparse.Namespace,
) -> tuple[float, float]:
    photon = tifffile.imread(str(photon_path))
    flim_gray = tifffile.imread(str(flim_path))

    pc_min, pc_max, flim_ns_min, flim_ns_max, flim_ns, _flim_median, _median_meta = compute_used_ranges(
        photon=photon,
        flim_gray=flim_gray,
        flim_unit_scale=args.flim_unit_scale,
        pc_min=args.pc_min,
        pc_max=args.pc_max,
        pc_pmin=args.pc_pmin,
        pc_pmax=args.pc_pmax,
        flim_ns_min=args.flim_ns_min,
        flim_ns_max=args.flim_ns_max,
    )

    rgb = apply_flim_coloring(
        photon=photon,
        flim_ns=flim_ns,
        pc_min=pc_min,
        pc_max=pc_max,
        flim_ns_min=flim_ns_min,
        flim_ns_max=flim_ns_max,
        invalid_mask=(ensure_2d(np.asarray(flim_gray)) <= 0),
    )

    # Save TIFF
    write_rgb_tiff(outfile, rgb, compress=args.compress, compression_level=args.compression_level)
    print(f"Saved: {outfile}")

    return flim_ns_min, flim_ns_max


# ------------------------ Main ---------------------------------------------

def main() -> int:
    args = parse_args()
    inputs = [Path(x) for x in args.inputs]

    if len(inputs) == 1 and inputs[0].is_dir():
        folder = inputs[0]

        # If user points directly at an RGB output folder created by this script, skip it.
        if looks_like_rgb_output_folder(folder, args.suffix):
            print(f"Input folder looks like an RGB output folder ({folder.name}); skipping.")
            return 0

        # If user points at a lifetime folder directly, but its parent looks like an export root,
        # process only that lifetime folder.
        if is_lifetime_folder_name(folder.name, args.intensity_folder, args.suffix) and detect_export_root(folder.parent, args.intensity_folder, args.suffix):
            root = folder.parent
            process_export_root(root, [folder], args)
            return 0

        # Export-root mode (preferred)
        if detect_export_root(folder, args.intensity_folder, args.suffix):
            lifetime_folders = find_lifetime_folders(folder, args.intensity_folder, args.suffix, only_prefixes=args.only_prefix)
            process_export_root(folder, lifetime_folders, args)
            return 0

        # Legacy "folder-pairing" mode
        outdir = infer_outdir_for_folder(folder, args.outdir, args)
        outdir.mkdir(parents=True, exist_ok=True)
        pairs = find_pairs_in_folder(folder)
        if not pairs:
            raise SystemExit(
                "Folder mode: no ch-based pairs found, and folder did not look like an export root.\n"
                "Provide an export-root folder (containing 06_Intensity_raw/ and lifetime folders),\n"
                "or a legacy folder with files named ...ch0..., ...ch1..., etc."
            )

        # Use the configured FLIM min/max for the single folder legend
        legend_vmin, legend_vmax = args.flim_ns_min, args.flim_ns_max

        for photon_p, flim_p, prefix, pair_idx in pairs:
            base = prefix if prefix else photon_p.stem
            outfile = outdir / f"{base}_FLIM_RGB_ch{pair_idx}.tif"
            process_pair(photon_p, flim_p, outdir, outfile, args)

        legend_png = outdir / f"{folder.name}_FLIM_RGB_legend.png"
        save_legend_png(legend_png, vmin=legend_vmin, vmax=legend_vmax)
        print(f"Saved legend: {legend_png}")
        return 0

    elif len(inputs) == 2 and inputs[0].is_file() and inputs[1].is_file():
        # Two-file mode
        photon_p, flim_p = inputs
        outdir, outfile = infer_paths_for_pair(photon_p, args.output, args)
        outdir.mkdir(parents=True, exist_ok=True)
        flim_min, flim_max = process_pair(photon_p, flim_p, outdir, outfile, args)

        legend_png = outfile.with_suffix("")
        legend_png = legend_png.parent / f"{legend_png.stem}_legend.png"
        save_legend_png(legend_png, vmin=flim_min, vmax=flim_max)
        print(f"Saved legend: {legend_png}")
        return 0

    else:
        raise SystemExit("Provide either a folder, or two files: <photon> <flim>.")


if __name__ == "__main__":
    raise SystemExit(main())
