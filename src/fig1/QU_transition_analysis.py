#!/usr/bin/env python3
"""Estimate Q-to-U transition rates from single-molecule dwell-time data.

Reads a table of per-molecule dwell times (xlsx or csv), groups by force, and
reports the maximum-likelihood exponential rate with exact Poisson confidence
intervals. Right-censored observations contribute exposure time but no events.

The header row and the column roles are inferred, so title blocks, notes, and
illustration rows above the header are tolerated, and column names need not
match exactly. The inferred layout is always reported; override any part of it
with --header-row / --force-col / --time-col / --event-col.

Usage:
    python qu_transition_analysis.py input.xlsx
    python qu_transition_analysis.py data.csv -o summary.csv
    python qu_transition_analysis.py messy.xlsx --sheet "Run 3" --dry-run
    python qu_transition_analysis.py data.csv --header-row 4 --event-col censored
    python qu_transition_analysis.py runs.xlsx --event-values QU --censor-values END,BL,BR
"""

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

# ----------------------------------------------------------------------------
# CONFIG
# ----------------------------------------------------------------------------

EXCEL_SUFFIXES = {".xlsx", ".xlsm", ".xltx", ".xls"}
CSV_SUFFIXES = {".csv", ".tsv", ".tab", ".txt"}
OUTPUT_SUFFIX = "_QU_transition_summary"

# Used when the input has no identifier column at all.
SINGLE_MOLECULE_ID = "molecule_1"

# Canonical role -> name patterns, checked in order. Names are normalized to
# lowercase with runs of non-alphanumerics collapsed to "_" before matching.
NAME_PATTERNS = {
    "force_pN": [
        "force_pn", "force", "f_pn", "tension", "load", "applied_force",
        "force_level", "pn", "f",
    ],
    "time_Q_s": [
        "time_q_s", "time_q", "t_q", "dwell_time_s", "dwell_time", "dwell_s",
        "dwell", "residence_time", "waiting_time", "survival_time",
        "lifetime_s", "lifetime", "duration_s", "duration", "time_s", "time",
        "t_s", "t",
    ],
    "event": [
        "event", "event_flag", "is_event", "events", "status", "transition",
        "transitioned", "uncensored", "observed", "censored", "is_censored",
        "censor", "right_censored",
    ],
    "molecule_id": [
        "molecule_id", "molecule", "mol_id", "mol", "trace_id", "trace",
        "particle_id", "particle", "roi_id", "roi", "molecule_name", "id",
        "name",
    ],
}

# Names whose 1 means "still in state Q" rather than "transitioned".
CENSOR_SEMANTICS = re.compile(r"censor")
NOT_CENSOR_SEMANTICS = re.compile(r"un_?censor|non_?censor|not_?censor")

# Categorical outcome codes, for files that record why an observation ended
# rather than a 0/1 flag. 1 means the Q->U transition was seen; 0 means the
# observation was right-censored for some other reason. Edit this to match your
# lab's conventions, or pass --event-values / --censor-values per run. A code
# that appears in the data but not here is never silently censored: those rows
# are dropped and reported.
EVENT_CODES = {
    "qu": 1,    # quenched -> unquenched state transition
    "end": 0,   # end of planned session
    "bl": 0,    # photobleached
    "br": 0,    # tether break
}

EVENT_TRUE = {"1", "true", "t", "yes", "y", "event", "transitioned",
              "transition", "observed", "uncensored", "complete", "dissociated"}
EVENT_FALSE = {"0", "false", "f", "no", "n", "censored", "censor",
               "right_censored", "no_event", "none", "incomplete", "survived"}

# A cell that is a number possibly followed by a unit ("5 pN", "2.3s").
NUMERIC_WITH_UNIT = re.compile(
    r"^\s*[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?\s*[a-zA-Z\u00b5%/\^\-\.]*\s*$"
)
LEADING_NUMBER = re.compile(r"^\s*[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


# ----------------------------------------------------------------------------
# Reading raw, un-headered tables
# ----------------------------------------------------------------------------

def normalize_name(value):
    text = re.sub(r"[^0-9a-zA-Z]+", "_", str(value).strip().lower())
    return text.strip("_")


def is_blank(value):
    if value is None:
        return True
    if isinstance(value, float) and np.isnan(value):
        return True
    return str(value).strip() == ""


def cell_is_numeric(value):
    if is_blank(value):
        return False
    if isinstance(value, bool):
        return False
    if isinstance(value, (int, float, np.integer, np.floating)):
        return True
    return bool(NUMERIC_WITH_UNIT.match(str(value)))


def coerce_numeric(series):
    """Numeric coercion that also accepts a number with a trailing unit."""
    out = pd.to_numeric(series, errors="coerce")
    missing = out.isna() & series.notna()
    if missing.any():
        salvaged = series[missing].astype(str).apply(
            lambda s: float(LEADING_NUMBER.match(s).group())
            if NUMERIC_WITH_UNIT.match(s) and LEADING_NUMBER.match(s)
            else np.nan
        )
        out.loc[missing] = salvaged
    return out


def coerce_event(series, event_values=None, censor_values=None):
    """Map event flags or outcome codes onto 1 (transition) and 0 (censored).

    Resolution order: values named on the command line, then the EVENT_CODES
    registry, then plain numbers, then generic spellings such as yes/no. An
    unrecognized value becomes NaN so that the row is dropped and reported,
    never silently counted as censored.
    """
    named_event = {normalize_name(v) for v in (event_values or [])}
    named_censor = {normalize_name(v) for v in (censor_values or [])}

    def one(value):
        if is_blank(value):
            return np.nan
        if isinstance(value, bool):
            return float(value)
        token = normalize_name(value)
        if token in named_event:
            return 1.0
        if token in named_censor:
            return 0.0
        if token in EVENT_CODES:
            return float(EVENT_CODES[token])
        if isinstance(value, (int, float, np.integer, np.floating)):
            return float(value)
        if token in EVENT_TRUE:
            return 1.0
        if token in EVENT_FALSE:
            return 0.0
        return pd.to_numeric(value, errors="coerce")

    return series.apply(one).astype(float)


def read_raw(path, sheet=0, sep=None):
    """Return the sheet as an all-object DataFrame with no header applied."""
    suffix = path.suffix.lower()

    if suffix in EXCEL_SUFFIXES:
        if isinstance(sheet, str) and sheet.isdigit():
            sheet = int(sheet)
        raw = pd.read_excel(path, sheet_name=sheet, header=None, dtype=object)

    elif suffix in CSV_SUFFIXES:
        if sep is None:
            sep = "\t" if suffix in {".tsv", ".tab"} else ","
        with open(path, newline="", encoding="utf-8-sig") as handle:
            rows = list(csv.reader(handle, delimiter=sep))
        if not rows:
            raise SystemExit("Input file is empty: %s" % path)
        width = max(len(r) for r in rows)
        rows = [r + [""] * (width - len(r)) for r in rows]  # ragged rows are fine
        raw = pd.DataFrame(rows, dtype=object)

    else:
        raise SystemExit(
            "Unrecognized input extension '%s'. Use one of: %s"
            % (suffix, ", ".join(sorted(EXCEL_SUFFIXES | CSV_SUFFIXES)))
        )

    # Row indices are left as they are in the file so that --header-row means
    # the same thing the spreadsheet shows. Blank columns carry no data ever.
    keep = [c for c in raw.columns if not all(is_blank(v) for v in raw[c])]
    raw = raw[keep]
    if raw.empty or not len(raw.columns):
        raise SystemExit("No data found in %s" % path)
    raw.columns = range(raw.shape[1])
    return raw


# ----------------------------------------------------------------------------
# Header inference
# ----------------------------------------------------------------------------

def name_role(normalized):
    """Canonical role for a normalized column name, or None."""
    for role, patterns in NAME_PATTERNS.items():
        if normalized in patterns:
            return role
    for role, patterns in NAME_PATTERNS.items():
        for pattern in patterns:
            if len(pattern) > 2 and pattern in normalized:
                return role
    return None


def numeric_fraction(row):
    cells = [v for v in row if not is_blank(v)]
    if not cells:
        return 0.0
    return sum(cell_is_numeric(v) for v in cells) / len(cells)


def score_header_row(raw, i):
    row = raw.iloc[i]
    width = len(row)
    filled = [v for v in row if not is_blank(v)]
    if len(filled) < max(2, 0.5 * width):
        return None  # title lines and stray notes, not a header

    roles = {name_role(normalize_name(v)) for v in filled}
    roles.discard(None)
    textish = 1.0 - numeric_fraction(row)

    below = raw.iloc[i + 1: i + 6]
    below_numeric = (
        float(np.mean([numeric_fraction(below.iloc[j]) for j in range(len(below))]))
        if len(below) else 0.0
    )

    score = 5 * len(roles) + 3 * textish + 2 * below_numeric + len(filled) / width
    return score, len(roles), textish


def detect_header(raw, max_scan=25):
    """Return the index of the most header-like row, or None if headerless."""
    best = None
    for i in range(min(max_scan, len(raw))):
        scored = score_header_row(raw, i)
        if scored is None:
            continue
        score, n_roles, textish = scored
        if best is None or score > best[0]:
            best = (score, i, n_roles, textish)

    if best is None:
        return None
    _, index, n_roles, textish = best
    if n_roles == 0 and textish < 0.5:
        return None  # looks like data all the way up; treat as headerless
    return index


def drop_blank_rows(data):
    keep = ~data.apply(lambda r: all(is_blank(v) for v in r), axis=1)
    return data.loc[keep].reset_index(drop=True)


def apply_header(raw, header_row):
    if header_row is None:
        data = raw.copy()
        data.columns = ["col_%d" % i for i in range(raw.shape[1])]
        return drop_blank_rows(data)

    names = []
    for i, value in enumerate(raw.iloc[header_row]):
        name = str(value).strip() if not is_blank(value) else "col_%d" % i
        while name in names:
            name += "_dup"
        names.append(name)

    data = raw.iloc[header_row + 1:].copy()
    data.columns = names
    return drop_blank_rows(data)


# ----------------------------------------------------------------------------
# Column inference
# ----------------------------------------------------------------------------

def match_column(data, requested):
    """Resolve a user-supplied column reference: exact, normalized, or index."""
    if requested in data.columns:
        return requested
    target = normalize_name(requested)
    for column in data.columns:
        if normalize_name(column) == target:
            return column
    if re.fullmatch(r"\d+", str(requested)):
        index = int(requested)
        if 0 <= index < data.shape[1]:
            return data.columns[index]
    raise SystemExit(
        "Column '%s' not found. Available: %s"
        % (requested, ", ".join(map(str, data.columns)))
    )


def profile_column(series):
    numeric = coerce_numeric(series)
    non_null = numeric.dropna()
    values = series[~series.apply(is_blank)]
    return {
        "numeric_frac": len(non_null) / len(values) if len(values) else 0.0,
        "n": len(non_null),
        "unique": non_null.nunique(),
        "unique_frac": non_null.nunique() / len(non_null) if len(non_null) else 0.0,
        "all_positive": bool(len(non_null)) and bool((non_null > 0).all()),
        "binary": bool(len(non_null)) and set(non_null.unique()) <= {0.0, 1.0},
        "frac_fractional": (
            float((non_null % 1 != 0).mean()) if len(non_null) else 0.0
        ),
    }


def resolve_columns(data, overrides, verbose=True):
    """Map columns onto the canonical roles by name, then by content."""
    mapping = {}
    how = {}

    for role, requested in overrides.items():
        if requested:
            mapping[role] = match_column(data, requested)
            how[role] = "specified"

    # Pass 1: column names.
    for column in data.columns:
        role = name_role(normalize_name(column))
        if role and role not in mapping and column not in mapping.values():
            mapping[role] = column
            how[role] = "name"

    # Pass 2: content of the columns nothing has claimed yet.
    unclaimed = [c for c in data.columns if c not in mapping.values()]
    profiles = {c: profile_column(data[c]) for c in unclaimed}
    numeric_cols = [
        c for c in unclaimed
        if profiles[c]["numeric_frac"] >= 0.7 and profiles[c]["n"] > 0
    ]

    if "event" not in mapping:
        binary = [c for c in numeric_cols if profiles[c]["binary"]]
        if binary:
            pick = binary[-1]  # event flags conventionally sit last
            mapping["event"] = pick
            how["event"] = "content (only 0/1 values)"
            numeric_cols.remove(pick)

    need_both = "force_pN" not in mapping and "time_Q_s" not in mapping
    if need_both and len(numeric_cols) >= 2:
        # The grouping variable repeats; dwell times rarely do.
        by_cardinality = sorted(numeric_cols, key=lambda c: profiles[c]["unique"])
        force, time = by_cardinality[0], by_cardinality[-1]
        mapping["force_pN"] = force
        how["force_pN"] = "content (few repeated values)"
        mapping["time_Q_s"] = time
        how["time_Q_s"] = "content (many distinct positive values)"
        numeric_cols = [c for c in numeric_cols if c not in (force, time)]
    else:
        if "force_pN" not in mapping and numeric_cols:
            pick = min(numeric_cols, key=lambda c: profiles[c]["unique"])
            mapping["force_pN"] = pick
            how["force_pN"] = "content (few repeated values)"
            numeric_cols.remove(pick)
        if "time_Q_s" not in mapping and numeric_cols:
            candidates = [c for c in numeric_cols
                          if profiles[c]["all_positive"]] or numeric_cols
            pick = max(
                candidates,
                key=lambda c: profiles[c]["unique_frac"] + profiles[c]["frac_fractional"],
            )
            mapping["time_Q_s"] = pick
            how["time_Q_s"] = "content (many distinct positive values)"
            numeric_cols.remove(pick)

    missing = [r for r in ("force_pN", "time_Q_s", "event") if r not in mapping]
    if missing:
        raise SystemExit(
            "Could not infer column(s): %s\nColumns available: %s\n"
            "Name them explicitly with --force-col / --time-col / --event-col, "
            "or check --header-row."
            % (", ".join(missing), ", ".join(map(str, data.columns)))
        )

    # Does the event column count events, or count censoring?
    event_name = normalize_name(mapping["event"])
    invert = bool(
        CENSOR_SEMANTICS.search(event_name)
        and not NOT_CENSOR_SEMANTICS.search(event_name)
    )

    if verbose:
        print("Column mapping:", file=sys.stderr)
        for role in ("force_pN", "time_Q_s", "event", "molecule_id"):
            if role in mapping:
                print("  %-12s <- %-24s [%s]"
                      % (role, mapping[role], how[role]), file=sys.stderr)
            else:
                print("  %-12s <- (absent, all rows treated as one molecule)"
                      % role, file=sys.stderr)
        if invert:
            print("  event column reads as censoring, so inverting 0/1",
                  file=sys.stderr)

    return mapping, invert


# ----------------------------------------------------------------------------
# Cleaning and fitting
# ----------------------------------------------------------------------------

def build_frame(data, mapping, invert_event, censored_flag,
                event_values=None, censor_values=None, verbose=True):
    out = pd.DataFrame(index=data.index)
    if "molecule_id" in mapping:
        out["molecule_id"] = data[mapping["molecule_id"]]
    else:
        # No identifier in the file: every row belongs to the same molecule.
        out["molecule_id"] = SINGLE_MOLECULE_ID
    out["force_pN"] = coerce_numeric(data[mapping["force_pN"]])
    out["time_Q_s"] = coerce_numeric(data[mapping["time_Q_s"]])

    raw_event = data[mapping["event"]]
    out["event"] = coerce_event(raw_event, event_values, censor_values)

    if invert_event or censored_flag:
        out["event"] = 1 - out["event"]

    if verbose:
        report_event_decoding(raw_event, out["event"])

    n_total = len(out)
    bad_force = out["force_pN"].isna()
    bad_time = out["time_Q_s"].isna() | (out["time_Q_s"] <= 0)
    bad_event = out["event"].isna() | ~out["event"].isin([0, 1])

    valid = out.loc[~(bad_force | bad_time | bad_event)].copy()
    valid["event"] = valid["event"].astype(int)

    if verbose and len(valid) < n_total:
        for label, mask, column in (
            ("force not a number", bad_force, mapping["force_pN"]),
            ("dwell time not a positive number", bad_time, mapping["time_Q_s"]),
            ("outcome not recognized", bad_event, mapping["event"]),
        ):
            if mask.any():
                offenders = sorted(
                    {str(v).strip() for v in raw_event[mask]} if label.startswith("outcome")
                    else {str(v).strip() for v in data[column][mask]}
                )[:8]
                print("  %d row(s) dropped, %s in '%s': %s"
                      % (mask.sum(), label, column, ", ".join(offenders) or "blank"),
                      file=sys.stderr)

    return valid, n_total


def report_event_decoding(raw_event, coded):
    """Show which outcome codes were read as events and which as censored."""
    pairs = pd.DataFrame({"code": raw_event.astype(str).str.strip(), "coded": coded})
    lines = []
    for value, label in ((1.0, "event"), (0.0, "censored")):
        subset = pairs.loc[pairs["coded"] == value, "code"].value_counts()
        if len(subset):
            lines.append("%s -> %s (%d)"
                         % (", ".join("%s x%d" % (k, v) for k, v in subset.items()),
                            label, int(subset.sum())))
    unknown = pairs.loc[pairs["coded"].isna(), "code"].value_counts()
    if len(unknown):
        lines.append("unrecognized: %s"
                     % ", ".join("%s x%d" % (k, v) for k, v in unknown.items()))
    if lines:
        print("Outcome decoding: " + "; ".join(lines), file=sys.stderr)


def fit_group(force, g, alpha):
    N = len(g)
    N_ids = int(g["molecule_id"].nunique())
    n = int(g["event"].sum())
    N_censored = N - n
    T_exp = float(g["time_Q_s"].sum())

    if n > 0:
        k_hat = n / T_exp
        tau_hat = T_exp / n

        lambda_low = 0.5 * chi2.ppf(alpha / 2, 2 * n)
        lambda_high = 0.5 * chi2.ppf(1 - alpha / 2, 2 * (n + 1))

        k_low = lambda_low / T_exp
        k_high = lambda_high / T_exp

        tau_low = 1 / k_high
        tau_high = 1 / k_low

        tau_lower_bound = np.nan
        zero_event_flag = False

    else:
        k_hat = np.nan
        tau_hat = np.nan

        k_low = 0.0
        k_high = -np.log(alpha) / T_exp  # one-sided upper bound

        tau_low = np.nan
        tau_high = np.nan

        tau_lower_bound = T_exp / (-np.log(alpha))
        zero_event_flag = True

    ci = int(round(100 * (1 - alpha)))

    return {
        "force_pN": force,
        "N_intervals": N,
        "N_molecule_ids": N_ids,
        "N_events": n,
        "N_censored": N_censored,
        "T_exp_s": T_exp,
        "k_QU_s_inv": k_hat,
        "k_CI%d_low_s_inv" % ci: k_low,
        "k_CI%d_high_s_inv" % ci: k_high,
        "tau_Q_s": tau_hat,
        "tau_CI%d_low_s" % ci: tau_low,
        "tau_CI%d_high_s" % ci: tau_high,
        "tau_lower_bound_s": tau_lower_bound,
        "zero_event_flag": zero_event_flag,
    }


# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------

def default_output(input_path):
    suffix = input_path.suffix.lower()
    if suffix in EXCEL_SUFFIXES:
        out_suffix = ".xlsx"
    else:
        out_suffix = suffix if suffix in CSV_SUFFIXES else ".csv"
    return input_path.with_name(input_path.stem + OUTPUT_SUFFIX + out_suffix)


def write_table(df, path):
    suffix = path.suffix.lower()
    if suffix in EXCEL_SUFFIXES:
        df.to_excel(path, index=False)
    elif suffix in CSV_SUFFIXES:
        sep = "\t" if suffix in {".tsv", ".tab"} else ","
        df.to_csv(path, sep=sep, index=False)
    else:
        raise SystemExit(
            "Unrecognized output extension '%s'. Use one of: %s"
            % (suffix, ", ".join(sorted(EXCEL_SUFFIXES | CSV_SUFFIXES)))
        )


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Estimate Q-to-U transition rates per force from dwell-time data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input", type=Path,
                   help="Input table (.xlsx/.xlsm/.xls or .csv/.tsv/.txt).")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output file; format follows its extension. "
                        "Default: <input_stem>%s<input_ext>." % OUTPUT_SUFFIX)
    p.add_argument("-a", "--alpha", type=float, default=0.05,
                   help="Significance level for confidence intervals.")
    p.add_argument("-s", "--sheet", default=0,
                   help="Sheet name or 0-based index (Excel input only).")
    p.add_argument("--sep", default=None,
                   help="Delimiter for text input. Default: tab for .tsv/.tab, "
                        "comma otherwise.")

    layout = p.add_argument_group("layout overrides")
    layout.add_argument("--header-row", type=int, default=None,
                        help="0-based header row as it sits in the file, i.e. "
                             "the Excel row number minus 1. Default: inferred.")
    layout.add_argument("--no-header", action="store_true",
                        help="Treat every row as data.")
    layout.add_argument("--max-scan", type=int, default=25,
                        help="Rows to scan when inferring the header.")
    layout.add_argument("--force-col", default=None, help="Column for force_pN.")
    layout.add_argument("--time-col", default=None, help="Column for time_Q_s.")
    layout.add_argument("--event-col", default=None, help="Column for event.")
    layout.add_argument("--id-col", default=None, help="Column for molecule_id.")
    layout.add_argument("--event-values", default=None,
                        help="Comma-separated outcome codes meaning the "
                             "transition was observed, e.g. 'QU'. Extends the "
                             "EVENT_CODES registry at the top of this script.")
    layout.add_argument("--censor-values", default=None,
                        help="Comma-separated codes meaning right-censored, "
                             "e.g. 'END,BL,BR'.")
    layout.add_argument("--censored", action="store_true",
                        help="The event column marks censoring, so invert it. "
                             "Inferred from names containing 'censor'.")

    p.add_argument("--dry-run", action="store_true",
                   help="Report the inferred layout and row counts, write nothing.")
    p.add_argument("-q", "--quiet", action="store_true",
                   help="Do not print the summary table to stdout.")
    return p.parse_args(argv)


def split_codes(text):
    if not text:
        return []
    return [part.strip() for part in text.split(",") if part.strip()]


def main(argv=None):
    args = parse_args(argv)

    if not args.input.exists():
        raise SystemExit("Input file not found: %s" % args.input)
    if not 0 < args.alpha < 1:
        raise SystemExit("--alpha must be between 0 and 1, got %s" % args.alpha)

    raw = read_raw(args.input, sheet=args.sheet, sep=args.sep)
    print("Read %d row(s) x %d column(s) from %s"
          % (raw.shape[0], raw.shape[1], args.input), file=sys.stderr)

    if args.no_header:
        header_row = None
    elif args.header_row is not None:
        if not 0 <= args.header_row < len(raw):
            raise SystemExit("--header-row %d is outside the sheet (%d rows)."
                             % (args.header_row, len(raw)))
        header_row = args.header_row
    else:
        header_row = detect_header(raw, max_scan=args.max_scan)

    if header_row is None:
        print("Header: none found, treating every row as data", file=sys.stderr)
    else:
        preview = ", ".join(str(v) for v in raw.iloc[header_row]
                            if not is_blank(v))
        print("Header: file row %d, i.e. row %d in Excel (%s)"
              % (header_row, header_row + 1, preview), file=sys.stderr)
        if header_row > 0:
            print("Skipped %d row(s) above it" % header_row, file=sys.stderr)

    data = apply_header(raw, header_row)

    overrides = {
        "force_pN": args.force_col,
        "time_Q_s": args.time_col,
        "event": args.event_col,
        "molecule_id": args.id_col,
    }
    mapping, invert = resolve_columns(data, overrides)

    event_values = split_codes(args.event_values)
    censor_values = split_codes(args.censor_values)
    valid, n_rows = build_frame(data, mapping, invert, args.censored,
                                event_values, censor_values)
    n_excluded = n_rows - len(valid)
    print("Kept %d of %d data row(s)%s"
          % (len(valid), n_rows,
             "; dropped %d that were not usable numbers" % n_excluded
             if n_excluded else ""),
          file=sys.stderr)

    if valid.empty:
        raise SystemExit(
            "No usable rows. Check the layout and the drop reasons above. If the "
            "outcome column holds codes rather than 0/1, name them with "
            "--event-values / --censor-values, or add them to EVENT_CODES at the "
            "top of this script."
        )

    results = [fit_group(f, g, args.alpha) for f, g in valid.groupby("force_pN")]
    summary = pd.DataFrame(results).sort_values("force_pN")

    if args.dry_run:
        print("Dry run: %d force group(s) would be fit, nothing written"
              % len(summary), file=sys.stderr)
        if not args.quiet:
            print(summary.to_string(index=False))
        return 0

    output = args.output or default_output(args.input)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_table(summary, output)
    print("Fit %d force group(s), wrote %s" % (len(summary), output),
          file=sys.stderr)

    if not args.quiet:
        print(summary.to_string(index=False))

    return 0


if __name__ == "__main__":
    sys.exit(main())
