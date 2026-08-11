"""Smoke test for the reproducibility package.

Two groups of checks:

* Downstream (plotting): run every src/<fig>/plot_*.py and confirm it writes the
  expected figures/<Panel>.svg. This validates that each committed source-data
  CSV regenerates its panel.
* Upstream (the analysis steps that produced those CSVs): run each
  examples/<step>/check_*.py. Each has a small committed input and its expected
  output; see the README in each folder. A check whose optional dependencies are
  absent reports SKIP, not a failure.

The force-trace extraction for Fig. 1f-j and S2c/S2d has no example: those traces
were extracted from LUMICKS C-Trap HDF5 upstream of this package, and the
committed CSVs are the extracted traces themselves (see
docs/NON_SCRIPTABLE_STEPS.md).

Usage:
    python examples/run_examples.py
"""

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "src"
FIG = REPO / "figures"
EXAMPLES = Path(__file__).resolve().parent


def panel_stem(script: Path) -> str:
    # src/fig1/plot_fig1d.py -> Fig1d
    name = script.stem.replace("plot_", "")
    return name[0].upper() + name[1:]


def run_panels() -> list[Path]:
    """Regenerate every panel; return the list of failures."""
    scripts = sorted(p for p in SRC.glob("*/plot_*.py"))
    if not scripts:
        print("no plot scripts found under src/")
        return [Path("<no plot scripts>")]
    failures = []
    for script in scripts:
        rel = script.relative_to(REPO)
        r = subprocess.run([sys.executable, str(script)], cwd=REPO,
                           capture_output=True, text=True)
        stem = panel_stem(script)
        produced = (FIG / f"{stem}.svg").exists()
        ok = r.returncode == 0 and produced
        print(f"[{'OK ' if ok else 'FAIL'}] {rel} -> figures/{stem}.svg")
        if not ok:
            failures.append(rel)
            if r.returncode != 0:
                print(r.stderr.strip()[-500:])
    print()
    print(f"{len(scripts) - len(failures)}/{len(scripts)} panels regenerated")
    return failures


def run_step_checks() -> list[Path]:
    """Run each upstream step check; return the list of failures."""
    checks = sorted(EXAMPLES.glob("*/check_*.py"))
    if not checks:
        print("no step checks found under examples/")
        return []
    failures = []
    for check in checks:
        rel = check.relative_to(REPO)
        print(f"--- {rel}")
        r = subprocess.run([sys.executable, str(check)], cwd=REPO,
                           capture_output=True, text=True)
        out = r.stdout.strip()
        skipped = "[SKIP]" in out
        ok = r.returncode == 0
        for line in out.splitlines():
            print(f"    {line}")
        if not ok:
            print(r.stderr.strip()[-500:])
            failures.append(rel)
        status = "SKIP" if (ok and skipped) else ("OK " if ok else "FAIL")
        print(f"[{status}] {rel}")
        print()
    print(f"{len(checks) - len(failures)}/{len(checks)} step checks passed")
    return failures


def main():
    print("=" * 72)
    print("Upstream analysis steps (examples/<step>/check_*.py)")
    print("=" * 72)
    step_failures = run_step_checks()
    print()
    print("=" * 72)
    print("Figure panels from the committed source data (src/<fig>/plot_*.py)")
    print("=" * 72)
    panel_failures = run_panels()

    print()
    if step_failures or panel_failures:
        if step_failures:
            print(f"FAILED step checks: {[str(p) for p in step_failures]}")
        if panel_failures:
            print(f"FAILED panels: {[str(p) for p in panel_failures]}")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
