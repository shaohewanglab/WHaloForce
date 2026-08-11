"""Regenerate every figure panel from the committed source-data CSVs.

Discovers every src/<fig>/plot_*.py script, runs each in turn, and writes the
panels into figures/. Any new plot_*.py added later (for example the pending
stretch and C. elegans panels) is picked up automatically.

Also runs every src/<fig>/stats_*.py script -- panels whose p-values are
reported in the legend rather than drawn on the panel (Fig. 6b) compute them in
a separate stats script, and those need to be regenerated alongside the
figures.

Usage:
    python make_all_figures.py
"""

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent
SRC = REPO / "src"


def main():
    scripts = sorted(p for pattern in ("*/plot_*.py", "*/stats_*.py")
                     for p in SRC.glob(pattern))
    if not scripts:
        print("no plot or stats scripts found under src/")
        return 1
    failures = []
    for script in scripts:
        rel = script.relative_to(REPO)
        print(f"=== {rel} ===")
        result = subprocess.run([sys.executable, str(script)], cwd=REPO)
        if result.returncode != 0:
            failures.append(rel)
    print()
    print(f"ran {len(scripts)} scripts, {len(failures)} failed")
    for f in failures:
        print(f"  FAILED: {f}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
