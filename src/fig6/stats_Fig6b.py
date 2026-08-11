"""
Fig. 6b statistics: two-sided Mann-Whitney U tests with Holm-Bonferroni
correction on spermathecal local lifetime contrast (tau_Spt(norm)).

The published Fig. 6b panel carries no p-values -- its legend reports only
mean +/- SD and n -- but the Results text makes two claims about it:

  (1) "CeLam-WF, but not the corresponding cpHalo or C-terminal WHaloForce
      controls, exhibited higher fluorescence lifetimes at the spermatheca
      than in neighboring tissues" (peak ovulation, day-2 adults);
  (2) "The elevation was absent in L4 animals before ovulation began and
      became much less pronounced in aged adults after ovulation had ceased."

Claim (1) is tested here as CeLam-WF vs each control construct within each
age, and claim (2) as CeLam-WF at D2 vs CeLam-WF at every other age. The two
claims are separate analyses, so family-wise error is controlled within each
one separately (Holm across 12 construct comparisons; Holm across 5 age
comparisons) rather than across all 17 tests pooled.

The unit of analysis is the worm: each row of Fig6b.csv is one animal, and
tau_Spt(norm) is already a within-animal contrast (spermatheca minus the mean
of proximal germline and uterus), so the two constructs being compared at a
given age are independent samples and Mann-Whitney U applies.

Run:
    python stats_Fig6b.py

Writes Fig6b_stats.csv next to this file and prints the same table.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --- shared Fig. 6 statistics ---
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _stats import holm, mwu  # noqa: E402

# --- house style ---
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from style import DATA_DIR  # noqa: E402


VALUE_COL = "tau_Spt_norm_ns"
AGE_ORDER = ["L4", "D1", "D2", "D3", "D4", "D5"]
SENSOR = "CeLam-WF"
CONTROLS = ["CeLam-cpH", "CeLam-C-WF"]
REF_AGE = "D2"                  # peak ovulation, the reference for the age series


def describe(vals):
    return dict(n=vals.size, mean=vals.mean(), sd=vals.std(ddof=1),
                median=np.median(vals))


def main():
    df = pd.read_csv(DATA_DIR / "Fig6" / "Fig6b.csv")

    def vals(construct, age):
        m = (df["construct"] == construct) & (df["age_day"] == age)
        return df.loc[m, VALUE_COL].to_numpy()

    rows = []

    # --- Family 1: sensor vs each control construct, within each age ---
    for age in AGE_ORDER:
        a = vals(SENSOR, age)
        for control in CONTROLS:
            b = vals(control, age)
            U, p, eff = mwu(a, b)
            da, db = describe(a), describe(b)
            rows.append(dict(
                family="construct", comparison=f"{SENSOR} vs {control}",
                age=age, group_a=SENSOR, group_b=control,
                n_a=da["n"], n_b=db["n"],
                mean_a=da["mean"], mean_b=db["mean"],
                sd_a=da["sd"], sd_b=db["sd"],
                median_a=da["median"], median_b=db["median"],
                U=U, effect_size=eff, p_raw=p,
            ))

    # --- Family 2: sensor at peak ovulation vs sensor at every other age ---
    a = vals(SENSOR, REF_AGE)
    for age in [x for x in AGE_ORDER if x != REF_AGE]:
        b = vals(SENSOR, age)
        U, p, eff = mwu(a, b)
        da, db = describe(a), describe(b)
        rows.append(dict(
            family="age", comparison=f"{SENSOR} {REF_AGE} vs {age}",
            age=age, group_a=f"{SENSOR} {REF_AGE}", group_b=f"{SENSOR} {age}",
            n_a=da["n"], n_b=db["n"],
            mean_a=da["mean"], mean_b=db["mean"],
            sd_a=da["sd"], sd_b=db["sd"],
            median_a=da["median"], median_b=db["median"],
            U=U, effect_size=eff, p_raw=p,
        ))

    out = pd.DataFrame(rows)

    # Holm within each family, not across the two pooled.
    out["p_holm"] = np.nan
    out["n_in_family"] = 0
    for family, idx in out.groupby("family").groups.items():
        out.loc[idx, "p_holm"] = holm(out.loc[idx, "p_raw"])
        out.loc[idx, "n_in_family"] = len(idx)
    out["significant"] = out["p_holm"] < 0.05

    out_path = Path(__file__).resolve().parent / "Fig6b_stats.csv"
    out.to_csv(out_path, index=False)

    # Printed report: raw and adjusted p side by side, so the cost of the
    # correction is visible per comparison.
    pd.set_option("display.width", 200)
    for family, label in [("construct", f"{SENSOR} vs control constructs, by age"),
                          ("age", f"{SENSOR}, {REF_AGE} vs other ages")]:
        sub = out[out["family"] == family]
        m = len(sub)
        print(f"\n[Fig6b] {label}  (Mann-Whitney U, two-sided; "
              f"Holm-Bonferroni across m={m})")
        for _, r in sub.iterrows():
            star = "*" if r["significant"] else " "
            print(f"  {star} {r['comparison']:<32s} {r['age']:>3s}  "
                  f"n={int(r['n_a']):>2d}/{int(r['n_b']):<2d}  "
                  f"U={r['U']:>6.1f}  PS={r['effect_size']:.3f}  "
                  f"p={r['p_raw']:.3g}  p_Holm={r['p_holm']:.3g}")

    print(f"\nWrote {out_path.name}")


if __name__ == "__main__":
    main()
