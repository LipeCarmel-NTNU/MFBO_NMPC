"""Parse the phi(z) surrogate vintage records into a flat CSV.

Reads results/<case>/registry/vintages/vintage_*.json and writes one row per
vintage to <storage>/surrogate_vintages.csv. This is the Python side of the
preprocessing: initial_preprocessing.m loads the CSV as the `surrogate`
variable and joins it to `evals` by (case, phi_vintage).

Each vintage governs a block of BO iterations (governs_iterations), so the fit
wall time is a property of the vintage, not of a single evaluation.

Usage:
    python parse_registry.py                 # cases from cases.txt
    python parse_registry.py --cases results_case1 results_case2   # override
"""

import argparse
import csv
import glob
import json
import os


def read_cases(path):
    """Read case subfolder names from cases.txt (one per line, # comments)."""
    cases = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                cases.append(line)
    return cases


COLUMNS = [
    "case", "vintage", "created_at", "gov_iter_start", "gov_iter_end",
    "n_doe_rows", "n_opt_rows_used", "n_runs_used", "n_samples",
    "SSE_a", "SSE_b", "SSE_lambda", "SSE_loss",
    "SSdU_a", "SSdU_b", "SSdU_lambda", "SSdU_loss",
    "fit_wall_s", "SSE_fit_wall_s", "SSdU_fit_wall_s",
]


def _governs(ctx):
    gi = ctx.get("governs_iterations") or [None, None]
    start = gi[0] if len(gi) > 0 else None
    end = gi[1] if len(gi) > 1 else None
    return start, end


def parse_vintage(case, path):
    with open(path, "r", encoding="utf-8") as fh:
        d = json.load(fh)
    ctx = d.get("context", {})
    tgt = d.get("targets", {})
    sse = tgt.get("SSE", {})
    ssdu = tgt.get("SSdU", {})
    wall = ctx.get("wall_s", {}).get("per_target", {})
    start, end = _governs(ctx)
    return {
        "case": case,
        "vintage": d.get("vintage"),
        "created_at": d.get("created_at"),
        "gov_iter_start": start,
        "gov_iter_end": end,
        "n_doe_rows": ctx.get("n_doe_rows"),
        "n_opt_rows_used": ctx.get("n_opt_rows_used"),
        "n_runs_used": sse.get("n_runs_used"),
        "n_samples": sse.get("n_samples"),
        "SSE_a": sse.get("a"),
        "SSE_b": sse.get("b"),
        "SSE_lambda": sse.get("lam"),
        "SSE_loss": sse.get("loss"),
        "SSdU_a": ssdu.get("a"),
        "SSdU_b": ssdu.get("b"),
        "SSdU_lambda": ssdu.get("lam"),
        "SSdU_loss": ssdu.get("loss"),
        "fit_wall_s": ctx.get("fit_wall_s"),
        "SSE_fit_wall_s": wall.get("SSE", {}).get("total"),
        "SSdU_fit_wall_s": wall.get("SSdU", {}).get("total"),
    }


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(here)              # Result analysis -> repo root
    default_storage = os.path.join(here, "storage")

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cases", nargs="+", default=None,
                    help="Override cases.txt with an explicit list.")
    ap.add_argument("--cases-file", default=os.path.join(here, "cases.txt"))
    ap.add_argument("--results-root", default=os.path.join(repo_root, "results"))
    ap.add_argument("--storage", default=default_storage)
    args = ap.parse_args()

    cases = args.cases if args.cases else read_cases(args.cases_file)

    os.makedirs(args.storage, exist_ok=True)
    rows = []
    for case in cases:
        vdir = os.path.join(args.results_root, case, "registry", "vintages")
        files = sorted(glob.glob(os.path.join(vdir, "vintage_*.json")))
        if not files:
            print(f"  [warn] no vintage files under {vdir}")
            continue
        for f in files:
            rows.append(parse_vintage(case, f))
        print(f"  {case}: {len(files)} vintages")

    out_csv = os.path.join(args.storage, "surrogate_vintages.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} rows -> {out_csv}")


if __name__ == "__main__":
    main()
