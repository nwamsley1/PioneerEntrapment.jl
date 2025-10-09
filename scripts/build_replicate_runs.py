#!/usr/bin/env python3
"""
Generate a Julia commands file (replicate_runs.txt) with back-to-back
run_efdr_replicate_plots(...) calls, grouped by condition from a TSV summary.

Expected TSV headers (case-insensitive):
  - library_json: path to a small JSON file that contains a key like
    library_precursors_path / library / lib / path
  - dataset
  - rep_dir: directory containing precursors_long.arrow
  - q_value
  - first_pass
  - mbr

Usage:
  python3 scripts/build_replicate_runs.py \
    --base /Users/nathanwamsley/Data/MS_DATA/EntrapmentTests \
    --summary run_summary_09-03-2025.tsv

This writes: <base>/replicate_runs.txt
"""

import argparse
import csv
import json
import os
from collections import defaultdict


REQUIRED_COLS = [
    "library_json",
    "dataset",
    "rep_dir",
    "q_value",
    "first_pass",
    "mbr",
]


def read_tsv(path):
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise SystemExit(f"No rows found in TSV: {path}")
    # normalize header keys to lowercase
    norm_rows = []
    for r in rows:
        norm_rows.append({k.lower(): v for k, v in r.items()})
    headers = set(norm_rows[0].keys())
    missing = [c for c in REQUIRED_COLS if c not in headers]
    if missing:
        raise SystemExit(f"Missing required headers in TSV: {missing}")
    return norm_rows


def resolve_path(base, p):
    if os.path.isabs(p):
        return p
    return os.path.normpath(os.path.join(base, p))


def read_library_path(json_file):
    with open(json_file, "r", encoding="utf-8") as f:
        data = json.load(f)
    for key in ("library_precursors_path", "library", "lib", "path"):
        if key in data and isinstance(data[key], str) and data[key]:
            return data[key]
    raise SystemExit(f"Could not find library path in JSON: {json_file}")


def jstr(s: str) -> str:
    """Escape for Julia string literal."""
    return s.replace("\\", "\\\\").replace("\"", "\\\"")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="EntrapmentTests base folder")
    ap.add_argument("--summary", required=True, help="Run summary TSV (path or name under base)")
    args = ap.parse_args()

    base = os.path.abspath(args.base)
    summary = args.summary
    if not os.path.isabs(summary):
        summary = os.path.join(base, summary)
    if not os.path.isfile(summary):
        raise SystemExit(f"Run summary not found: {summary}")

    rows = read_tsv(summary)

    # Group rows by condition
    groups = defaultdict(list)
    for r in rows:
        key = (
            resolve_path(base, r["library_json"]),
            r["dataset"],
            str(r["q_value"]),
            str(r["first_pass"]),
            str(r["mbr"]),
        )
        groups[key].append(resolve_path(base, r["rep_dir"]))

    out_root = os.path.join(base, "efdr_compare")
    os.makedirs(out_root, exist_ok=True)
    out_file = os.path.join(base, "replicate_runs.txt")

    lines = ["# Paste these into the Julia REPL (package activated)\n"]

    for (lib_json, dataset, qv, fp, mbr), rep_dirs in sorted(groups.items()):
        lib_tag = os.path.basename(lib_json).replace(".json", "")
        try:
            lib_path = read_library_path(lib_json)
        except SystemExit as e:
            print(f"Skipping {dataset}, q={qv}, fp={fp}, mbr={mbr}: {e}")
            continue

        # Build replicate named tuples
        entries = []
        for rd in rep_dirs:
            prec = os.path.join(rd, "precursors_long.arrow")
            if not os.path.isfile(prec):
                print(f"Warning: missing {prec}; skipping rep")
                continue
            entries.append(
                f"(precursor_results_path=\"{jstr(prec)}\", "
                f"library_precursors_path=\"{jstr(lib_path)}\", "
                f"label=\"{jstr(os.path.basename(rd))}\")"
            )
        if not entries:
            continue

        outdir = os.path.join(out_root, f"{dataset}_q{qv}_fp{fp}_mbr{mbr}_{lib_tag}")
        cmd = (
            f"run_efdr_replicate_plots([" + ", ".join(entries) + "];\n"
            f"  output_dir=\"{jstr(outdir)}\",\n"
            f"  r_lib=1.0,\n"
            f"  paired_stride=10,\n"
            f"  plot_formats=[:png, :pdf],\n"
            f"  verbose=true,\n"
            f")\n"
        )
        lines.append(cmd + "\n")

    with open(out_file, "w", encoding="utf-8") as f:
        f.writelines(lines)

    print(f"Wrote: {out_file}")


if __name__ == "__main__":
    main()

