# PioneerEntrapment.jl

Entrapment-based empirical FDR analysis for DIA proteomics (Pioneer/DIA-NN).

- Inputs: precursor-level Pioneer results + spectral library (for precursor EFDR), or protein-level outputs (Pioneer/DIA-NN) for protein EFDR.
- Outputs: Arrow tables with EFDR columns, Markdown reports, and plots (PNG/PDF).
- CLI: `pioneer-entrapment` supports `precursor`, `protein`, and `both` modes.

Quick start:

```bash
# Precursor-level (requires library)
pioneer-entrapment --mode precursor \
  --precursor-results path/to/precursors.arrow \
  --library path/to/library.arrow \
  --outdir results/precursor \
  --paired-step 5

# Protein-level
pioneer-entrapment --mode protein \
  --protein-results path/to/proteins.arrow \
  --outdir results/protein \
  --paired-step 5

# Both
pioneer-entrapment --mode both \
  --precursor-results path/to/precursors.arrow \
  --library path/to/library.arrow \
  --protein-results path/to/proteins.arrow \
  --outdir results \
  --paired-step 5
```

From Julia:

```julia
using PioneerEntrapment
run_efdr_analysis("precursors.arrow", "library.arrow"; output_dir="out")
run_protein_efdr_analysis("proteins.arrow"; output_dir="out")
run_both_analyses(; precursor_results_path="precursors.arrow", library_precursors_path="library.arrow", protein_results_path="proteins.arrow", output_dir="out")
```

Batch analysis:

```bash
# Scan a directory tree for Pioneer outputs and run EFDR per folder
bin/pioneer-entrapment-batch \
  --base-dir /Users/nathanwamsley/Data/MS_DATA/ARROW/MTAC3PAstralStandard/MTAC3P_entrapR1_feature-fix-mbr-dee81328-F_08-29-2025 \
  --library /Users/nathanwamsley/Data/SPEC_LIBS/entrapment_08-20-2025/altimeter_3P_len7o40_ch2o3_mc1_MTACAstral_Jul312025_entrapR1.poin \
  --out-name efdr_output \
  --paired-step 5

# Dry-run to preview commands without executing
bin/pioneer-entrapment-batch --base-dir <BASE> --library <LIB> --dry-run
```

License: MIT

## Example: Julia REPL on local data

This example runs EFDR from a Julia REPL using a .poin library directory (auto-detects `precursors_table.arrow`) and a results folder with `precursors_long.arrow`.

```bash
cd /Users/nathanwamsley/Projects/PioneerEntrapment.jl
julia --project=. -q
```

```julia
using Pkg; Pkg.instantiate()
using PioneerEntrapment

library = "/Users/nathanwamsley/Data/SPEC_LIBS/entrapment_08-20-2025/altimeter_yeast_len7o40_ch2o3_mc1_MTACAstral_Aug302025_entrapR1.poin/altimeter_yeast_len7o40_ch2o3_mc1_MTACAstral_Aug302025_entrapR1.poin.poin"
precursors = "/Users/nathanwamsley/Data/MS_DATA/MTAC_Y_Astral/YEAST_3MIN/MTAC_Y_entrapR1_feature-fix-mbr_D_08-30-2025/precursors_long.arrow"
# If TSV instead, use:
# precursors = "/Users/nathanwamsley/Data/MS_DATA/MTAC_Y_Astral/YEAST_3MIN/MTAC_Y_entrapR1_feature-fix-mbr_D_08-30-2025/precursors_long.tsv"
out = "/Users/nathanwamsley/Data/MS_DATA/MTAC_Y_Astral/YEAST_3MIN/MTAC_Y_entrapR1_feature-fix-mbr_D_08-30-2025/efdr_out"

run_efdr_analysis(precursors, library;
    output_dir=out,
    r_lib=1.0,
    paired_stride=10,      # sampling stride for paired EFDR
    plot_formats=[:png, :pdf],
    verbose=true,
)
```

Notes
- Passing the library as the `.poin.poin` directory is supported; the loader resolves `precursors_table.arrow` internally.
- Prefer `precursors_long.arrow` when available; `.tsv` also works.

## Using Revise for live reloads

To avoid restarting Julia on code changes:

```julia
using Pkg; Pkg.add("Revise")  # one-time
using Revise
using PioneerEntrapment

# Run an example script (re-runs on each includet call)
Revise.includet("scripts/replicate_plot_example.txt")

# Edit files under src/ and save; then re-run
Revise.includet("scripts/replicate_plot_example.txt")
```

## Multi-run comparison (replicates)

Compare EFDR curves from multiple runs on a single plot. Colors indicate method (blue=Combined, red=Paired); each replicate is an additional solid line.

```julia
using PioneerEntrapment

library = "/path/to/library/.poin/.poin"  # resolves precursors_table.arrow inside
rep1 = "/path/to/run1/precursors_long.arrow"
rep2 = "/path/to/run2/precursors_long.arrow"

replicates = [
  (precursor_results_path=rep1, library_precursors_path=library, label="run1"),
  (precursor_results_path=rep2, library_precursors_path=library, label="run2"),
]

run_efdr_replicate_plots(replicates;
  output_dir="./efdr_compare",
  r_lib=1.0,
  paired_stride=10,
  plot_formats=[:png, :pdf],
  verbose=true,
)
```

Outputs are written to the `output_dir` for each score (e.g., `efdr_comparison_replicates_global_prob.png`).

CLI via config

- JSON:
  - bin/pioneer-entrapment --mode replicates \
    --replicates-config scripts/replicates_example.json \
    --outdir ./efdr_compare --paired-step 10
- YAML (requires YAML.jl):
  - bin/pioneer-entrapment --mode replicates \
    --replicates-config scripts/replicates_example.yaml \
    --outdir ./efdr_compare --paired-step 10

JSON expects an array of objects. YAML expects a top-level key `replicates:` with a list. Each item needs:
- `precursor_results_path`
- `library_precursors_path`
- `label` (optional)
