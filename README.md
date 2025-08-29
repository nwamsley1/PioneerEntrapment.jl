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

License: MIT
