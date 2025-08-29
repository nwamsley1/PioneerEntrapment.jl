# PioneerEntrapment.jl

Entrapment-based empirical FDR analysis for DIA proteomics.

## Installation

```julia
] add https://github.com/your-org-or-user/PioneerEntrapment.jl
```

## Quick Start

```julia
using PioneerEntrapment
run_efdr_analysis("precursors.arrow", "library.arrow"; output_dir="out")
run_protein_efdr_analysis("proteins.arrow"; output_dir="out")
run_both_analyses(; precursor_results_path="precursors.arrow", library_precursors_path="library.arrow", protein_results_path="proteins.arrow", output_dir="out")
```

## CLI

```bash
pioneer-entrapment --mode precursor --precursor-results prec.arrow --library lib.arrow --outdir out
pioneer-entrapment --mode protein --protein-results protein.arrow --outdir out
pioneer-entrapment --mode both --precursor-results prec.arrow --library lib.arrow --protein-results protein.arrow --outdir out
```

## API Reference

- `run_efdr_analysis`
- `run_protein_efdr_analysis`
- `run_both_analyses`
