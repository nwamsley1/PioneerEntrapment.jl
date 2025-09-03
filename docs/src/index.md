# PioneerEntrapment.jl

Entrapment-based empirical FDR analysis for DIA proteomics.

## Installation

```julia
] add https://github.com/nwamsley1/PioneerEntrapment.jl
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

```@autodocs
Modules = [PioneerEntrapment]
Order = [:function, :type]
```

## Example: Julia REPL on local data

Run EFDR from a Julia REPL using a .poin library directory (auto-detects `precursors_table.arrow`) and a results folder with `precursors_long.arrow`.

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

## Live Reload with Revise

```julia
using Pkg; Pkg.add("Revise")
using Revise
using PioneerEntrapment

Revise.includet("scripts/replicate_plot_example.txt")

# after edits in src/, just call again
Revise.includet("scripts/replicate_plot_example.txt")
```

## Multi-run comparison (replicates)

Plot EFDR from multiple runs in a single figure (one color per method, solid lines per replicate).

```julia
using PioneerEntrapment

library = "/path/to/library/.poin/.poin"
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

See the output directory for files named like `efdr_comparison_replicates_global_prob.png`.
