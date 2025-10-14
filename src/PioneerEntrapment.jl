module PioneerEntrapment

using DataFrames
using Arrow
using Tables
using Printf
using Dates
using CSV
using Plots
using Statistics
using ProgressBars
using Dictionaries
using Markdown

# Core, analysis, plotting, and API
include("core/efdr_methods.jl")
include("core/paired_fast.jl")  # Must come after efdr_methods.jl (defines PairedEFDR) but before protein_efdr.jl
include("core/entrapment_pairing.jl")
include("core/protein_entrapment_pairing.jl")
include("core/scoring.jl")
include("core/protein_scoring.jl")
include("core/protein_efdr.jl")
include("analysis/efdr_analysis.jl")
include("analysis/calibration.jl")
include("plotting/efdr_plots.jl")
include("api.jl")
include("cli.jl")

# Re-exports (public API)
export run_efdr_analysis, run_protein_efdr_analysis, run_both_analyses, run_efdr_plots, run_efdr_replicate_plots

# Types
export EFDRMethod, CombinedEFDR, PairedEFDR

# Core functions
export calculate_efdr, add_efdr_columns!
export add_entrap_pair_ids!
export assign_protein_entrapment_pairs!, add_protein_entrap_pair_ids!
export add_original_target_scores!, get_complement_score
export add_original_target_protein_scores!, create_global_protein_results_df
export add_protein_efdr_columns!, compare_protein_efdr_methods, calculate_protein_efdr_calibration_error
export getModKey

# Analysis
export compare_efdr_methods, calculate_efdr_calibration_error
export analyze_efdr_at_threshold, print_efdr_comparison_table

# Plotting
export plot_efdr_comparison, plot_efdr_vs_qval
export plot_multiple_efdr_comparisons, save_efdr_plots, plot_efdr_comparison_replicates, save_efdr_replicate_plots

end # module
