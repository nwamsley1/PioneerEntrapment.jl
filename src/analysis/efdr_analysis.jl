using DataFrames
using Statistics

function analyze_efdr_at_threshold(df::DataFrame, qval_col::Symbol, efdr_col::Symbol, threshold::Float64, library_precursors::DataFrame)
    sorted_indices = sortperm(df[!, qval_col])
    sorted_df = df[sorted_indices, :]
    passing_mask = sorted_df[!, qval_col] .<= threshold
    passing_df = sorted_df[passing_mask, :]
    if nrow(passing_df) == 0
        return (threshold=threshold, n_passing=0, empirical_fdr=0.0, n_targets=0, n_entrapments=0)
    end
    entrap_labels = [library_precursors.entrapment_group_id[pid] for pid in passing_df.precursor_idx]
    n_targets = sum(entrap_labels .== 0)
    n_entrapments = sum(entrap_labels .> 0)
    empirical_fdr = passing_df[end, efdr_col]
    return (threshold = threshold, n_passing = nrow(passing_df), empirical_fdr = empirical_fdr, n_targets = n_targets, n_entrapments = n_entrapments)
end

function compare_efdr_methods(df::DataFrame, qval_col::Symbol, score_col::Symbol, library_precursors::DataFrame; thresholds::Vector{Float64}=[0.001, 0.01, 0.05, 0.1])
    combined_col = Symbol(String(score_col) * "_combined_efdr")
    paired_col = Symbol(String(score_col) * "_paired_efdr")
    results = DataFrame()
    for threshold in thresholds
        qval_sorted = sortperm(df[!, qval_col])
        first_failing = findfirst(df[qval_sorted, qval_col] .> threshold)
        qval_passing = isnothing(first_failing) ? df[qval_sorted, :] : df[qval_sorted[1:first_failing-1], :]
        qval_entrap_labels = [library_precursors.entrapment_group_id[pid] for pid in qval_passing.precursor_idx]
        combined_result = analyze_efdr_at_threshold(df, qval_col, combined_col, threshold, library_precursors)
        paired_result = analyze_efdr_at_threshold(df, qval_col, paired_col, threshold, library_precursors)
        push!(results, (
            threshold = threshold,
            qval_n = nrow(qval_passing),
            qval_actual_fdr = sum(qval_entrap_labels .> 0) / max(1, length(qval_entrap_labels)),
            combined_n = combined_result.n_passing,
            combined_efdr = combined_result.empirical_fdr,
            paired_n = paired_result.n_passing,
            paired_efdr = paired_result.empirical_fdr
        ))
    end
    return results
end

function print_efdr_comparison_table(comparison_df::DataFrame)
    println("\nEFDR Method Comparison")
    println("="^80)
    println("Threshold | Q-val IDs | Actual FDR | Combined IDs | Combined EFDR | Paired IDs | Paired EFDR")
    println("-"^80)
    for row in eachrow(comparison_df)
        @printf("%9.3f | %9d | %10.4f | %12d | %13.4f | %10d | %11.4f\n",
                row.threshold, row.qval_n, row.qval_actual_fdr,
                row.combined_n, row.combined_efdr,
                row.paired_n, row.paired_efdr)
    end
    println("="^80)
end

export analyze_efdr_at_threshold, compare_efdr_methods, print_efdr_comparison_table

