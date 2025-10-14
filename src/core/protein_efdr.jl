"""
Protein-level EFDR helpers.

These utilities mirror precursor-level EFDR computation but operate on protein
aggregates, respecting optional `:species` partitions and file identity.
"""

"""
    add_protein_efdr_columns!(protein_results; method_types=[CombinedEFDR,PairedEFDR], score_qval_pairs=[(:pg_score,:qval)], r=1.0, paired_stride=5, use_fast_paired=true, entrap_labels_override=nothing) -> Nothing

Compute EFDR for protein-level results and attach columns named like
`"<score>_combined_efdr"` and `"<score>_paired_efdr"`.

If `<score>_original_target` is missing, PairedEFDR for that score is skipped.
When `entrap_labels_override == nothing`, `:entrap_id` is used to label target
vs entrapment proteins (0 vs non-zero).

Parameters
- use_fast_paired: If true (default), use fast O(n log n) implementation for PairedEFDR. If false, use standard O(n²) implementation.
"""
function add_protein_efdr_columns!(protein_results::DataFrame;
                                  method_types::Vector=[CombinedEFDR, PairedEFDR],
                                  score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:pg_score, :qval)],
                                  r::Float64=1.0,
                                  paired_stride::Int=5,
                                  use_fast_paired::Bool=true,
                                  entrap_labels_override::Union{Nothing,AbstractVector}=nothing)
    entrap_labels = if entrap_labels_override === nothing
        if !hasproperty(protein_results, :entrap_id)
            error("DataFrame must have :entrap_id column (or provide entrap_labels_override)")
        end
        protein_results.entrap_id
    else
        collect(Int, entrap_labels_override)
    end
    for (score_col, qval_col) in score_qval_pairs
        if hasproperty(protein_results, :entrap_id)
            sort!(protein_results, [score_col, :entrap_id], rev = [true, false])
        else
            sort!(protein_results, [score_col], rev = [true])
        end
        if !hasproperty(protein_results, score_col)
            @warn "Column $score_col not found in DataFrame, skipping..."; continue
        end
        if !hasproperty(protein_results, qval_col)
            @warn "Column $qval_col not found in DataFrame, skipping..."; continue
        end
        original_target_col = Symbol(String(score_col) * "_original_target")
        scores = Float64.(protein_results[!, score_col])
        have_original = hasproperty(protein_results, original_target_col)
        original_target_scores = have_original ? Float64.(protein_results[!, original_target_col]) : copy(scores)
        qvals = Float64.(protein_results[!, qval_col])
        for method_type in method_types
            method_name = method_type == CombinedEFDR ? "combined" : method_type == PairedEFDR ? "paired" : lowercase(string(method_type))
            efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
            if method_type == PairedEFDR && !have_original
                @warn "Missing $original_target_col; skipping PairedEFDR for $score_col"
                continue
            end
            method = method_type(scores, original_target_scores, entrap_labels, qvals; r=r)
            efdr_values = if method_type == PairedEFDR
                if use_fast_paired
                    calculate_efdr_fast(method; cuts_mode=:all)
                else
                    calculate_efdr(method; stride=paired_stride)
                end
            else
                calculate_efdr(method)
            end
            protein_results[!, efdr_col] = efdr_values
        end
    end
    return nothing
end

"""
    compare_protein_efdr_methods(protein_results, qval_col, score_col; entrap_labels_override=nothing, include_paired=true) -> DataFrame

Summarize method calibration by comparing selections at a set of EFDR/q-value
thresholds against the actual entrapment fraction.

Returns one row per threshold with counts and realized FDR for q-value-based
selection, combined EFDR, and optionally paired EFDR.
"""
function compare_protein_efdr_methods(protein_results::DataFrame, qval_col::Symbol, score_col::Symbol; entrap_labels_override::Union{Nothing,AbstractVector}=nothing, include_paired::Bool=true)
    combined_efdr_col = Symbol(String(score_col) * "_combined_efdr")
    paired_efdr_col = Symbol(String(score_col) * "_paired_efdr")
    required_cols = [qval_col, combined_efdr_col]
    missing_cols = [col for col in required_cols if !hasproperty(protein_results, col)]
    if !isempty(missing_cols)
        error("Missing required columns: $missing_cols")
    end
    entrap_labels = entrap_labels_override === nothing ? (hasproperty(protein_results, :entrap_id) ? protein_results.entrap_id : error(":entrap_id missing; provide entrap_labels_override")) : collect(Int, entrap_labels_override)
    thresholds = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]
    results = []
    for threshold in thresholds
        qval_mask = protein_results[!, qval_col] .<= threshold
        qval_n = sum(qval_mask)
        qval_actual_fdr = qval_n > 0 ? sum(qval_mask .& (entrap_labels .!= 0)) / qval_n : 0.0
        combined_mask = protein_results[!, combined_efdr_col] .<= threshold
        combined_n = sum(combined_mask)
        combined_actual_fdr = combined_n > 0 ? sum(combined_mask .& (entrap_labels .!= 0)) / combined_n : 0.0
        if include_paired && hasproperty(protein_results, paired_efdr_col)
            paired_mask = protein_results[!, paired_efdr_col] .<= threshold
            paired_n = sum(paired_mask)
            paired_actual_fdr = paired_n > 0 ? sum(paired_mask .& (entrap_labels .!= 0)) / paired_n : 0.0
            push!(results, (
                threshold = threshold,
                qval_n = qval_n,
                qval_actual_fdr = qval_actual_fdr,
                combined_n = combined_n,
                combined_efdr = threshold,
                combined_actual_fdr = combined_actual_fdr,
                paired_n = paired_n,
                paired_efdr = threshold,
                paired_actual_fdr = paired_actual_fdr
            ))
        else
            push!(results, (
                threshold = threshold,
                qval_n = qval_n,
                qval_actual_fdr = qval_actual_fdr,
                combined_n = combined_n,
                combined_efdr = threshold,
                combined_actual_fdr = combined_actual_fdr
            ))
        end
    end
    return DataFrame(results)
end

"""
    calculate_protein_efdr_calibration_error(protein_results, qval_col, efdr_col; entrap_labels_override=nothing) -> (DataFrame, Float64)

Compute realized FDR across thresholds 0.001:0.001:0.1 for a given EFDR column
and return a table with the absolute error vs the nominal threshold, along with
the mean absolute error.
"""
function calculate_protein_efdr_calibration_error(protein_results::DataFrame, qval_col::Symbol, efdr_col::Symbol; entrap_labels_override::Union{Nothing,AbstractVector}=nothing)
    required_cols = [qval_col, efdr_col]
    missing_cols = [col for col in required_cols if !hasproperty(protein_results, col)]
    if !isempty(missing_cols)
        error("Missing required columns: $missing_cols")
    end
    entrap_labels = entrap_labels_override === nothing ? (hasproperty(protein_results, :entrap_id) ? protein_results.entrap_id : error(":entrap_id missing; provide entrap_labels_override")) : collect(Int, entrap_labels_override)
    thresholds = 0.001:0.001:0.1
    calibration_data = []
    for threshold in thresholds
        mask = protein_results[!, efdr_col] .<= threshold
        n_selected = sum(mask)
        if n_selected > 0
            actual_fdr = sum(mask .& (entrap_labels .!= 0)) / n_selected
            error = abs(actual_fdr - threshold)
            push!(calibration_data, (
                threshold = threshold,
                estimated_fdr = threshold,
                actual_fdr = actual_fdr,
                n_selected = n_selected,
                error = error
            ))
        end
    end
    cal_df = DataFrame(calibration_data)
    mean_error = isempty(cal_df) ? 0.0 : mean(cal_df.error)
    return cal_df, mean_error
end
