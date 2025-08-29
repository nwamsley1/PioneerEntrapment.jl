"""
Protein-specific EFDR calculation functions.
"""

function add_protein_efdr_columns!(protein_results::DataFrame;
                                  method_types::Vector=[CombinedEFDR, PairedEFDR],
                                  score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:pg_score, :qval)],
                                  r::Float64=1.0,
                                  paired_stride::Int=5)
    if !hasproperty(protein_results, :entrap_id)
        error("DataFrame must have :entrap_id column")
    end
    entrap_labels = protein_results.entrap_id
    for (score_col, qval_col) in score_qval_pairs
        sort!(protein_results, [score_col, :entrap_id], rev = [true, false])
        if !hasproperty(protein_results, score_col)
            @warn "Column $score_col not found in DataFrame, skipping..."; continue
        end
        if !hasproperty(protein_results, qval_col)
            @warn "Column $qval_col not found in DataFrame, skipping..."; continue
        end
        original_target_col = Symbol(String(score_col) * "_original_target")
        if !hasproperty(protein_results, original_target_col)
            @warn "Column $original_target_col not found. Make sure to run add_original_target_protein_scores! first."; continue
        end
        scores = Float64.(protein_results[!, score_col])
        original_target_scores = Float64.(protein_results[!, original_target_col])
        qvals = Float64.(protein_results[!, qval_col])
        for method_type in method_types
            method_name = method_type == CombinedEFDR ? "combined" : method_type == PairedEFDR ? "paired" : lowercase(string(method_type))
            efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
            method = method_type(scores, original_target_scores, entrap_labels, qvals; r=r)
            efdr_values = method_type == PairedEFDR ? calculate_efdr(method; stride=paired_stride) : calculate_efdr(method)
            protein_results[!, efdr_col] = efdr_values
        end
    end
    return nothing
end

function compare_protein_efdr_methods(protein_results::DataFrame, qval_col::Symbol, score_col::Symbol)
    combined_efdr_col = Symbol(String(score_col) * "_combined_efdr")
    paired_efdr_col = Symbol(String(score_col) * "_paired_efdr")
    required_cols = [qval_col, combined_efdr_col, paired_efdr_col, :entrap_id]
    missing_cols = [col for col in required_cols if !hasproperty(protein_results, col)]
    if !isempty(missing_cols)
        error("Missing required columns: $missing_cols")
    end
    entrap_labels = protein_results.entrap_id
    thresholds = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]
    results = []
    for threshold in thresholds
        qval_mask = protein_results[!, qval_col] .<= threshold
        qval_n = sum(qval_mask)
        qval_actual_fdr = qval_n > 0 ? sum(qval_mask .& (entrap_labels .!= 0)) / qval_n : 0.0
        combined_mask = protein_results[!, combined_efdr_col] .<= threshold
        combined_n = sum(combined_mask)
        combined_actual_fdr = combined_n > 0 ? sum(combined_mask .& (entrap_labels .!= 0)) / combined_n : 0.0
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
    end
    return DataFrame(results)
end

function calculate_protein_efdr_calibration_error(protein_results::DataFrame, qval_col::Symbol, efdr_col::Symbol)
    required_cols = [qval_col, efdr_col, :entrap_id]
    missing_cols = [col for col in required_cols if !hasproperty(protein_results, col)]
    if !isempty(missing_cols)
        error("Missing required columns: $missing_cols")
    end
    entrap_labels = protein_results.entrap_id
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
