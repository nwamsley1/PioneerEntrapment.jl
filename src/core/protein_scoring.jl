using DataFrames
using Dictionaries

function add_original_target_protein_scores!(protein_results::DataFrame; score_col=:pg_score)
    required_cols = [:protein, :entrap_id]
    missing_cols = [col for col in required_cols if !hasproperty(protein_results, col)]
    if !isempty(missing_cols)
        error("DataFrame missing required columns: $missing_cols")
    end
    if !hasproperty(protein_results, score_col)
        error("DataFrame must have :$score_col column.")
    end
    file_col = if hasproperty(protein_results, :ms_file_idx)
        :ms_file_idx
    elseif hasproperty(protein_results, :file_name)
        :file_name
    else
        error("DataFrame must have either :ms_file_idx or :file_name column")
    end
    original_target_col = Symbol(String(score_col) * "_original_target")
    # Keys: (file, species?, protein). If :species missing, use empty string.
    protein_to_target = Dictionary{Tuple{Any, String, String}, Float32}()
    for row in eachrow(protein_results)
        if row.entrap_id == 0 && !ismissing(row[score_col])
            species = hasproperty(protein_results, :species) ? String(row.species) : ""
            key = (row[file_col], species, row.protein)
            if haskey(protein_to_target, key)
                error("Duplicate target protein found: protein '$(row.protein)' appears multiple times with entrap_id=0 in file '$(row[file_col])'.")
            end
            insert!(protein_to_target, key, Float32(row[score_col]))
        end
    end
    original_target_scores = Float32[]
    for row in eachrow(protein_results)
        if !ismissing(row[score_col])
            if row.entrap_id == 0
                push!(original_target_scores, Float32(row[score_col]))
            else
                species = hasproperty(protein_results, :species) ? String(row.species) : ""
                key = (row[file_col], species, row.protein)
                if haskey(protein_to_target, key)
                    push!(original_target_scores, protein_to_target[key])
                else
                    push!(original_target_scores, -1.0f0)
                end
            end
        else
            push!(original_target_scores, -1.0f0)
        end
    end
    protein_results[!, original_target_col] = original_target_scores
    return nothing
end

function add_original_target_protein_scores!(protein_results::DataFrame, score_cols::Vector{Symbol})
    for score_col in score_cols
        add_original_target_protein_scores!(protein_results; score_col=score_col)
    end
    return nothing
end

function create_global_protein_results_df(protein_results::DataFrame; score_col::Symbol=:global_pg_score)
    required_cols = [:protein, score_col]
    missing_cols = [col for col in required_cols if !hasproperty(protein_results, col)]
    if !isempty(missing_cols)
        error("DataFrame missing required columns: $missing_cols")
    end
    file_col = if hasproperty(protein_results, :ms_file_idx)
        :ms_file_idx
    elseif hasproperty(protein_results, :file_name)
        :file_name
    else
        # Default to file_name if present, else synthesize
        file_col = :file_name
    end
    protein_results_copy = copy(protein_results)
    grouped = hasproperty(protein_results_copy, :species) ? groupby(protein_results_copy, [:species, :protein]) : groupby(protein_results_copy, [:protein])
    global_df = combine(grouped) do group
        valid_rows = group[.!ismissing.(group[!, score_col]), :]
        if nrow(valid_rows) == 0
            return similar(group, 0)
        end
        best_idx = argmax(valid_rows[!, score_col])
        return valid_rows[best_idx:best_idx, :]
    end
    if nrow(global_df) > 0
        if hasproperty(protein_results, :ms_file_idx)
            global_df[!, :ms_file_idx] = fill(0, nrow(global_df))
        elseif hasproperty(protein_results, :file_name)
            global_df[!, :file_name] = fill("global", nrow(global_df))
        else
            global_df[!, :file_name] = fill("global", nrow(global_df))
        end
    end
    return global_df
end
