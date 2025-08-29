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
    if file_col == :ms_file_idx
        protein_to_target = Dictionary{Tuple{Int, String, String}, Float32}()
    else
        protein_to_target = Dictionary{Tuple{String, String, String}, Float32}()
    end
    for row in eachrow(protein_results)
        if row.entrap_id == 0 && !ismissing(row[score_col])
            key = (row[file_col], row.species, row.protein)
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
                key = (row[file_col], row.species, row.protein)
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
    required_cols = [:species,:protein,:file_name,score_col]
    missing_cols = [col for col in required_cols if !hasproperty(protein_results, col)]
    if !isempty(missing_cols)
        error("DataFrame missing required columns: $missing_cols")
    end
    file_col = if hasproperty(protein_results, :ms_file_idx)
        :ms_file_idx
    elseif hasproperty(protein_results, :file_name)
        :file_name
    else
        error("DataFrame must have either :ms_file_idx or :file_name column")
    end
    protein_results_copy = copy(protein_results)
    grouped = groupby(protein_results_copy, [:species,:protein])
    global_df = combine(grouped) do group
        valid_rows = group[.!ismissing.(group[!, score_col]), :]
        if nrow(valid_rows) == 0
            return similar(group, 0)
        end
        best_idx = argmax(valid_rows[!, score_col])
        return valid_rows[best_idx:best_idx, :]
    end
    if nrow(global_df) > 0
        if file_col == :ms_file_idx
            global_df[!, file_col] .= 0
        else
            global_df[!, file_col] .= "global"
        end
    end
    return global_df
end

