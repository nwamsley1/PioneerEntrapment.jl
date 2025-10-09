using DataFrames
using Dictionaries

"""
    add_original_target_protein_scores!(protein_results::DataFrame; score_col=:pg_score) -> Nothing

Add a `"<score>_original_target"` column capturing, for each protein row,
the score of the original target protein (entrap_id == 0) with the same
`file` and `species`. For protein group identifiers that contain multiple
entries separated by `;`, the first entry (i.e., `split(protein,';')[1]`) is
used as the group key so targets/entrapments match on the same canonical ID.

If the row itself is the original target, its own score is used. If there is no
matching original target, the value is set to `-1.0f0`.

Required columns
- protein_results: `:protein`, `:entrap_id`, `:species`, `score_col`, and either `:ms_file_idx` or `:file_name`.
"""
function add_original_target_protein_scores!(protein_results::DataFrame; score_col=:pg_score)
    required_cols = [:protein, :entrap_id, :species]
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
    # Keys: (file, species, protein). Species column is always required.
    protein_to_target = Dictionary{Tuple{Any, String, String}, Float32}()
    for row in eachrow(protein_results)
        if row.entrap_id == 0 && !ismissing(row[score_col])
            species = String(row.species)
            protkey = first(split(String(row.protein), ';'))
            key = (row[file_col], species, protkey)
            if haskey(protein_to_target, key)
                error("Duplicate target protein found: protein key '$(protkey)' appears multiple times with entrap_id=0 in file '$(row[file_col])'.")
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
                species = String(row.species)
                protkey = first(split(String(row.protein), ';'))
                key = (row[file_col], species, protkey)
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

"""
    add_original_target_protein_scores!(protein_results::DataFrame, score_cols::Vector{Symbol}) -> Nothing

Vectorized overload to compute original target columns for multiple score fields.
"""
function add_original_target_protein_scores!(protein_results::DataFrame, score_cols::Vector{Symbol})
    for score_col in score_cols
        add_original_target_protein_scores!(protein_results; score_col=score_col)
    end
    return nothing
end

"""
    create_global_protein_results_df(protein_results::DataFrame; score_col::Symbol=:global_pg_score) -> DataFrame

Select the best row across files per grouping of `(species?, protein, entrap_id?)`.

- Groups by `:species` if present, always by `:protein`, and by `:entrap_id` if present.
- For each group, picks the row with the maximum `score_col`.
- If `:ms_file_idx` exists, sets it to 0 for returned rows; otherwise assigns `:file_name => "global"`.

Notes
- Decoys: upstream APIs filter out decoys (`:target == false`) before calling this. If a `:target` column is present, rows with `target == false` are dropped defensively.

Required columns
- protein_results: `:protein`, `score_col`, and either `:ms_file_idx` or `:file_name`.
"""
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
    protein_results_copy = if hasproperty(protein_results, :target)
        protein_results[protein_results.target .== true, :]
    else
        copy(protein_results)
    end
    # Build grouping keys: [:species]? , :protein , [:entrap_id]?
    keys = Symbol[]
    if hasproperty(protein_results_copy, :species)
        push!(keys, :species)
    end
    push!(keys, :protein)
    if hasproperty(protein_results_copy, :entrapment_group_id)
        push!(keys, :entrapment_group_id)
    elseif hasproperty(protein_results_copy, :entrap_id)
        push!(keys, :entrap_id)
    end
    grouped = groupby(protein_results_copy, keys)
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
