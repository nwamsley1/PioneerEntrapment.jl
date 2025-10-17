using DataFrames
using Dictionaries

"""
    add_original_target_protein_scores!(protein_results::DataFrame; score_col=:pg_score) -> Nothing

Add a `"<score>_original_target"` column capturing, for each protein row,
the score of the original target protein (entrap_id == 0) that has the maximum
protein overlap with the entrapment group.

Pairing Strategy:
For each entrapment protein group, find the target group (in the same file and species)
that shares the most protein accessions. For example:
- Entrapment: "PROT1;PROT2;PROT3"
- Target candidates: "PROT1" (overlap=1) vs "PROT2;PROT3;PROT4" (overlap=2)
- Result: Pair with "PROT2;PROT3;PROT4" (best match)

For global results (file="global"), the full protein string is used as a unique key
for exact matching.

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

    # Build index of target groups per (file, species)
    # Structure: Dict{(file_id, species), Vector{(row_idx, protein_group_string, protein_set, score)}}
    target_index = Dict{Tuple{Any, String}, Vector{Tuple{Int, String, Set{String}, Float32}}}()

    for (idx, row) in enumerate(eachrow(protein_results))
        if row.entrap_id == 0 && !ismissing(row[score_col])
            file_id = row[file_col]
            species = String(row.species)
            protein_group = String(row.protein)

            # For global results, use full protein string; for per-file, parse into set
            is_global = file_id == "global" || (file_id isa Integer && file_id == 0)

            if is_global
                # Global: exact match only, store full string as single-element set
                protein_set = Set([protein_group])
            else
                # Per-file: parse semicolon-delimited proteins into set
                protein_set = Set(String(strip(p)) for p in split(protein_group, ';'))
            end

            key = (file_id, species)
            if !haskey(target_index, key)
                target_index[key] = Tuple{Int, String, Set{String}, Float32}[]
            end

            push!(target_index[key], (idx, protein_group, protein_set, Float32(row[score_col])))
        end
    end

    # For each row, find best matching target
    original_target_scores = Float32[]

    for row in eachrow(protein_results)
        if ismissing(row[score_col])
            push!(original_target_scores, -1.0f0)
            continue
        end

        if row.entrap_id == 0
            # Target rows pair with themselves
            push!(original_target_scores, Float32(row[score_col]))
        else
            # Entrapment rows: find best matching target
            file_id = row[file_col]
            species = String(row.species)
            protein_group = String(row.protein)

            is_global = file_id == "global" || (file_id isa Integer && file_id == 0)

            key = (file_id, species)
            candidates = get(target_index, key, Tuple{Int, String, Set{String}, Float32}[])

            if isempty(candidates)
                # No targets in this file/species
                push!(original_target_scores, -1.0f0)
            elseif is_global
                # Global: exact string match
                match_idx = findfirst(c -> c[2] == protein_group, candidates)
                if match_idx !== nothing
                    push!(original_target_scores, candidates[match_idx][4])
                else
                    push!(original_target_scores, -1.0f0)
                end
            else
                # Per-file: find target with maximum protein overlap
                entrap_proteins = Set(String(strip(p)) for p in split(protein_group, ';'))

                best_score = -1.0f0
                max_overlap = 0

                for (_, target_group, target_set, target_score) in candidates
                    overlap = length(intersect(entrap_proteins, target_set))
                    if overlap > max_overlap
                        max_overlap = overlap
                        best_score = target_score
                    end
                end

                push!(original_target_scores, best_score)
            end
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
