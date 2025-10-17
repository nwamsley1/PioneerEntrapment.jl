using DataFrames
using Dictionaries

"""
    get_best_protein_representative(protein_group::String, peptide_indices::Union{Missing, Vector{UInt32}},
                                     precursors_library::DataFrame) -> String

Given a semicolon-delimited protein group and its supporting peptides, return the
single protein that appears most frequently in those peptides' protein_groups.

Arguments
- protein_group: e.g., "P12345;Q67890;O11111"
- peptide_indices: Vector of precursor_idx values from :peptides column (can contain missing)
- precursors_library: Must have :protein_groups column (or :protein, depending on schema)

Returns
- Single protein ID with maximum support, or first protein if library lookup fails

Algorithm
1. Parse protein_group into candidate proteins
2. For each peptide in peptide_indices:
   - Look up peptide's :protein_groups in precursors_library
   - For each candidate protein, check if it appears in peptide's groups
   - Increment count for matching candidates
3. Return candidate with highest count

Edge Cases
- Single-protein groups: Return immediately
- Empty peptide_indices: Fall back to first protein
- Missing peptide_indices: Fall back to first protein
- Invalid precursor_idx: Skip and continue
- Tie in counts: Use first occurrence order
"""
function get_best_protein_representative(protein_group::String,
                                          peptide_indices::Union{Missing, Vector{UInt32}},
                                          precursors_library::DataFrame)
    # Split the group into candidate proteins
    candidates = split(protein_group, ';')

    # Handle single-protein groups
    if length(candidates) == 1
        return String(strip(candidates[1]))
    end

    # Handle missing/empty peptides - fall back to first
    if ismissing(peptide_indices) || isempty(peptide_indices)
        return String(strip(candidates[1]))
    end

    # Initialize count dictionary for each candidate
    protein_counts = Dict{String, Int}()
    for protein in candidates
        protein_counts[String(strip(protein))] = 0
    end

    # Determine column name (check :accession_numbers first, then :protein_groups, then :protein)
    protein_col = if hasproperty(precursors_library, :accession_numbers)
        :accession_numbers
    elseif hasproperty(precursors_library, :protein_groups)
        :protein_groups
    elseif hasproperty(precursors_library, :protein)
        :protein
    else
        @warn "Library missing :accession_numbers, :protein_groups, or :protein column, falling back to first protein"
        return String(strip(candidates[1]))
    end

    # Query each peptide's protein_groups in library
    for prec_idx in peptide_indices
        # Skip missing entries in peptide_indices
        if ismissing(prec_idx)
            continue
        end

        # Validate index bounds
        if prec_idx < 1 || prec_idx > nrow(precursors_library)
            continue
        end

        # Get peptide's protein group annotation
        peptide_protein_group = precursors_library[prec_idx, protein_col]

        # Handle missing values
        if ismissing(peptide_protein_group)
            continue
        end

        peptide_proteins = split(String(peptide_protein_group), ';')

        # Increment count for any candidate found in this peptide
        for pep_prot in peptide_proteins
            pep_prot_clean = String(strip(pep_prot))
            if haskey(protein_counts, pep_prot_clean)
                protein_counts[pep_prot_clean] += 1
            end
        end
    end

    # Return the protein with maximum count (ties broken by first occurrence)
    best_protein = candidates[1]  # default fallback
    max_count = 0
    for protein in candidates  # iterate in order to break ties consistently
        count = protein_counts[String(strip(protein))]
        if count > max_count
            max_count = count
            best_protein = protein
        end
    end

    return String(strip(best_protein))
end

"""
    add_original_target_protein_scores!(protein_results::DataFrame, precursors_library::DataFrame; score_col=:pg_score) -> Nothing

Add a `"<score>_original_target"` column capturing, for each protein row,
the score of the original target protein (entrap_id == 0) with the same
`file` and `species`. For protein group identifiers that contain multiple
entries separated by `;`, the protein with the most peptide evidence is
used as the group key.

If the row itself is the original target, its own score is used. If there is no
matching original target, the value is set to `-1.0f0`.

Required columns
- protein_results: `:protein`, `:entrap_id`, `:species`, `:peptides`, `score_col`, and either `:ms_file_idx` or `:file_name`.
- precursors_library: DataFrame with `:protein_groups` or `:protein` column for evidence-based selection.
"""
function add_original_target_protein_scores!(protein_results::DataFrame, precursors_library::DataFrame; score_col=:pg_score)
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
    # For global results (file="global"), use full protein string. For per-file, use first from semicolon-split.
    protein_to_target = Dictionary{Tuple{Any, String, String}, Float32}()
    for row in eachrow(protein_results)
        if row.entrap_id == 0 && !ismissing(row[score_col])
            species = String(row.species)
            # Global results: use full protein string. Per-file: use evidence-based selection.
            is_global = row[file_col] == "global" || (row[file_col] isa Integer && row[file_col] == 0)
            protkey = if is_global
                String(row.protein)
            else
                # Use evidence-based selection
                peptide_indices = hasproperty(protein_results, :peptides) ? row.peptides : missing
                get_best_protein_representative(String(row.protein), peptide_indices, precursors_library)
            end
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
                # Global results: use full protein string. Per-file: use evidence-based selection.
                is_global = row[file_col] == "global" || (row[file_col] isa Integer && row[file_col] == 0)
                protkey = if is_global
                    String(row.protein)
                else
                    # Use evidence-based selection
                    peptide_indices = hasproperty(protein_results, :peptides) ? row.peptides : missing
                    get_best_protein_representative(String(row.protein), peptide_indices, precursors_library)
                end
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
    add_original_target_protein_scores!(protein_results::DataFrame, precursors_library::DataFrame, score_cols::Vector{Symbol}) -> Nothing

Vectorized overload to compute original target columns for multiple score fields.
"""
function add_original_target_protein_scores!(protein_results::DataFrame, precursors_library::DataFrame, score_cols::Vector{Symbol})
    for score_col in score_cols
        add_original_target_protein_scores!(protein_results, precursors_library; score_col=score_col)
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
