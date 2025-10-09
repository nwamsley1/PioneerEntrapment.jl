"""
Protein-level entrapment pairing helpers.

These functions define how proteins are paired across target/entrapment species
so EFDR can be computed at the protein level in a manner analogous to
precursor-level pairing.
"""

"""
    assign_protein_entrapment_pairs!(df::DataFrame) -> Nothing

Assign a stable `:entrapment_pair_id` for each protein across target and entrapment entries.

Rules (per unique `:protein`):
- Split rows by `:entrap_id == 0` (original target) vs `:entrap_id != 0` (entrapment groups).
- If a protein lacks at least one target AND one entrapment, it is skipped.
- For each target row, assign a new pair ID and round‑robin assign that same
  ID to exactly one row from each entrapment group for that protein. This avoids
  bias when there are unequal numbers of entries per group.

Required columns
- df: `:protein`, `:entrap_id`

Side effects
- Adds/fills `:entrapment_pair_id::Union{Missing,UInt32}` in `df`.
"""
function assign_protein_entrapment_pairs!(df::DataFrame)
    required_cols = [:protein, :entrap_id]
    missing_cols = [col for col in required_cols if !hasproperty(df, col)]
    if !isempty(missing_cols)
        error("DataFrame missing required columns: $missing_cols")
    end
    if !hasproperty(df, :entrapment_pair_id)
        df[!, :entrapment_pair_id] = Vector{Union{Missing, UInt32}}(missing, nrow(df))
    end
    pair_counter = UInt32(1)
    grouped = groupby(df, :protein)
    for (_, group_df) in pairs(grouped)
        # Partition the group into original target rows (entrap_id=0) and others
        group_0_indices = findall(group_df.entrap_id .== 0)
        other_indices = findall(group_df.entrap_id .!= 0)
        if isempty(group_0_indices) || isempty(other_indices)
            continue
        end
        # Bucket entrapment rows by entrapment group ID
        entrap_groups_dict = Dict{UInt8, Vector{Int}}()
        for idx in other_indices
            entrap_id = group_df.entrap_id[idx]
            if !haskey(entrap_groups_dict, entrap_id)
                entrap_groups_dict[entrap_id] = Int[]
            end
            push!(entrap_groups_dict[entrap_id], idx)
        end
        unique_entrap_groups = sort(collect(keys(entrap_groups_dict)))
        for (idx_0_pos, idx_0) in enumerate(group_0_indices)
            # Anchor the pair by the position of the original target row
            group_df[idx_0, :entrapment_pair_id] = pair_counter
            for entrap_id in unique_entrap_groups
                group_member_indices = entrap_groups_dict[entrap_id]
                # Round‑robin assign to handle unequal counts across groups
                member_pos = ((idx_0_pos - 1) % length(group_member_indices)) + 1
                member_idx = group_member_indices[member_pos]
                group_df[member_idx, :entrapment_pair_id] = pair_counter
            end
            pair_counter += UInt32(1)
        end
    end
    return nothing
end

"""
    add_protein_entrap_pair_ids!(protein_results::DataFrame, protein_library::DataFrame) -> Nothing

Propagate preassigned `:entrapment_pair_id` from a protein library table to a results table.

This maps on `:protein`. If the library contains multiple rows per protein, the
first encountered `:entrapment_pair_id` is used.

Required columns
- protein_results: `:protein`
- protein_library: `:protein`, `:entrapment_pair_id`

Side effects
- Adds `:entrapment_pair_id` column to `protein_results`.
"""
function add_protein_entrap_pair_ids!(protein_results::DataFrame, protein_library::DataFrame)
    if !hasproperty(protein_results, :protein)
        error("protein_results must have :protein column")
    end
    if !hasproperty(protein_library, :entrapment_pair_id)
        error("protein_library must have :entrapment_pair_id column. Run assign_protein_entrapment_pairs! first.")
    end
    if !hasproperty(protein_library, :protein)
        error("protein_library must have :protein column")
    end
    # Build a mapping from protein -> first seen entrapment_pair_id
    protein_to_pair = Dict{String, Union{Missing, UInt32}}()
    for row in eachrow(protein_library)
        if !haskey(protein_to_pair, row.protein)
            protein_to_pair[row.protein] = row.entrapment_pair_id
        end
    end
    protein_results[!, :entrapment_pair_id] = [get(protein_to_pair, protein, missing) for protein in protein_results.protein]
    return nothing
end
