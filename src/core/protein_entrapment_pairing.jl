"""
Protein-level entrapment pairing functions for EFDR analysis.
"""

function assign_protein_entrapment_pairs!(df::DataFrame)
    required_cols = [:protein, :entrap_id]
    missing_cols = [col for col in required_cols if !hasproperty(df, col)]
    if !isempty(missing_cols)
        error("DataFrame missing required columns: $missing_cols")
    end
    if !hasproperty(df, :entrap_pair_id)
        df[!, :entrap_pair_id] = Vector{Union{Missing, UInt32}}(missing, nrow(df))
    end
    pair_counter = UInt32(1)
    grouped = groupby(df, :protein)
    for (_, group_df) in pairs(grouped)
        group_0_indices = findall(group_df.entrap_id .== 0)
        other_indices = findall(group_df.entrap_id .!= 0)
        if isempty(group_0_indices) || isempty(other_indices)
            continue
        end
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
            group_df[idx_0, :entrap_pair_id] = pair_counter
            for entrap_id in unique_entrap_groups
                group_member_indices = entrap_groups_dict[entrap_id]
                member_pos = ((idx_0_pos - 1) % length(group_member_indices)) + 1
                member_idx = group_member_indices[member_pos]
                group_df[member_idx, :entrap_pair_id] = pair_counter
            end
            pair_counter += UInt32(1)
        end
    end
    return nothing
end

function add_protein_entrap_pair_ids!(protein_results::DataFrame, protein_library::DataFrame)
    if !hasproperty(protein_results, :protein)
        error("protein_results must have :protein column")
    end
    if !hasproperty(protein_library, :entrap_pair_id)
        error("protein_library must have :entrap_pair_id column. Run assign_protein_entrapment_pairs! first.")
    end
    if !hasproperty(protein_library, :protein)
        error("protein_library must have :protein column")
    end
    protein_to_pair = Dict{String, Union{Missing, UInt32}}()
    for row in eachrow(protein_library)
        if !haskey(protein_to_pair, row.protein)
            protein_to_pair[row.protein] = row.entrap_pair_id
        end
    end
    protein_results[!, :entrap_pair_id] = [get(protein_to_pair, protein, missing) for protein in protein_results.protein]
    return nothing
end

