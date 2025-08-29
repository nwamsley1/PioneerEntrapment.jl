"""
    getModKey(mod_string::AbstractString)
"""
function getModKey(mod_string::AbstractString)
    mod_pattern = r"\(\d+,[A-Z],([^)]+)\)"
    mod_names = String[]
    for match in eachmatch(mod_pattern, mod_string)
        push!(mod_names, match.captures[1])
    end
    sort!(mod_names)
    return join(mod_names, ";")
end

"""
    assign_entrapment_pairs!(df::DataFrame)
"""
function assign_entrapment_pairs!(df::DataFrame)
    if !hasproperty(df, :entrap_pair_id)
        df[!, :entrap_pair_id] = Vector{Union{Missing, UInt32}}(missing, nrow(df))
    end
    pair_counter = UInt32(1)
    grouped = groupby(df, [:base_pep_id, :prec_charge, :is_decoy, :mod_key])
    for (_, group_df) in pairs(grouped)
        group_0_indices = findall(group_df.entrapment_group_id .== 0)
        other_indices = findall(group_df.entrapment_group_id .!= 0)
        if isempty(group_0_indices) || isempty(other_indices)
            continue
        end
        entrap_groups_dict = Dict{Int, Vector{Int}}()
        for idx in other_indices
            entrap_id = group_df.entrapment_group_id[idx]
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

"""
    add_entrap_pair_ids!(prec_results::DataFrame, library_precursors::DataFrame)
"""
function add_entrap_pair_ids!(prec_results::DataFrame, library_precursors::DataFrame)
    if !hasproperty(prec_results, :precursor_idx)
        error("prec_results must have :precursor_idx column")
    end
    if !hasproperty(library_precursors, :entrap_pair_id)
        error("library_precursors must have :entrap_pair_id column")
    end
    prec_results[!, :entrap_pair_id] = [library_precursors[!, :entrap_pair_id][pid] for pid in prec_results[!, :precursor_idx]]
    return nothing
end

