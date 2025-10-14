"""
    getModKey(mod_string::AbstractString) -> String

Extract a canonical modification key from a modification annotation string.

The function searches for substrings of the form `(pos,AA,MODNAME)` and returns
the sorted, semicolon-delimited list of `MODNAME` values. This allows grouping
precursors that differ only by order or position of identical modifications.

Arguments
- mod_string: Freeform modification annotation containing tokens like `(5,M,Oxidation)`.

Returns
- A string such as `"Acetyl;Oxidation"`. Returns `""` if no mods are found.

Example
    julia> getModKey("(5,M,Oxidation)(1,K,Acetyl)")
    "Acetyl;Oxidation"
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

## assign_entrapment_pairs! removed — pairing must be precomputed and provided as :entrapment_pair_id.

"""
    add_entrap_pair_ids!(prec_results::DataFrame, library_precursors::DataFrame) -> Nothing

Attach precomputed entrapment pair IDs to a precursor-level results table.

The `library_precursors` table must contain `:entrapment_pair_id` (precomputed during
library construction). This function copies that ID onto each row of
`prec_results` by indexing via `:precursor_idx`.

Required columns
- prec_results: `:precursor_idx`
- library_precursors: `:entrapment_pair_id`

Side effects
- Adds a `:entrapment_pair_id::Union{Missing,UInt32}` column to `prec_results`.
"""
function add_entrap_pair_ids!(prec_results::DataFrame, library_precursors::DataFrame)
    if !hasproperty(prec_results, :precursor_idx)
        error("prec_results must have :precursor_idx column")
    end
    if !hasproperty(library_precursors, :entrapment_pair_id)
        error("library_precursors must have :entrapment_pair_id column")
    end
    # Map each results row to its precursor in the library and carry over the pair ID.
    prec_results[!, :entrapment_pair_id] = [library_precursors[!, :entrapment_pair_id][pid] for pid in prec_results[!, :precursor_idx]]
    return nothing
end
