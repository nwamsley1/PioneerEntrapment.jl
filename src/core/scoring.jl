"""
    get_complement_score(prec_results::DataFrame, row_idx::Int; score_col=:score) -> Float64

Return the best score of the complementary entry within the same entrapment pair.

Given a row index, this scans all rows with the same `:entrapment_pair_id` and
returns the maximum score among those rows excluding `row_idx` itself. If the
pair ID is missing or no complement exists, returns `0.0`.

Required columns
- prec_results: `:entrapment_pair_id`, `score_col`
"""
function get_complement_score(prec_results::DataFrame, row_idx::Int; score_col=:score)
    if !hasproperty(prec_results, :entrapment_pair_id)
        error("DataFrame must have :entrapment_pair_id column. Run add_entrap_pair_ids! first.")
    end
    if !hasproperty(prec_results, score_col)
        error("DataFrame must have :$score_col column.")
    end
    if row_idx < 1 || row_idx > nrow(prec_results)
        error("row_idx $row_idx is out of bounds. DataFrame has $(nrow(prec_results)) rows.")
    end
    pair_id = prec_results[row_idx, :entrapment_pair_id]
    if ismissing(pair_id)
        return 0.0
    end
    same_pair = coalesce.(prec_results.entrapment_pair_id .== pair_id, false)
    not_self = (1:nrow(prec_results)) .!= row_idx
    complement_mask = same_pair .& not_self
    complement_scores = prec_results[complement_mask, score_col]
    return isempty(complement_scores) ? 0.0 : maximum(skipmissing(complement_scores), init=0.0)
end

"""
    add_original_target_scores!(prec_results::DataFrame, library_precursors::DataFrame; score_col=:score) -> Nothing

Attach the original target score for each row based on entrapment pairing.

For each row in `prec_results` with `:entrapment_pair_id`, we determine the
original target (the library entry with `:entrapment_group_id == 0`) within the
same pair and MS file, and write that score to a new column
`Symbol(String(score_col) * "_original_target")`.

When a row itself is the original target, the new column equals its own score.
When the original target does not exist (e.g., missing partner), the value is
set to `-1.0` to flag absence.

Required columns
- prec_results: `:entrapment_pair_id`, `:precursor_idx`, `score_col`, and either `:ms_file_idx` or implicitly use 0
- library_precursors: `:entrapment_group_id`
"""
function add_original_target_scores!(prec_results::DataFrame, library_precursors::DataFrame; score_col=:score)
    if !hasproperty(prec_results, :entrapment_pair_id)
        error("DataFrame must have :entrapment_pair_id column. Run add_entrap_pair_ids! first.")
    end
    if !hasproperty(prec_results, score_col)
        error("DataFrame must have :$score_col column.")
    end
    if !hasproperty(prec_results, :precursor_idx)
        error("DataFrame must have :precursor_idx column.")
    end
    if !hasproperty(library_precursors, :entrapment_group_id)
        error("library_precursors must have :entrapment_group_id column.")
    end
    original_target_col = Symbol(String(score_col) * "_original_target")
    # Map (ms_file_idx, entrapment_pair_id) -> (target_row_idx, target_score)
    pair_to_target = Dictionary{Int, Dictionary{UInt32, @NamedTuple{target_row::UInt32,target_score::Float64}}}()
    for (idx, row) in enumerate(eachrow(prec_results))
        if !ismissing(row.entrapment_pair_id) && !ismissing(row[score_col])
            pair_id = row.entrapment_pair_id
            precursor_idx = row.precursor_idx
            ms_file_idx = hasproperty(prec_results, :ms_file_idx) ? row.ms_file_idx : 0
            entrap_group = library_precursors.entrapment_group_id[precursor_idx]
            if entrap_group == 0
                if !haskey(pair_to_target, ms_file_idx)
                    insert!(pair_to_target, ms_file_idx, Dictionary{UInt32, @NamedTuple{target_row::UInt32,target_score::Float32}}())
                end
                insert!(pair_to_target[ms_file_idx], pair_id, (target_row = UInt32(idx), target_score = Float64(row[score_col])))
            end
        end
    end
    original_target_scores = zeros(Float64, nrow(prec_results))
    for (idx, row) in enumerate(eachrow(prec_results))
        if !ismissing(row.entrapment_pair_id) && !ismissing(row[score_col])
            pair_id = row.entrapment_pair_id
            precursor_idx = row.precursor_idx
            ms_file_idx = hasproperty(prec_results, :ms_file_idx) ? row.ms_file_idx : 0
            entrap_group = library_precursors.entrapment_group_id[precursor_idx]
            if entrap_group == 0
                original_target_scores[idx] = row[score_col]
            else
                if haskey(pair_to_target, ms_file_idx) && haskey(pair_to_target[ms_file_idx], pair_id)
                    _, target_score = pair_to_target[ms_file_idx][pair_id]
                    original_target_scores[idx] = target_score
                else
                    original_target_scores[idx] = -1.0
                end
            end
        end
    end
    prec_results[!, original_target_col] = original_target_scores
    return nothing
end

function add_original_target_scores!(prec_results::DataFrame, library_precursors::DataFrame, score_cols::Vector{Symbol})
    for score_col in score_cols
        add_original_target_scores!(prec_results, library_precursors; score_col=score_col)
    end
    return nothing
end

"""
    create_global_results_df(prec_results::DataFrame; score_col::Symbol=:global_prob) -> DataFrame

Build a per-precursor global table by selecting, for each `:precursor_idx`, the
row with the best available `score_col` across files.

If `:ms_file_idx` exists, the returned rows have `:ms_file_idx = 0`. If only
`:file_name` exists, a synthetic `"global"` file name is assigned.

Required columns
- prec_results: `:precursor_idx`, `:ms_file_idx` (or `:file_name`), `score_col`
"""
function create_global_results_df(prec_results::DataFrame; score_col::Symbol=:global_prob)
    required_cols = [:precursor_idx, :ms_file_idx, score_col]
    missing_cols = [col for col in required_cols if !hasproperty(prec_results, col)]
    if !isempty(missing_cols)
        error("DataFrame missing required columns: $missing_cols")
    end
    prec_results_copy = copy(prec_results)
    grouped = groupby(prec_results_copy, :precursor_idx)
    global_df = combine(grouped) do group
        valid_rows = group[.!ismissing.(group[!, score_col]), :]
        if nrow(valid_rows) == 0
            return similar(group, 0)
        end
        best_idx = argmax(valid_rows[!, score_col])
        return valid_rows[best_idx:best_idx, :]
    end
    if nrow(global_df) > 0
        global_df[!, :ms_file_idx] .= 0
    end
    return global_df
end
