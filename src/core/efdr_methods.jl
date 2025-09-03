using DataFrames

# Base type for EFDR methods
abstract type EFDRMethod end

# EFDR method structs with data fields
struct CombinedEFDR{T<:Real, I<:Integer} <: EFDRMethod
    score::Vector{T}
    original_target_score::Vector{T}
    entrapment_label::Vector{I}
    qval::Vector{T}
    r::T
    function CombinedEFDR(score::Vector{T}, original_target_score::Vector{T}, entrapment_label::Vector{I}, qval::Vector{T}, r::T) where {T<:Real, I<:Integer}
        n = length(score)
        if length(original_target_score) != n || length(entrapment_label) != n || length(qval) != n
            error("All input vectors must have the same length")
        end
        new{T,I}(score, original_target_score, entrapment_label, qval, r)
    end
end

struct PairedEFDR{T<:Real, I<:Integer} <: EFDRMethod
    score::Vector{T}
    original_target_score::Vector{T}
    entrapment_label::Vector{I}
    qval::Vector{T}
    r::T
    function PairedEFDR(score::Vector{T}, original_target_score::Vector{T}, entrapment_label::Vector{I}, qval::Vector{T}, r::T) where {T<:Real, I<:Integer}
        n = length(score)
        if length(original_target_score) != n || length(entrapment_label) != n || length(qval) != n
            error("All input vectors must have the same length")
        end
        new{T,I}(score, original_target_score, entrapment_label, qval, r)
    end
end

# Convenience constructors
CombinedEFDR(score, original_target_score, entrapment_label, qval; r=1.0) = CombinedEFDR(convert(Vector{Float64}, score), convert(Vector{Float64}, original_target_score), convert(Vector{Int}, entrapment_label), convert(Vector{Float64}, qval), Float64(r))
PairedEFDR(score, original_target_score, entrapment_label, qval; r=1.0) = PairedEFDR(convert(Vector{Float64}, score), convert(Vector{Float64}, original_target_score), convert(Vector{Int}, entrapment_label), convert(Vector{Float64}, qval), Float64(r))

"""
    calculate_efdr(method::EFDRMethod)

Calculate empirical FDR using the data stored in the method struct.

Methods:
- CombinedEFDR: Standard combined empirical FDR calculation
- PairedEFDR: Paired empirical FDR that considers score relationships

Returns a vector of EFDR values aligned to the input order.
"""
function calculate_efdr(method::CombinedEFDR)
    # Sort by q-value (ascending) first, then by score (descending) to break ties
    sort_indices = sortperm(collect(zip(-method.score, method.qval)))

    entrapment_fdr = zeros(eltype(method.qval), length(method.qval))
    Nτ = 0
    Nϵ = 0
    for i in sort_indices
        # Determine if this is a target or entrapment
        is_original_target = iszero(method.entrapment_label[i])

        if is_original_target
            Nτ += 1
        else
            Nϵ += 1
        end

        # Avoid division by zero and cap at 1.0
        if Nϵ + Nτ > 0
            entrapment_fdr[i] = min(1.0, (Nϵ * (1 + 1/method.r)) / (Nϵ + Nτ))
        else
            entrapment_fdr[i] = 0.0
        end
    end

    # Enforce monotonicity over the sorted order (cumulative minimum in reverse)
    fdr = Inf
    @inbounds @fastmath for i in reverse(sort_indices)
        if entrapment_fdr[i] > fdr
            entrapment_fdr[i] = fdr
        else
            fdr = entrapment_fdr[i]
        end
    end

    return entrapment_fdr
end

function calculate_efdr(method::PairedEFDR; stride::Int=5)
    # Sort by q-value (ascending) first, then by score (descending) to break ties
    sort_indices = sortperm(collect(zip(-method.score, method.qval)))

    entrapment_fdr = zeros(eltype(method.qval), length(method.qval))
    n = length(sort_indices)

    # Normalize stride: 1 means compute at every index (full O(n^2)); >1 samples
    stride = stride <= 0 ? 1 : stride
    ks = collect(1:stride:n)
    if ks[end] != n
        push!(ks, n)
    end

    total_ops = sum(ks) # approximate work for progress bar
    pb = ProgressBar(total=total_ops)
    completed_ops = 0
    last_update = 0

    last_k = 0
    last_value = zero(eltype(method.qval))
    for k in ks
        Nτ = 0
        Nϵ = 0
        Nϵsτ = 0
        Nϵτs = 0
        s = method.score[sort_indices[k]]
        for j in 1:k
            # Determine if this is a target or entrapment
            is_original_target = method.entrapment_label[sort_indices[j]] == 0
            if is_original_target
                Nτ += 1
            else
                t = method.original_target_score[sort_indices[j]]
                e = method.score[sort_indices[j]]
                Nϵ += 1
                if (e >= s) & (t < s)
                    Nϵsτ += 1
                elseif (e > t) & (t >= s)
                    Nϵτs += 1
                end
            end
            completed_ops += 1
            if completed_ops % 1000 == 0
                update(pb, completed_ops - last_update)
                last_update = completed_ops
            end
        end

        # Compute EFDR at sampled position
        value = if Nϵ + Nτ > 0
            min(1.0, (Nϵ + Nϵsτ + 2*Nϵτs) / (Nϵ + Nτ))
        else
            0.0
        end

        # Fill forward for unsampled indices since last checkpoint
        for kk in (last_k + 1):k
            entrapment_fdr[sort_indices[kk]] = value
        end
        last_k = k
        last_value = value
    end

    # Enforce monotonicity over the sorted order (cumulative minimum in reverse)
    fdr = Inf
    @inbounds @fastmath for i in reverse(sort_indices)
        if entrapment_fdr[i] > fdr
            entrapment_fdr[i] = fdr
        else
            fdr = entrapment_fdr[i]
        end
    end

    return entrapment_fdr
end

function get_combined_efdr(score::AbstractVector{<:Real}, original_target_score::AbstractVector{<:Real}, entrapment_label::AbstractVector{<:Integer}, qval::AbstractVector{<:Real}, r::Real = 1.0)
    method = CombinedEFDR(score, original_target_score, entrapment_label, qval; r=r)
    return calculate_efdr(method)
end

function get_paired_efdr(score::AbstractVector{<:Real}, original_target_score::AbstractVector{<:Real}, entrapment_label::AbstractVector{<:Integer}, qval::AbstractVector{<:Real}, r::Real = 1.0)
    method = PairedEFDR(score, original_target_score, entrapment_label, qval; r=r)
    return calculate_efdr(method)
end

function add_efdr_columns!(df::DataFrame, library_precursors::DataFrame; method_types::Vector=[CombinedEFDR, PairedEFDR], score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:score, :qval)], r::Float64=1.0, paired_stride::Int=5)
    if !hasproperty(df, :precursor_idx)
        error("DataFrame must have :precursor_idx column")
    end
    if !hasproperty(library_precursors, :entrapment_group_id)
        error("library_precursors must have :entrapment_group_id column")
    end
    entrap_labels = [library_precursors.entrapment_group_id[pid] for pid in df.precursor_idx]
    for (score_col, qval_col) in score_qval_pairs
        if !hasproperty(df, score_col)
            @warn "Column $score_col not found in DataFrame, skipping..."; continue
        end
        if !hasproperty(df, qval_col)
            @warn "Column $qval_col not found in DataFrame, skipping..."; continue
        end
        original_target_col = Symbol(String(score_col) * "_original_target")
        if !hasproperty(df, original_target_col)
            @warn "Column $original_target_col not found. Make sure to run add_original_target_scores! first."; continue
        end
        scores = Float64.(df[!, score_col])
        original_target_scores = Float64.(df[!, original_target_col])
        qvals = Float64.(df[!, qval_col])
        for method_type in method_types
            method_name = method_type == CombinedEFDR ? "combined" : method_type == PairedEFDR ? "paired" : lowercase(string(method_type))
            efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
            method = method_type(scores, original_target_scores, entrap_labels, qvals, r)
            efdr_values = method_type == PairedEFDR ? calculate_efdr(method; stride=paired_stride) : calculate_efdr(method)
            df[!, efdr_col] = efdr_values
        end
    end
    return nothing
end
