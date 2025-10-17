"""
    calculate_efdr_fast(method::PairedEFDR; cuts_mode::Symbol=:targets_only)

Compute paired EFDR using an O(n log n + K) approach that aggregates counts
over score cutoffs rather than scanning quadratically per threshold.

cuts_mode controls which unique score cutoffs are used:
- :targets_only => unique(t) where t = original_target_score (non-missing)
- :all          => unique union of t and e (entrap scores)
"""
function calculate_efdr_fast(method::PairedEFDR; cuts_mode::Symbol=:targets_only)
    e = method.score
    t = method.original_target_score
    # Allow missing targets; treat them as NA in counting logic
    t_miss = isnan.(fill(NaN, length(t)))  # placeholder vector for type inference
    t_anymissing = any(ismissing, t)

    # Build cuts ascending
    cuts = if cuts_mode === :targets_only
        unique(skipmissing(t)) |> collect |> sort!
    elseif cuts_mode === :all
        vcat(collect(skipmissing(t)), collect(skipmissing(e))) |> unique! |> sort!
    else
        error("Unsupported cuts_mode: $cuts_mode")
    end
    K = length(cuts)
    if K == 0
        return zeros(Float64, length(e))
    end

    # Helper: interval index 0..K where idx = max{ j : cuts[j] <= value }
    idx_t_all = Vector{Union{Int,Missing}}(undef, length(t))
    idx_e_all = Vector{Union{Int,Missing}}(undef, length(e))
    @inbounds for i in eachindex(t)
        ti = t[i]
        if ismissing(ti)
            idx_t_all[i] = missing
        else
            idx_t_all[i] = searchsortedlast(cuts, ti)
        end
        ei = e[i]
        if ismissing(ei)
            idx_e_all[i] = missing
        else
            idx_e_all[i] = searchsortedlast(cuts, ei)
        end
    end

    # Utility: suffix sum over length K
    function suffsum!(arr)
        s = 0
        @inbounds for i in K:-1:1
            s += arr[i]
            arr[i] = s
        end
        return arr
    end

    # Utility: tabulate indices 1..K (ignore <=0 and missing)
    function ztab(idxs)
        freq = zeros(Int, K)
        @inbounds for v in idxs
            if v isa Int
                iv = v::Int
                if iv > 0
                    freq[iv] += 1
                end
            end
        end
        return freq
    end

    # n_targets(s): # { target rows with score >= s } for s = cuts[j]
    # Filter by entrapment_label == 0 (target rows) AND use their score (e)
    target_mask = method.entrapment_label .== 0
    idx_e_targets = idx_e_all[target_mask]
    idx_e_targets_valid = filter(!ismissing, idx_e_targets)
    T_asc = suffsum!(ztab(idx_e_targets_valid))

    # n_e(s): # { entrapment rows with score >= s }
    # Filter by entrapment_label != 0 (entrapment rows) AND use their score (e)
    entrap_mask = method.entrapment_label .!= 0
    idx_e_entraps = idx_e_all[entrap_mask]
    idx_e_entraps_valid = filter(!ismissing, idx_e_entraps)
    E_asc = suffsum!(ztab(idx_e_entraps_valid))

    # n_ets(s): entrap >= s AND target >= s AND entrap > target
    # Only consider entrapment rows
    both_mask = entrap_mask .& .!(ismissing.(t)) .& .!(ismissing.(e))
    ETS_asc = zeros(Int, K)
    if any(both_mask)
        idx_t_b = idx_t_all[both_mask]
        idx_e_b = idx_e_all[both_mask]
        gt = (e[both_mask] .> t[both_mask])
        if any(gt)
            min_idx = map((a,b)->begin
                if ismissing(a) || ismissing(b)
                    0
                else
                    min(a::Int,b::Int)
                end
            end, idx_t_b[gt], idx_e_b[gt])
            # suffix sum over min_idx
            ETS_asc = suffsum!(ztab(min_idx))
        end
    end

    # n_est(s): entrap >= s AND (target < s OR target is NA)
    # Only consider entrapment rows
    diff = zeros(Int, K+1)
    # target is NA & e present (entrapment rows only)
    e_only = entrap_mask .& ismissing.(t) .& .!(ismissing.(e))
    if any(e_only)
        idx_e_na = filter(x -> x isa Int && x > 0, idx_e_all[e_only])
        if !isempty(idx_e_na)
            # Each contributes to [1 .. idx_e]
            diff[1] += length(idx_e_na)
            # subtract 1 at (idx_e+1)
            offs = zeros(Int, K+1)
            @inbounds for v in idx_e_na
                offs[(v::Int)+1] += 1
            end
            diff .-= offs
        end
    end
    # both present: start = idx_t + 1 .. end = idx_e (entrapment rows only)
    # Reuse both_mask which already includes entrap_mask
    if any(both_mask)
        bt = idx_t_all[both_mask]
        be = idx_e_all[both_mask]
        starts = Int[]; ends = Int[]
        @inbounds for i in eachindex(bt)
            it = bt[i]; ie = be[i]
            if !(ismissing(it) || ismissing(ie))
                st = (it::Int) + 1
                en = (ie::Int)
                if en > 0 && st <= en
                    push!(starts, st); push!(ends, en)
                end
            end
        end
        if !isempty(starts)
            @inbounds for st in starts
                diff[st] += 1
            end
            @inbounds for en in ends
                diff[en+1] -= 1
            end
        end
    end
    EST_asc = cumsum(diff)[1:K]

    # efdr over ascending cuts
    efdr_asc = Vector{Float64}(undef, K)
    @inbounds for i in 1:K
        num = E_asc[i] + EST_asc[i] + 2*ETS_asc[i]
        den = T_asc[i] + E_asc[i]
        efdr_asc[i] = den > 0 ? num/den : NaN
    end

    # Map back to per-row EFDR by the entrap score threshold index
    entrapment_fdr = zeros(Float64, length(e))
    @inbounds for i in eachindex(e)
        ie = idx_e_all[i]
        if ie isa Int && (ie::Int) > 0
            entrapment_fdr[i] = efdr_asc[ie::Int]
        else
            entrapment_fdr[i] = 0.0
        end
    end

    # Enforce monotonicity over the sorted (-score, qval) order
    sort_indices = sortperm(collect(zip(-method.score, method.qval)))
    
    fdr = Inf
    @inbounds @fastmath for i in reverse(sort_indices)
        v = entrapment_fdr[i]
        if !isfinite(v)
            entrapment_fdr[i] = fdr < Inf ? fdr : 0.0
        elseif v > fdr
            entrapment_fdr[i] = fdr
        else
            fdr = v
        end
    end
    
    return entrapment_fdr
end

