using Test
using PioneerEntrapment

@testset "Fast Paired EFDR matches quadratic" begin
    # Construct a deterministic dataset covering cases:
    # e >= s and t < s; e >= s and t >= s with e > t; mixed orders, and missing targets
    e = Float64[0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50]
    t = Any[  0.90, 0.92, 0.80, 0.82, 0.70, 0.72, missing, 0.58, missing, 0.40]
    t = Float64[ismissing(x) ? NaN : x for x in t] # preserve as Float64; NaN will be treated as missing upstream
    # Labels: 0=target row (original), 1=entrap row
    labs = Int[0, 1, 0, 1, 0, 1, 1, 0, 1, 0]
    q = Float64[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10]

    # Replace NaN back to missing for method
    t_m = Vector{Union{Float64,Missing}}(undef, length(t))
    for i in eachindex(t)
        t_m[i] = isnan(t[i]) ? missing : t[i]
    end

    quad = calculate_efdr(PairedEFDR(e, t_m, labs, q, 1.0))
    fast = calculate_efdr_fast(PairedEFDR(e, t_m, labs, q, 1.0))
    @test length(quad) == length(fast)
    @test all(isfinite.(fast))
    @test maximum(abs.(quad .- fast)) < 1e-12

    # Shuffle multiple times to ensure order-independence
    for _ in 1:10
        perm = randperm(length(e))
        e2 = e[perm]; t2 = t_m[perm]; l2 = labs[perm]; q2 = q[perm]
        quad2 = calculate_efdr(PairedEFDR(e2, t2, l2, q2, 1.0))
        fast2 = calculate_efdr_fast(PairedEFDR(e2, t2, l2, q2, 1.0))
        @test maximum(abs.(quad2 .- fast2)) < 1e-12
    end
end

