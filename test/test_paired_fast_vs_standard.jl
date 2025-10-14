using Test
using PioneerEntrapment
using Random

# Import the fast implementation (not exported)
import PioneerEntrapment: calculate_efdr_fast

"""
Comprehensive test suite comparing the fast O(n log n) paired EFDR algorithm
with the standard O(n²) implementation (stride=1) to ensure identical results.
"""

@testset "Paired EFDR: Fast vs Standard Implementation" begin

    @testset "Small synthetic dataset" begin
        # Test with small dataset covering various edge cases
        e = Float64[0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50]
        t = Float64[0.90, 0.92, 0.80, 0.82, 0.70, 0.72, 0.50, 0.58, 0.45, 0.40]
        labs = Int[0, 1, 0, 1, 0, 1, 1, 0, 1, 0]
        q = Float64[0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]

        method = PairedEFDR(e, t, labs, q, 1.0)
        standard = calculate_efdr(method; stride=1)  # Exact quadratic
        fast = calculate_efdr_fast(method; cuts_mode=:all)

        @test length(standard) == length(fast)
        @test all(isfinite.(fast))
        @test maximum(abs.(standard .- fast)) < 1e-12
    end

    @testset "Larger dataset (n=100)" begin
        Random.seed!(42)
        n = 100

        # Generate realistic scores
        e = sort(rand(n), rev=true)  # Descending scores

        # Target scores: some lower, some higher than entrap
        t = similar(e)
        for i in 1:n
            if rand() < 0.5  # 50% lower than entrap
                t[i] = e[i] - rand() * 0.2
            else  # 50% higher than entrap
                t[i] = e[i] + rand() * 0.1
            end
        end

        labs = rand([0, 1], n)  # Random mix of targets and entrapments
        q = sort(rand(n))  # Ascending q-values

        method = PairedEFDR(e, t, labs, q, 1.0)
        standard = calculate_efdr(method; stride=1)
        fast = calculate_efdr_fast(method; cuts_mode=:all)

        @test length(standard) == length(fast)
        @test maximum(abs.(standard .- fast)) < 1e-10
        println("  Max difference (n=100): ", maximum(abs.(standard .- fast)))
    end

    @testset "All targets (no entrapments)" begin
        n = 50
        e = sort(rand(n), rev=true)
        t = e .- rand(n) .* 0.1
        labs = zeros(Int, n)  # All targets
        q = sort(rand(n))

        method = PairedEFDR(e, t, labs, q, 1.0)
        standard = calculate_efdr(method; stride=1)
        fast = calculate_efdr_fast(method; cuts_mode=:all)

        @test all(standard .≈ fast)
    end

    @testset "All entrapments (no targets)" begin
        n = 50
        e = sort(rand(n), rev=true)
        t = e .- rand(n) .* 0.1  # All slightly lower
        labs = ones(Int, n)  # All entrapments
        q = sort(rand(n))

        method = PairedEFDR(e, t, labs, q, 1.0)
        standard = calculate_efdr(method; stride=1)
        fast = calculate_efdr_fast(method; cuts_mode=:all)

        @test maximum(abs.(standard .- fast)) < 1e-10
    end

    @testset "Entrapment scores always > target scores" begin
        n = 50
        Random.seed!(123)
        t = sort(rand(n), rev=true)
        e = t .+ rand(n) .* 0.2  # Always higher
        labs = rand([0, 1], n)
        q = sort(rand(n))

        method = PairedEFDR(e, t, labs, q, 1.0)
        standard = calculate_efdr(method; stride=1)
        fast = calculate_efdr_fast(method; cuts_mode=:all)

        @test maximum(abs.(standard .- fast)) < 1e-10
    end

    @testset "Entrapment scores always < target scores" begin
        n = 50
        Random.seed!(456)
        t = sort(rand(n), rev=true)
        e = t .- rand(n) .* 0.2  # Always lower
        labs = rand([0, 1], n)
        q = sort(rand(n))

        method = PairedEFDR(e, t, labs, q, 1.0)
        standard = calculate_efdr(method; stride=1)
        fast = calculate_efdr_fast(method; cuts_mode=:all)

        @test maximum(abs.(standard .- fast)) < 1e-10
    end

    @testset "Identical scores (ties)" begin
        n = 30
        e = repeat([0.9, 0.7, 0.5], 10)
        t = repeat([0.85, 0.65, 0.45], 10)
        labs = repeat([0, 1, 1], 10)
        q = collect(range(0.01, 0.3, length=n))

        method = PairedEFDR(e, t, labs, q, 1.0)
        standard = calculate_efdr(method; stride=1)
        fast = calculate_efdr_fast(method; cuts_mode=:all)

        @test all(standard .≈ fast)
    end

    @testset "Order independence (shuffled data)" begin
        Random.seed!(789)
        n = 80
        e = sort(rand(n), rev=true)
        t = e .+ randn(n) .* 0.1
        labs = rand([0, 1], n)
        q = sort(rand(n))

        # Test with original order
        method1 = PairedEFDR(e, t, labs, q, 1.0)
        standard1 = calculate_efdr(method1; stride=1)
        fast1 = calculate_efdr_fast(method1; cuts_mode=:all)
        @test maximum(abs.(standard1 .- fast1)) < 1e-10

        # Test with shuffled order (multiple times)
        for trial in 1:5
            perm = randperm(n)
            e2 = e[perm]
            t2 = t[perm]
            labs2 = labs[perm]
            q2 = q[perm]

            method2 = PairedEFDR(e2, t2, labs2, q2, 1.0)
            standard2 = calculate_efdr(method2; stride=1)
            fast2 = calculate_efdr_fast(method2; cuts_mode=:all)

            @test maximum(abs.(standard2 .- fast2)) < 1e-10
        end
    end

    @testset "Different r values (target:entrapment ratio)" begin
        Random.seed!(999)
        n = 60
        e = sort(rand(n), rev=true)
        t = e .+ randn(n) .* 0.05
        labs = rand([0, 1], n)
        q = sort(rand(n))

        for r_val in [0.5, 1.0, 2.0]
            method = PairedEFDR(e, t, labs, q, r_val)
            standard = calculate_efdr(method; stride=1)
            fast = calculate_efdr_fast(method; cuts_mode=:all)

            @test maximum(abs.(standard .- fast)) < 1e-10
        end
    end

    @testset "Edge case: single element" begin
        e = Float64[0.9]
        t = Float64[0.8]
        labs = [0]
        q = Float64[0.01]

        method = PairedEFDR(e, t, labs, q, 1.0)
        standard = calculate_efdr(method; stride=1)
        fast = calculate_efdr_fast(method; cuts_mode=:all)

        @test standard == fast
    end

    @testset "Edge case: two elements" begin
        e = Float64[0.9, 0.8]
        t = Float64[0.85, 0.75]
        labs = [0, 1]
        q = Float64[0.01, 0.02]

        method = PairedEFDR(e, t, labs, q, 1.0)
        standard = calculate_efdr(method; stride=1)
        fast = calculate_efdr_fast(method; cuts_mode=:all)

        @test all(standard .≈ fast)
    end

    @testset "Stress test: Large dataset (n=500)" begin
        Random.seed!(2025)
        n = 500

        # Realistic scenario
        e = sort(rand(n), rev=true)
        t = e .+ randn(n) .* 0.08
        labs = rand([0, 1], n)
        q = sort(rand(n))

        method = PairedEFDR(e, t, labs, q, 1.0)

        # This might take a moment for the standard O(n²) method
        println("  Running standard method (n=500, O(n²))...")
        standard = calculate_efdr(method; stride=1)

        println("  Running fast method (n=500, O(n log n))...")
        fast = calculate_efdr_fast(method; cuts_mode=:all)

        max_diff = maximum(abs.(standard .- fast))
        @test max_diff < 1e-9
        println("  Max difference (n=500): ", max_diff)

        # Verify monotonicity is preserved
        sort_indices = sortperm(collect(zip(-e, q)))
        @test issorted(standard[sort_indices], rev=true) || all(standard .== 0.0)
        @test issorted(fast[sort_indices], rev=true) || all(fast .== 0.0)
    end
end
