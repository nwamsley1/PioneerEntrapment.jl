using Test
using PioneerEntrapment
using DataFrames
using Arrow

@testset "PioneerEntrapment" begin
    include("test_entrapment_pairing.jl")

    @testset "EFDR Method Tests" begin
        scores = [0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45]
        original_target_scores = [0.9, 0.9, 0.8, 0.8, 0.7, 0.7, 0.6, 0.6, 0.5, 0.5]
        entrapment_labels = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
        qvals = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
        combined_method = CombinedEFDR(scores, original_target_scores, entrapment_labels, qvals, 1.0)
        combined_efdr = calculate_efdr(combined_method)
        @test length(combined_efdr) == length(scores)
        @test all(0 .<= combined_efdr .<= 1)
        @test issorted(combined_efdr)
        paired_method = PairedEFDR(scores, original_target_scores, entrapment_labels, qvals, 1.0)
        paired_efdr = calculate_efdr(paired_method)
        @test length(paired_efdr) == length(scores)
        @test all(0 .<= paired_efdr .<= 1)
    end

    @testset "Scoring Function Tests" begin
        test_df = DataFrame(
            entrap_pair_id = [1, 1, 2, 2, missing],
            score = [0.9, 0.8, 0.7, 0.6, 0.5],
            precursor_idx = [1, 2, 3, 4, 5]
        )
        library_df = DataFrame(entrapment_group_id = [0, 1, 0, 1, 0])
        @test get_complement_score(test_df, 1; score_col=:score) == 0.8
        @test get_complement_score(test_df, 2; score_col=:score) == 0.9
        @test get_complement_score(test_df, 5; score_col=:score) == 0.0
        add_original_target_scores!(test_df, library_df; score_col=:score)
        @test hasproperty(test_df, :score_original_target)
        @test test_df.score_original_target[1] == 0.9
        @test test_df.score_original_target[2] == 0.9
    end

    include("test_protein_entrapment.jl")
end

println("\nAll tests completed!")
