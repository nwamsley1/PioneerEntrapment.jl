@testset "Precursor EFDR comparison with no passing rows" begin
    df = DataFrame(
        precursor_idx=UInt32[1, 2], qval=Float32[0.005, 0.008],
        score=Float32[0.9, 0.8], score_combined_efdr=Float64[0, 0.5],
        score_paired_efdr=Float64[0, 0.25]
    )
    lib = DataFrame(entrapment_group_id=UInt8[0, 1])
    for (paired, has_paired) in ((true, true), (false, true), (true, false)),
        labels in (nothing, [0, 1])
        report_df = has_paired ? df : select(df, Not(:score_paired_efdr))
        result = PioneerEntrapment.compare_efdr_methods(report_df, :qval, :score, lib;
            include_paired=paired, entrap_labels_override=labels)
        @test result.threshold == [0.001, 0.01, 0.05, 0.1]
        @test result.qval_n == [0, 2, 2, 2]
        @test result.qval_actual_fdr == [0.0, 0.5, 0.5, 0.5]
        @test result.combined_n == [0, 2, 2, 2]
        @test result.combined_efdr == [0.0, 0.5, 0.5, 0.5]
        @test hasproperty(result, :paired_n) == (paired && has_paired)
        @test hasproperty(result, :paired_efdr) == (paired && has_paired)
        if paired && has_paired
            @test result.paired_n == [0, 2, 2, 2]
            @test result.paired_efdr == [0.0, 0.25, 0.25, 0.25]
        end

        empty_result = PioneerEntrapment.compare_efdr_methods(report_df[1:0, :], :qval, :score, lib;
            include_paired=paired, entrap_labels_override=isnothing(labels) ? nothing : Int[])
        @test empty_result.threshold == result.threshold
        @test names(empty_result) == names(result)
        @test all(iszero, empty_result.qval_n)
        @test all(iszero, empty_result.qval_actual_fdr)
        @test all(iszero, empty_result.combined_n)
        @test all(iszero, empty_result.combined_efdr)
        if paired && has_paired
            @test all(iszero, empty_result.paired_n)
            @test all(iszero, empty_result.paired_efdr)
        end
    end

    @testset "Partial selections and label overrides" begin
        # Unsorted rows ensure overrides follow the q-value sort and selection.
        unsorted_df = df[[2, 1], :]
        thresholds = [0.001, Float64(df.qval[1]), 0.01]
        for paired in (true, false),
            (labels, expected_fdr) in ((nothing, [0.0, 0.0, 0.5]),
                                      ([0, 2], [0.0, 1.0, 0.5]),
                                      ([2, 3], [0.0, 1.0, 1.0]))
            result = PioneerEntrapment.compare_efdr_methods(unsorted_df, :qval, :score, lib;
                thresholds=thresholds, include_paired=paired, entrap_labels_override=labels)
            @test result.threshold == thresholds
            @test result.qval_n == [0, 1, 2]
            @test result.qval_actual_fdr == expected_fdr
            @test result.combined_n == [0, 1, 2]
            @test result.combined_efdr == [0.0, 0.0, 0.5]
            if paired
                @test result.paired_n == [0, 1, 2]
                @test result.paired_efdr == [0.0, 0.0, 0.25]
            end
        end
    end
end
