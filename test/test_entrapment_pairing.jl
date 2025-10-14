using Test
using DataFrames
using PioneerEntrapment

@testset "getModKey tests" begin
    @test getModKey("(5,M,x)(5,M,Unimod:4)(5,M,Unimod:1)") == "Unimod:1;Unimod:4;x"
    @test getModKey("(5,M,Unimod:4)(5,M,Unimod:35)") == "Unimod:35;Unimod:4"
    @test getModKey("(5,M,Unimod:35)(5,M,Unimod:4)") == "Unimod:35;Unimod:4"
    @test getModKey("(5,M,x)") == "x"
end

@testset "add_entrap_pair_ids! mapping" begin
    library_precursors = DataFrame(
        entrapment_pair_id = Union{Missing,UInt32}[1, 1, missing, 2],
    )
    prec_results = DataFrame(
        precursor_idx = [1, 2, 3, 4]
    )
    add_entrap_pair_ids!(prec_results, library_precursors)
    @test hasproperty(prec_results, :entrapment_pair_id)
    @test all(prec_results.entrapment_pair_id .=== Union{Missing,UInt32}[1, 1, missing, 2])
end
