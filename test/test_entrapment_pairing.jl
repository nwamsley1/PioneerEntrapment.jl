using Test
using DataFrames
using PioneerEntrapment

@testset "getModKey tests" begin
    @test getModKey("(5,M,x)(5,M,Unimod:4)(5,M,Unimod:1)") == "Unimod:1;Unimod:4;x"
    @test getModKey("(5,M,Unimod:4)(5,M,Unimod:35)") == "Unimod:35;Unimod:4"
    @test getModKey("(5,M,Unimod:35)(5,M,Unimod:4)") == "Unimod:35;Unimod:4"
    @test getModKey("(5,M,x)") == "x"
end

@testset "assign_entrapment_pairs! tests" begin
    test_df1 = DataFrame(
        base_pep_id = [1, 1, 1],
        prec_charge = [2, 2, 2],
        is_decoy = [false, false, false],
        mod_key = ["", "", ""],
        entrapment_group_id = [0, 1, 1],
        sequence = ["PEPTIDE", "PEPTIDE", "PEPTIDE"]
    )
    assign_entrapment_pairs!(test_df1)
    @test !any(ismissing.(test_df1.entrap_pair_id[1:2]))
    @test test_df1.entrap_pair_id[1] == test_df1.entrap_pair_id[2]
    @test ismissing(test_df1.entrap_pair_id[3])

    test_df2 = DataFrame(
        base_pep_id = fill(1, 6),
        prec_charge = fill(2, 6),
        is_decoy = fill(false, 6),
        mod_key = fill("", 6),
        entrapment_group_id = [0, 0, 1, 1, 2, 2],
        sequence = fill("PEPTIDE", 6)
    )
    assign_entrapment_pairs!(test_df2)
    @test length(unique(test_df2.entrap_pair_id[.!ismissing.(test_df2.entrap_pair_id)])) == 2
    @test count(test_df2.entrap_pair_id .== test_df2.entrap_pair_id[1]) == 3
    @test count(test_df2.entrap_pair_id .== test_df2.entrap_pair_id[2]) == 3

    test_df3 = DataFrame(
        base_pep_id = [1, 1, 1],
        prec_charge = [2, 2, 2],
        is_decoy = [false, false, false],
        mod_key = ["", "", ""],
        entrapment_group_id = [1, 1, 2],
        sequence = ["PEPTIDE", "PEPTIDE", "PEPTIDE"]
    )
    assign_entrapment_pairs!(test_df3)
    @test all(ismissing.(test_df3.entrap_pair_id))

    test_df4 = DataFrame(
        base_pep_id = [1, 1, 2, 2],
        prec_charge = [2, 2, 2, 2],
        is_decoy = [false, false, false, false],
        mod_key = ["", "", "", ""],
        entrapment_group_id = [0, 1, 0, 1],
        sequence = ["PEPTIDE", "PEPTIDE", "PROTEIN", "PROTEIN"]
    )
    assign_entrapment_pairs!(test_df4)
    @test length(unique(skipmissing(test_df4.entrap_pair_id))) == 2
    @test !ismissing(test_df4.entrap_pair_id[1]) && !ismissing(test_df4.entrap_pair_id[2])
    @test !ismissing(test_df4.entrap_pair_id[3]) && !ismissing(test_df4.entrap_pair_id[4])
    @test test_df4.entrap_pair_id[1] == test_df4.entrap_pair_id[2]
    @test test_df4.entrap_pair_id[3] == test_df4.entrap_pair_id[4]
    @test test_df4.entrap_pair_id[1] != test_df4.entrap_pair_id[3]

    test_df5 = DataFrame(
        base_pep_id = fill(1, 5),
        prec_charge = fill(2, 5),
        is_decoy = fill(false, 5),
        mod_key = fill("", 5),
        entrapment_group_id = [0, 0, 0, 1, 1],
        sequence = fill("PEPTIDE", 5)
    )
    assign_entrapment_pairs!(test_df5)
    @test length(unique(skipmissing(test_df5.entrap_pair_id))) == 3
end

