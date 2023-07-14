using ChemicalIdentifiers
using Test

@testset "ChemicalIdentifiers.jl" begin
    res1 = search_chemical("water",nothing)
    res2 =  search_chemical("SMILES=O",nothing)
    res3 = search_chemical("water (H2O)", nothing)
    @test res1.formula == res2.formula
    @test res1.formula == res3.formula
    @test ismissing(search_chemical("[3-(2,3-EPOXYPROXY)PROPYL]TRIMETHOXYSILANE",nothing))
end

@testset "issue #10" begin
    @test !ismissing(search_chemical("propyl    ethylether   ",nothing)) #spurious spaces
end