using ChemicalIdentifiers
using Test

@testset "ChemicalIdentifiers.jl" begin
    res1 = search_chemical("water",nothing)
    res2 =  search_chemical("SMILES=O",nothing)
    @test res1.formula == res2.formula
end
