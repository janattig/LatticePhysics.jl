using LatticePhysics
using Base.Test


# TEST EVERYTHING FROM HERE ON
latticephysics_testset = @testset "All LatticePhysics.jl tests" begin

    # UNITCELL STUFF
    include("unitcell_tests.jl")

    # LATTICE STUFF

end;
