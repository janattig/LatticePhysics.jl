using LatticePhysics
using Base.Test


# TEST EVERYTHING FROM HERE ON
latticephysics_testset = @testset "All LatticePhysics.jl tests" begin

    # UNITCELL STUFF
    include("unitcell_tests.jl")

    # LATTICE STUFF
    include("lattice_tests.jl")

    # PLAQUETTE STUFF
    include("plaquette_tests.jl")


    # SVG PLOTTING STUFF
    include("SVG_plotting_tests.jl")

end;
