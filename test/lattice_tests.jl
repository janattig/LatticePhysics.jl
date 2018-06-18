# ALL TESTS THAT CONCERN UNITCELLS

# include
using LatticePhysics
using Base.Test

# begin the Unitcell testblock
lattice_testset = @testset "Lattice tests" begin


################################################################################
#
#   TYPE DEFINITIONS OF OBJECT CLASSES IN JULIA
#
################################################################################

# begin the testset
@testset "General Lattice functions" begin

    # get a test unitcell
    unitcell = getUnitcellHoneycomb()
    # get a test lattice
    lattice = getLattice(unitcell, -4)


    @testset "File access with JLD" begin
        # test saving
        @test typeof(saveLattice(lattice))==String
        # test loading
        @test typeof(loadLattice(lattice))==String
        @test typeof(loadLattice(lattice.filename))==Lattice
    end;


    @testset "Test and Information functions" begin
        # test printing information
        @test printInfo(lattice)==nothing
        @test printInfo(lattice, detailed=true)==nothing
        # test testing method
        #@test typeof(testUnitcell(unitcell, 2, 2))==Bool
    end;




    @testset "Connectivity information functions" begin
        # testing connectivity information
        @test typeof(getConnectivityList(lattice)) == Array{Array,1}
        # testing connection information
        @test typeof(getConnectionList(lattice)) == Array{Array,1}
        # testing connection strength information
        @test typeof(getConnectionStrengthList(lattice)) == Array{Any,1}
    end

# end the testset here
end;









# end the lattice test block here
end;
