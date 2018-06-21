# ALL TESTS THAT CONCERN MODFICIATION

# include
using LatticePhysics
using Base.Test

# begin the Modification testblock
modification_testset = @testset "Modification tests" begin


################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT LATTICE MODIFICATION FUNCTIONS
#   (MOST OF THEM ALSO WORK FOR UNITCELLS)
#
################################################################################

################################################################################
#
#   1) CONNECTIONS ADDING REMOVING
#       - Add a new connection
#       - remove connections (based on indices)
#       - remove connections (based on strength)
#
################################################################################

# begin the testset
@testset "Connection adding & removing" begin

    # First, Unitcell
    @testset "Adding and removing in Unitcells (honeycomb)" begin

        # get a test unitcell
        unitcell = getUnitcellHoneycomb(4)

        # First, try to add connections
        @testset "Adding connections" begin
            # try adding a random connection and testing the unitcell afterwards
            @test addConnection!(unitcell, 1, 1, "t1", (1,1)) == nothing
            @test testUnitcell(unitcell, 2,2)
            # test whether there are now 8 connections
            @test length(unitcell.connections) == 8
            # test whether there are now 4 connection strengths
            @test length(getConnectionStrengthList(unitcell)) == 4
            # test if one adds up the same connection again, the numbers stay the same
            @test addConnection!(unitcell, 1, 1, "t1", (1,1)) == nothing
            @test testUnitcell(unitcell, 2,2)
            @test length(unitcell.connections) == 8
            @test length(getConnectionStrengthList(unitcell)) == 4
            # test if one adds up the another connection but with different wrap, the numbers change
            @test addConnection!(unitcell, 1, 1, "t1", (1,-1)) == nothing
            @test testUnitcell(unitcell, 2,2)
            @test length(unitcell.connections) == 10
            @test length(getConnectionStrengthList(unitcell)) == 4
        end;

        # get a fresh test unitcell
        unitcell = getUnitcellHoneycomb(4)

        # Then, try to remove connections (based on index&strength)
        @testset "Removing connections based on index & strength" begin
            # try adding a random connection and testing the unitcell afterwards
            @test removeConnections!(unitcell, 1, 2, "tx") == nothing
            @test testUnitcell(unitcell, 2,2)
            # test whether there are now 4 connections
            @test length(unitcell.connections) == 4
            # test whether there are now 2 connection strengths
            @test length(getConnectionStrengthList(unitcell)) == 2
            # test if one adds up the same connection again, the numbers stay the same
            @test removeConnections!(unitcell, 1, 2, "tx") == nothing
            @test testUnitcell(unitcell, 2,2)
            @test length(unitcell.connections) == 4
            @test length(getConnectionStrengthList(unitcell)) == 2
        end;

        # get a fresh test unitcell
        unitcell = getUnitcellSquareOctagon(4)

        # Then, try to remove connections (based on index)
        @testset "Removing connections based on index" begin
            # try adding a random connection and testing the unitcell afterwards
            @test removeConnections!(unitcell, 1, 2) == nothing
            @test testUnitcell(unitcell, 2,2)
            # test whether there are now 10 connections
            @test length(unitcell.connections) == 10
            # test whether there are still 3 connection strengths
            @test length(getConnectionStrengthList(unitcell)) == 3
            # test if one adds up the same connection again, the numbers stay the same
            @test removeConnections!(unitcell, 1, 2) == nothing
            @test testUnitcell(unitcell, 2,2)
            @test length(unitcell.connections) == 10
            @test length(getConnectionStrengthList(unitcell)) == 3
        end;

        # get a fresh test unitcell
        unitcell = getUnitcellSquareOctagon(4)

        # Then, try to remove connections (based on strength)
        @testset "Removing connections based on index" begin
            # try adding a random connection and testing the unitcell afterwards
            @test removeConnections!(unitcell, "tx") == nothing
            @test testUnitcell(unitcell, 2,2)
            # test whether there are now 8 connections
            @test length(unitcell.connections) == 8
            # test whether there are still 2 connection strengths
            @test length(getConnectionStrengthList(unitcell)) == 2
            # test if one adds up the same connection again, the numbers stay the same
            @test removeConnections!(unitcell, 1, 2) == nothing
            @test testUnitcell(unitcell, 2,2)
            @test length(unitcell.connections) == 8
            @test length(getConnectionStrengthList(unitcell)) == 2
        end;

    # end the Unitcell testset here
    end;


    # Then, Lattice
    @testset "Adding and removing in Lattices (honeycomb 4x2 periodic)" begin

        # get a test unitcell
        unitcell = getUnitcellHoneycomb(4)
        # get a test lattice
        lattice = getLattice(unitcell, [-4, -2])


        # First, try to add connections
        @testset "Adding connections" begin
            # try adding a random connection and testing the unitcell afterwards
            @test addConnection!(lattice, 2, 7, "t1", (1,1)) == nothing
            @test testLattice(lattice, 2,2)
            # test whether there are now 50 connections
            @test length(lattice.connections) == 50
            # test whether there are now 4 connection strengths
            @test length(getConnectionStrengthList(lattice)) == 4
            # test if one adds up the same connection again, the numbers stay the same
            @test addConnection!(lattice, 2, 7, "t1", (1,1)) == nothing
            @test testLattice(lattice, 2,2)
            @test length(lattice.connections) == 50
            @test length(getConnectionStrengthList(lattice)) == 4
            # test if one adds up the another connection but with different wrap, the numbers change
            @test addConnection!(lattice, 1, 1, "t1", (1,-1)) == nothing
            @test testLattice(lattice, 2,2)
            @test length(lattice.connections) == 52
            @test length(getConnectionStrengthList(lattice)) == 4
        end;

        # get a test unitcell
        unitcell = getUnitcellHoneycomb(4)
        # get a test lattice
        lattice = getLattice(unitcell, [-4, -2])

        # Then, try to remove connections (based on index&strength)
        @testset "Removing connections based on index & strength" begin
            # try adding a random connection and testing the unitcell afterwards
            @test removeConnections!(lattice, 1, 2, "tx") == nothing
            @test testLattice(lattice, 2,2)
            # test whether there are now 46 connections
            @test length(lattice.connections) == 46
            # test whether there are now 3 connection strengths
            @test length(getConnectionStrengthList(lattice)) == 3
            # test if one adds up the same connection again, the numbers stay the same
            @test removeConnections!(lattice, 1, 2, "tx") == nothing
            @test testLattice(lattice, 2,2)
            @test length(lattice.connections) == 46
            @test length(getConnectionStrengthList(lattice)) == 3
        end;

        # get a test unitcell
        unitcell = getUnitcellHoneycomb(4)
        # get a test lattice
        lattice = getLattice(unitcell, [-4, -2])

        # Then, try to remove connections (based on index)
        @testset "Removing connections based on index" begin
            # try adding a random connection and testing the unitcell afterwards
            @test removeConnections!(lattice, 6, 7) == nothing
            @test testLattice(lattice, 2,2)
            # test whether there are now 46 connections
            @test length(lattice.connections) == 46
            # test whether there are still 3 connection strengths
            @test length(getConnectionStrengthList(lattice)) == 3
            # test if one adds up the same connection again, the numbers stay the same
            @test removeConnections!(lattice, 6, 7) == nothing
            @test testLattice(lattice, 2,2)
            @test length(lattice.connections) == 46
            @test length(getConnectionStrengthList(lattice)) == 3
        end;

        # get a test unitcell
        unitcell = getUnitcellHoneycomb(4)
        # get a test lattice
        lattice = getLattice(unitcell, [-4, -2])

        # Then, try to remove connections (based on strength)
        @testset "Removing connections based on index" begin
            # try adding a random connection and testing the unitcell afterwards
            @test removeConnections!(lattice, "tx") == nothing
            @test testLattice(lattice, 2,2)
            # test whether there are now 32 connections
            @test length(lattice.connections) == 32
            # test whether there are still 2 connection strengths
            @test length(getConnectionStrengthList(lattice)) == 2
            # test if one adds up the same connection again, the numbers stay the same
            @test removeConnections!(lattice, 1, 2) == nothing
            @test testLattice(lattice, 2,2)
            @test length(lattice.connections) == 32
            @test length(getConnectionStrengthList(lattice)) == 2
        end;


    # end the Unitcell testset here
    end;

# end the testset here
end;





################################################################################
#
#   2) SITES ADDING REMOVING
#       - add a new site
#       - remove a site (based on index)
#
################################################################################

################################################################################
#
#   3) CONNECTION STRENGTH MODIFICATION
#       - all connections
#       - mapping of strengths
#       - evaluate strengths
#       - optimize connections
#
################################################################################



# end the Modification test block here
end;
