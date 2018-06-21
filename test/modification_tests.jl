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
@testset "Adding & removing connections" begin

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

    # end the Lattice testset here
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

# begin the testset
@testset "Adding & removing sites" begin

    # First, Unitcell
    @testset "Adding and removing in Unitcells (honeycomb)" begin

        # get a test unitcell
        unitcell = getUnitcellHoneycomb(4)

        # First, try to add site
        @testset "Adding sites" begin
            # try adding a random connection and testing the unitcell afterwards
            @test typeof(addSite!(unitcell, [1.0, 2.0])) == Int64
            @test testUnitcell(unitcell, 2,2)
        end;
        # Then, try to remove site
        @testset "Removing sites" begin
            # try removing the last (newly added) site and testing the unitcell afterwards
            @test removeSite!(unitcell, 3) == nothing
            @test testUnitcell(unitcell, 2,2)
            # try removing the first site and testing the unitcell afterwards
            @test removeSite!(unitcell, 1) == nothing
            @test testUnitcell(unitcell, 2,2)
        end;

    # end the Unitcell testset here
    end;

    # Then, Lattice
    @testset "Adding and removing in Lattices (honeycomb 4x2 periodic)" begin

        # get a test unitcell
        unitcell = getUnitcellHoneycomb(4)
        # get a test lattice
        lattice = getLattice(unitcell, [-4, -2])

        # First, try to add site
        @testset "Adding sites" begin
            # try adding a random connection and testing the unitcell afterwards
            @test typeof(addSite!(lattice, [1.0, 2.0])) == Int64
            @test testLattice(lattice, 2,2)
        end;
        # Then, try to remove site
        @testset "Removing sites" begin
            # try removing the last (newly added) site and testing the unitcell afterwards
            @test removeSite!(lattice, 17) == nothing
            @test testLattice(lattice, 2,2)
            # try removing some intermediate site and testing the unitcell afterwards
            @test removeSite!(lattice, 5) == nothing
            @test testLattice(lattice, 2,2)
        end;

    # end the Lattice testset here
    end;

# end the testset here
end;





################################################################################
#
#   3) CONNECTION STRENGTH MODIFICATION
#       - all connections
#       - mapping of strengths
#       - evaluate strengths
#       - optimize connections
#
################################################################################

# begin the testset
@testset "Modifying connection strengths" begin

    # setting all connection strengths
    @testset "Setting all connection strengths" begin

        # get a test unitcell
        unitcell = getUnitcellHoneycomb(4)
        # then, get a test lattice
        lattice  = getLattice(unitcell, [-4,-2])

        # First, Unitcell
        @testset "Setting in Unitcell" begin
            @test setAllInteractionStrengths!(unitcell, 1.0) == nothing
            @test testUnitcell(unitcell, 2,2)
            @test length(getConnectionStrengthList(unitcell)) == 1
        end;
        # Then, Lattice
        @testset "Setting in Lattice" begin
            @test setAllInteractionStrengths!(lattice, 1.0) == nothing
            @test testLattice(lattice, 2,2)
            @test length(getConnectionStrengthList(lattice)) == 1
        end;

    # end the testset here
    end;

    # Only mapping desired connection strengths
    @testset "Mapping connection strengths" begin

        # get a test unitcell
        unitcell = getUnitcellHoneycomb(4)
        # then, get a test lattice
        lattice  = getLattice(unitcell, [-4,-2])

        # First, Unitcell
        @testset "Mapping in Unitcell" begin
            # first, map only a part
            mapping = Dict("tx"=>"t", "ty"=>"t")
            @test mapInteractionStrengths!(unitcell, mapping) == nothing
            @test testUnitcell(unitcell, 2,2)
            @test length(getConnectionStrengthList(unitcell)) == 2
            # then, map inside strings
            mapping = Dict("t"=>"J")
            @test mapInteractionStrengths!(unitcell, mapping, replace_in_strings=true) == nothing
            @test testUnitcell(unitcell, 2,2)
            @test length(getConnectionStrengthList(unitcell)) == 2
            @test getConnectionStrengthList(unitcell)[1][1] == 'J'
            # then, map and evaluate
            mapping = Dict("J"=>1.0, "Jz"=>2.0)
            @test mapInteractionStrengths!(unitcell, mapping, replace_in_strings=false, evaluate=true) == nothing
            @test testUnitcell(unitcell, 2,2)
            @test length(getConnectionStrengthList(unitcell)) == 2
            @test typeof(getConnectionStrengthList(unitcell)[1]) == Float64
        end;
        # Then, Lattice
        @testset "Mapping in Lattice" begin
            # first, map only a part
            mapping = Dict("tx"=>"t", "ty"=>"t")
            @test mapInteractionStrengths!(lattice, mapping) == nothing
            @test testLattice(lattice, 2,2)
            @test length(getConnectionStrengthList(lattice)) == 2
            # then, map inside strings
            mapping = Dict("t"=>"J")
            @test mapInteractionStrengths!(lattice, mapping, replace_in_strings=true) == nothing
            @test testLattice(lattice, 2,2)
            @test length(getConnectionStrengthList(lattice)) == 2
            @test getConnectionStrengthList(lattice)[1][1] == 'J'
            # then, map and evaluate
            mapping = Dict("J"=>1.0, "Jz"=>2.0)
            @test mapInteractionStrengths!(lattice, mapping, replace_in_strings=false, evaluate=true) == nothing
            @test testLattice(lattice, 2,2)
            @test length(getConnectionStrengthList(lattice)) == 2
            @test typeof(getConnectionStrengthList(lattice)[1]) == Float64
        end;

    # end the testset here
    end;

    # Only mapping desired connection strengths
    @testset "Evaluating connection strengths" begin

        # get a test unitcell
        unitcell = getUnitcellHoneycomb(4)
        # then, get a test lattice
        lattice  = getLattice(unitcell, [-4,-2])

        # First, Unitcell
        @testset "Evaluating in Unitcell" begin
            # evaluate each kitaev parameter
            parameters = Dict("tx"=>1.0, "ty"=>1.0, "tz"=>2.0)
            @test evaluateInteractionStrengths!(unitcell, parameters) == nothing
            @test testUnitcell(unitcell, 2,2)
            @test length(getConnectionStrengthList(unitcell)) == 2
            @test typeof(getConnectionStrengthList(unitcell)[1]) == Float64
        end;
        # Then, Lattice
        @testset "Evaluating in Lattice" begin
            # evaluate each kitaev parameter
            parameters = Dict("tx"=>1.0, "ty"=>1.0, "tz"=>2.0)
            @test evaluateInteractionStrengths!(lattice, parameters) == nothing
            @test testLattice(lattice, 2,2)
            @test length(getConnectionStrengthList(lattice)) == 2
            @test typeof(getConnectionStrengthList(lattice)[1]) == Float64
        end;

    # end the Unitcell testset here
    end;

    # Optimizing all connection, first Float64
    @testset "Optimizing connections (Float64)" begin

        # get a test unitcell
        unitcell = getUnitcellHoneycomb(1)
        # then, get a test lattice
        lattice  = getLattice(unitcell, [-4,-2])

        # Optimizing, Unitcell
        @testset "Optimizing in Unitcell" begin
            # add a redundant connection
            @test addConnection!(unitcell, 1,2, 2.0, (0,-1), overwrite=true) == nothing
            @test testUnitcell(unitcell, 2,2)
            @test length(unitcell.connections) == 8
            # optimize
            @test optimizeConnections!(unitcell) == nothing
            @test testUnitcell(unitcell, 2,2)
            @test length(getConnectionStrengthList(unitcell)) == 2
            @test length(unitcell.connections) == 6
            @test typeof(getConnectionStrengthList(unitcell)[1]) == Float64
        end;
        # Optimizing, Lattice
        @testset "Optimizing in Lattice" begin
            # add a redundant connection
            @test addConnection!(lattice, 1,2, 2.0, (0,0), overwrite=true) == nothing
            @test testLattice(lattice, 2,2)
            @test length(lattice.connections) == 50
            # optimize
            @test optimizeConnections!(lattice) == nothing
            @test testLattice(lattice, 2,2)
            @test length(getConnectionStrengthList(lattice)) == 2
            @test length(lattice.connections) == 48
            @test typeof(getConnectionStrengthList(lattice)[1]) == Float64
        end;

    # end the Unitcell testset here
    end;

    # Optimizing all connection, second String
    @testset "Optimizing connections (String)" begin

        # get a test unitcell
        unitcell = getUnitcellHoneycomb(4)
        # then, get a test lattice
        lattice  = getLattice(unitcell, [-4,-2])

        # Optimizing, Unitcell
        @testset "Optimizing in Unitcell" begin
            # add a redundant connection
            @test addConnection!(unitcell, 1,2, "t", (0,-1), overwrite=true) == nothing
            @test testUnitcell(unitcell, 2,2)
            @test length(unitcell.connections) == 8
            # optimize
            @test optimizeConnections!(unitcell) == nothing
            @test testUnitcell(unitcell, 2,2)
            @test length(getConnectionStrengthList(unitcell)) == 3
            @test length(unitcell.connections) == 6
            @test typeof(getConnectionStrengthList(unitcell)[1]) == String
        end;
        # Optimizing, Lattice
        @testset "Optimizing in Lattice" begin
            # add a redundant connection
            @test addConnection!(lattice, 1,2, "t", (0,0), overwrite=true) == nothing
            @test testLattice(lattice, 2,2)
            @test length(lattice.connections) == 50
            # optimize
            @test optimizeConnections!(lattice) == nothing
            @test testLattice(lattice, 2,2)
            @test length(getConnectionStrengthList(lattice)) == 4
            @test length(lattice.connections) == 48
            @test typeof(getConnectionStrengthList(lattice)[1]) == String
        end;

    # end the Unitcell testset here
    end;

# end the testset here
end;






# end the Modification test block here
end;
