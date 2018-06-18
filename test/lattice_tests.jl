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
        @test typeof(testLattice(lattice, 2, 2))==Bool
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





################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT LATTICE CONSTRUCTION RELATED STUFF
#   AS WELL AS DIFFERNT LATTICE MODIFICATION STUFF
#
################################################################################

# begin the testset
@testset "Lattice construction functions" begin

    # get a test unitcell
    unitcell_2d = getUnitcellHoneycomb()
    unitcell_3d = getUnitcellDiamond()


    # BASED ON UNITCELLS
    @testset "Construction based on Unitcells" begin

        # PERIODIC
        @testset "Periodic" begin
            # SPECIFIC FUNCTION
            # - just getting lattice
            @test testLattice(getLatticePeriodic(unitcell_2d, [-4,-4]), 2, 2)
            @test testLattice(getLatticePeriodic(unitcell_3d, [-4,-4,-4]), 3, 3)
            # - getting lattice and saving
            @test testLattice(getLatticePeriodic(unitcell_2d, [-4,-4], save=true), 2, 2)
            @test testLattice(getLatticePeriodic(unitcell_3d, [-4,-4,-4], save=true), 3, 3)
            # - getting lattice and saving and loading
            @test testLattice(getLatticePeriodic(unitcell_2d, [-4,-4], save=true, load=true), 2, 2)
            @test testLattice(getLatticePeriodic(unitcell_3d, [-4,-4,-4], save=true, load=true), 3, 3)
            # GENERAL FUNCTION
            # - just getting lattice
            @test testLattice(getLattice(unitcell_2d, [-4,-4]), 2, 2)
            @test testLattice(getLattice(unitcell_3d, [-4,-4,-4]), 3, 3)
            # - getting lattice and saving
            @test testLattice(getLattice(unitcell_2d, [-4,-4], save=true), 2, 2)
            @test testLattice(getLattice(unitcell_3d, [-4,-4,-4], save=true), 3, 3)
            # - getting lattice and saving and loading
            @test testLattice(getLattice(unitcell_2d, [-4,-4], save=true, load=true), 2, 2)
            @test testLattice(getLattice(unitcell_3d, [-4,-4,-4], save=true, load=true), 3, 3)
        end;

        # SEMI PERIODIC
        @testset "Semi-Periodic" begin
            # SPECIFIC FUNCTION
            # - just getting lattice
            @test testLattice(getLatticeSemiperiodic(unitcell_2d, [4,-4]), 1, 2)
            @test testLattice(getLatticeSemiperiodic(unitcell_3d, [4,-4,-4]), 2, 3)
            # - getting lattice and saving
            @test testLattice(getLatticeSemiperiodic(unitcell_2d, [4,-4], save=true), 1, 2)
            @test testLattice(getLatticeSemiperiodic(unitcell_3d, [4,-4,-4], save=true), 2, 3)
            # - getting lattice and saving and loading
            @test testLattice(getLatticeSemiperiodic(unitcell_2d, [4,-4], save=true, load=true), 1, 2)
            @test testLattice(getLatticeSemiperiodic(unitcell_3d, [4,-4,-4], save=true, load=true), 2, 3)
            # GENERAL FUNCTION
            # - just getting lattice
            @test testLattice(getLattice(unitcell_2d, [4,-4]), 1, 2)
            @test testLattice(getLattice(unitcell_3d, [4,-4,-4]), 2, 3)
            # - getting lattice and saving
            @test testLattice(getLattice(unitcell_2d, [4,-4], save=true), 1, 2)
            @test testLattice(getLattice(unitcell_3d, [4,-4,-4], save=true), 2, 3)
            # - getting lattice and saving and loading
            @test testLattice(getLattice(unitcell_2d, [4,-4], save=true, load=true), 1, 2)
            @test testLattice(getLattice(unitcell_3d, [4,-4,-4], save=true, load=true), 2, 3)
        end;

        # OPEN
        @testset "Open" begin
            # SPECIFIC FUNCTION
            # - just getting lattice
            @test testLattice(getLatticeOpen(unitcell_2d, [4,4]), 0, 2)
            @test testLattice(getLatticeOpen(unitcell_3d, [4,4,4]), 0, 3)
            # - getting lattice and saving
            @test testLattice(getLatticeOpen(unitcell_2d, [4,4], save=true), 0, 2)
            @test testLattice(getLatticeOpen(unitcell_3d, [4,4,4], save=true), 0, 3)
            # - getting lattice and saving and loading
            @test testLattice(getLatticeOpen(unitcell_2d, [4,4], save=true, load=true), 0, 2)
            @test testLattice(getLatticeOpen(unitcell_3d, [4,4,4], save=true, load=true), 0, 3)
            # GENERAL FUNCTION
            # - just getting lattice
            @test testLattice(getLattice(unitcell_2d, [4,4]), 0, 2)
            @test testLattice(getLattice(unitcell_3d, [4,4,4]), 0, 3)
            # - getting lattice and saving
            @test testLattice(getLattice(unitcell_2d, [4,4], save=true), 0, 2)
            @test testLattice(getLattice(unitcell_3d, [4,4,4], save=true), 0, 3)
            # - getting lattice and saving and loading
            @test testLattice(getLattice(unitcell_2d, [4,4], save=true, load=true), 0, 2)
            @test testLattice(getLattice(unitcell_3d, [4,4,4], save=true, load=true), 0, 3)
        end;

    # end of unitcell based construction
    end;



    # BASED ON BOND DISTANCE
    @testset "Construction based on Bond Distance" begin

        # - just getting lattice
        @test testLattice(getLatticeByBondDistance(unitcell_2d, 4), 0, 2)
        @test testLattice(getLatticeByBondDistance(unitcell_3d, 4), 0, 3)
        @test testLattice(getLatticeByBondDistance(unitcell_2d, 4, origin=2), 0, 2)
        @test testLattice(getLatticeByBondDistance(unitcell_3d, 4, origin=2), 0, 3)
        # - getting lattice and saving
        @test testLattice(getLatticeByBondDistance(unitcell_2d, 4, save=true), 0, 2)
        @test testLattice(getLatticeByBondDistance(unitcell_3d, 4, save=true), 0, 3)
        # - getting lattice and saving and loading
        @test testLattice(getLatticeByBondDistance(unitcell_2d, 4, save=true, load=true), 0, 2)
        @test testLattice(getLatticeByBondDistance(unitcell_3d, 4, save=true, load=true), 0, 3)

    # end of unitcell based construction
    end;


    # BASED ON SHAPE
    @testset "Construction based on Shape" begin

        # SPHERE / CIRCLE
        @testset "Sphere" begin
            # - just getting lattice
            @test testLattice(getLatticeInSphere(unitcell_2d, 2.0), 0, 2)
            @test testLattice(getLatticeInSphere(unitcell_3d, 2.0), 0, 3)
            @test testLattice(getLatticeInSphere(unitcell_2d, 2.0, origin=2), 0, 2)
            @test testLattice(getLatticeInSphere(unitcell_3d, 2.0, origin=2), 0, 3)
            # - getting lattice and saving
            @test testLattice(getLatticeInSphere(unitcell_2d, 2.0, save=true), 0, 2)
            @test testLattice(getLatticeInSphere(unitcell_3d, 2.0, save=true), 0, 3)
            # - getting lattice and saving and loading
            @test testLattice(getLatticeInSphere(unitcell_2d, 2.0, save=true, load=true), 0, 2)
            @test testLattice(getLatticeInSphere(unitcell_3d, 2.0, save=true, load=true), 0, 3)
        # end testset
        end;

        # BOX
        @testset "Box" begin
            # defining a box that workds for 2d and 3d
            box = [4.0, 2.0, 3.0]
            # - just getting lattice
            @test testLattice(getLatticeInBox(unitcell_2d, box), 0, 2)
            @test testLattice(getLatticeInBox(unitcell_3d, box), 0, 3)
            @test testLattice(getLatticeInBox(unitcell_2d, box, origin=2), 0, 2)
            @test testLattice(getLatticeInBox(unitcell_3d, box, origin=2), 0, 3)
            # - getting lattice and saving
            @test testLattice(getLatticeInBox(unitcell_2d, box, save=true), 0, 2)
            @test testLattice(getLatticeInBox(unitcell_3d, box, save=true), 0, 3)
            # - getting lattice and saving and loading
            @test testLattice(getLatticeInBox(unitcell_2d, box, save=true, load=true), 0, 2)
            @test testLattice(getLatticeInBox(unitcell_3d, box, save=true, load=true), 0, 3)
        # end testset
        end;

        # GENERAL
        @testset "General shape" begin
            # defining a box that workds for 2d and 3d
            shape_function(point) = sum(point.*point) < 2.0
            # - just getting lattice
            @test testLattice(getLatticeInShape(unitcell_2d, shape_function), 0, 2)
            @test testLattice(getLatticeInShape(unitcell_3d, shape_function), 0, 3)
            @test testLattice(getLatticeInShape(unitcell_2d, shape_function, "myshape"), 0, 2)
            @test testLattice(getLatticeInShape(unitcell_3d, shape_function, "myshape"), 0, 3)
            @test testLattice(getLatticeInShape(unitcell_2d, shape_function, "myshape", origin=2), 0, 2)
            @test testLattice(getLatticeInShape(unitcell_3d, shape_function, "myshape", origin=2), 0, 3)
            # - getting lattice and saving
            @test testLattice(getLatticeInShape(unitcell_2d, shape_function, "myshape", save=true), 0, 2)
            @test testLattice(getLatticeInShape(unitcell_3d, shape_function, "myshape", save=true), 0, 3)
            # - getting lattice and saving and loading
            @test testLattice(getLatticeInShape(unitcell_2d, shape_function, "myshape", save=true, load=true), 0, 2)
            @test testLattice(getLatticeInShape(unitcell_3d, shape_function, "myshape", save=true, load=true), 0, 3)
        # end testset
        end;

    # end of unitcell based construction
    end;



# end the testset here
end;











# end the lattice test block here
end;
