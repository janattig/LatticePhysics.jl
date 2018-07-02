# ALL TESTS THAT CONCERN MODFICIATION

# include
using LatticePhysics
using PyPlot
using Base.Test

# begin the Modification testblock
fermi_surface_testset = @testset "Fermi surfaces" begin

################################################################################
#
#   METHODS FOR CONSTRUCTION OF BANDSTRUCTURES (ALONG PATHS)
#
################################################################################

################################################################################
#
#   1) TYPE BANDSTRUCTURE
#       - type definition
#       - TODO printInfo function
#
################################################################################

################################################################################
#
#   2) CALCULATION OF BAND STRUCTURES OF UNTICELL OBJECTS
#
################################################################################


# begin the testset
@testset "2D" begin

    # begin the testset
    @testset "honeycomb" begin

        # get a unitcell
        unitcell = getUnitcellHoneycomb()
        # get a path
        path = getDefaultPathTriangular(resolution=30)

        # test the calculations (bare)
        @test typeof(getBandStructureAlongPath(unitcell, path)) == Bandstructure

        # test the calculations (resolution set)
        @test typeof(getBandStructureAlongPath(unitcell, path, resolution=60)) == Bandstructure

        # test the calculations (hermitian)
        @test typeof(getBandStructureAlongPath(unitcell, path, enforce_hermitian=true)) == Bandstructure

    # end the testset here
    end;


    # begin the testset
    @testset "square" begin

        # get a unitcell
        unitcell = getUnitcellSquare()
        # get a path
        path = getDefaultPathSquare(resolution=30)

        # test the calculations (bare)
        @test typeof(getBandStructureAlongPath(unitcell, path)) == Bandstructure

        # test the calculations (resolution set)
        @test typeof(getBandStructureAlongPath(unitcell, path, resolution=60)) == Bandstructure

        # test the calculations (hermitian)
        @test typeof(getBandStructureAlongPath(unitcell, path, enforce_hermitian=true)) == Bandstructure

    # end the testset here
    end;

# end the testset here
end;






# end the testset here
end;
