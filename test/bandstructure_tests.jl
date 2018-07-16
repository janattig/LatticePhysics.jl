# ALL TESTS THAT CONCERN MODFICIATION

# include
using LatticePhysics
using PyPlot
using Base.Test

# begin the Modification testblock
bandstructure_testset = @testset "Bandstructures" begin

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
@testset "Calculation of bandstructures" begin

    # begin the testset
    @testset "2D (honeycomb)" begin

        # get a unitcell
        unitcell = getUnitcellHoneycomb()
        # get a path
        path = getDefaultPathTriangular(resolution=30)

        # test the calculations (bare)
        @test typeof(getBandStructure(unitcell, path)) == Bandstructure

        # test the calculations (resolution set)
        @test typeof(getBandStructure(unitcell, path, resolution=60)) == Bandstructure

        # test the calculations (hermitian)
        @test typeof(getBandStructure(unitcell, path, enforce_hermitian=true)) == Bandstructure

    # end the testset here
    end;


    # begin the testset
    @testset "2D (square)" begin

        # get a unitcell
        unitcell = getUnitcellSquare()
        # get a path
        path = getDefaultPathSquare(resolution=30)

        # test the calculations (bare)
        @test typeof(getBandStructure(unitcell, path)) == Bandstructure

        # test the calculations (resolution set)
        @test typeof(getBandStructure(unitcell, path, resolution=60)) == Bandstructure

        # test the calculations (hermitian)
        @test typeof(getBandStructure(unitcell, path, enforce_hermitian=true)) == Bandstructure

    # end the testset here
    end;

    # begin the testset
    @testset "3D (diamond)" begin

        # get a unitcell
        unitcell = getUnitcellDiamond()
        # get a path
        path = getDefaultPathFCC(resolution=60)

        # test the calculations (bare)
        @test typeof(getBandStructure(unitcell, path)) == Bandstructure

        # test the calculations (resolution set)
        @test typeof(getBandStructure(unitcell, path, resolution=120)) == Bandstructure

        # test the calculations (hermitian)
        @test typeof(getBandStructure(unitcell, path, enforce_hermitian=true)) == Bandstructure

    # end the testset here
    end;

# end the testset here
end;




################################################################################
#
#   3) PLOTTING OF BAND STRUCTURES
#       - plotting of Bandstructure objects
#       - plotting of bandstructures of unitcells along paths
#
################################################################################


# begin the testset
@testset "Plotting of bandstructures" begin

    # begin the testset
    @testset "2D (honeycomb)" begin

        # get a unitcell
        unitcell = getUnitcellHoneycomb()
        # get a path
        path = getDefaultPathTriangular(resolution=30)

        # get the bandstructure (bare)
        bandstructure = getBandStructure(unitcell, path)

        # test the plotting
        @test_nowarn close(plotBandstructure(bandstructure, showPlot=false))
        @test_nowarn close(plotBandstructure(bandstructure, limits_energy=[-1, 2.0], showPlot=false))
        @test_nowarn close(plotBandstructure(bandstructure, plot_title="myplot", showPlot=false))
        @test_nowarn close(plotBandstructure(bandstructure, plot_color="r", showPlot=false))
        @test_nowarn close(plotBandstructure(bandstructure, figsize=(10,8), showPlot=false))
        @test_nowarn close(plotBandstructure(bandstructure, save_filename="$(FOLDER_SPECTRA)testplot.pdf", showPlot=false))

        # test calculation and plotting
        @test_nowarn close(plotBandstructure(unitcell, path, showPlot=false))
        @test_nowarn close(plotBandstructure(unitcell, path, limits_energy=[-1, 2.0], showPlot=false))
        @test_nowarn close(plotBandstructure(unitcell, path, plot_title="myplot", showPlot=false))
        @test_nowarn close(plotBandstructure(unitcell, path, plot_color="r", showPlot=false))
        @test_nowarn close(plotBandstructure(unitcell, path, figsize=(10,8), showPlot=false))
        @test_nowarn close(plotBandstructure(unitcell, path, save_filename="$(FOLDER_SPECTRA)testplot.pdf", showPlot=false))

    # end the testset here
    end;


    # begin the testset
    @testset "2D (square)" begin

        # get a unitcell
        unitcell = getUnitcellSquare()
        # get a path
        path = getDefaultPathSquare(resolution=30)

        # get the bandstructure (bare)
        bandstructure = getBandStructure(unitcell, path)

        # test the plotting
        @test_nowarn close(plotBandstructure(bandstructure, showPlot=false))

    # end the testset here
    end;

    # begin the testset
    @testset "3D (diamond)" begin

        # get a unitcell
        unitcell = getUnitcellDiamond()
        # get a path
        path = getDefaultPathFCC(resolution=60)

        # get the bandstructure (bare)
        bandstructure = getBandStructure(unitcell, path)

        # test the plotting
        @test_nowarn close(plotBandstructure(bandstructure, showPlot=false))

    # end the testset here
    end;

# end the testset here
end;



# end the testset here
end;
