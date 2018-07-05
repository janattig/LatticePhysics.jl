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
#   1) CALCULATION OF FERMI SURFACE
#
#   2) PLOTTING OF FERMI SURFACE
#       - plotting from points
#       - plotting from unitcell
#
################################################################################

# begin the testset
@testset "2D" begin

    # begin the testset only for calculations
    @testset "honeycomb (calculations)" begin

        # get a unitcell
        unitcell = getUnitcellHoneycomb()

        # test the type of the Fermi surface calculations (bare)
        @test typeof(getFermiSurface(unitcell, 20)) == Array{Float64,2}
        @test size(getFermiSurface(unitcell, 20)) == (20,2)

        # test if modifying the fermi energy is possible
        @test typeof(getFermiSurface(unitcell, 20, fermi_energy=1.0)) == Array{Float64,2}

        # test more options
        @test typeof(getFermiSurface(unitcell, 20, enforce_hermitian=true)) == Array{Float64,2}
        @test typeof(getFermiSurface(unitcell, 20, refold_to_first_BZ=false)) == Array{Float64,2}

    # end the testset here
    end;

    # begin the testset for plotting
    @testset "honeycomb (plotting)" begin

        # get a unitcell
        unitcell = getUnitcellHoneycomb()

        # get the fermi surface
        fermi_surface = getFermiSurface(unitcell, 20)

        # test if modifying the fermi energy is possible
        @test_nowarn close(plotFermiSurface(fermi_surface, showPlot=false))
        @test_nowarn close(plotFermiSurface(fermi_surface, plot_title="myplot", showPlot=false))
        @test_nowarn close(plotFermiSurface(fermi_surface, plot_color="r", showPlot=false))
        @test_nowarn close(plotFermiSurface(fermi_surface, figsize=(10,8), showPlot=false))
        @test_nowarn close(plotFermiSurface(fermi_surface, save_filename="$(FOLDER_SPECTRA)testplot.pdf", showPlot=false))

        # test calculation at the same time
        @test_nowarn close(plotFermiSurface(unitcell, 20, showPlot=false))

    # end the testset here
    end;



    # begin the testset
    @testset "square (calculations)" begin

        # get a unitcell
        unitcell = getUnitcellSquare()

        # test the type of the Fermi surface calculations (bare)
        @test typeof(getFermiSurface(unitcell, 20)) == Array{Float64,2}
        @test size(getFermiSurface(unitcell, 20)) == (20,2)

    # end the testset here
    end;

    # begin the testset for plotting
    @testset "square (plotting)" begin

        # get a unitcell
        unitcell = getUnitcellSquare()

        # get the fermi surface
        fermi_surface = getFermiSurface(unitcell, 20)

        # test if modifying the fermi energy is possible
        @test_nowarn close(plotFermiSurface(fermi_surface, showPlot=false))
        @test_nowarn close(plotFermiSurface(fermi_surface, plot_title="myplot", showPlot=false))
        @test_nowarn close(plotFermiSurface(fermi_surface, plot_color="r", showPlot=false))
        @test_nowarn close(plotFermiSurface(fermi_surface, figsize=(10,8), showPlot=false))
        @test_nowarn close(plotFermiSurface(fermi_surface, save_filename="$(FOLDER_SPECTRA)testplot.pdf", showPlot=false))

        # test calculation at the same time
        @test_nowarn close(plotFermiSurface(unitcell, 20, showPlot=false))

    # end the testset here
    end;

# end the testset here
end;













# begin the testset
@testset "3D" begin

    # begin the testset only for calculations
    @testset "diamond (calculations)" begin

        # get a unitcell
        unitcell = getUnitcellDiamond()

        # test the type of the Fermi surface calculations (bare)
        @test typeof(getFermiSurface(unitcell, 20)) == Array{Float64,2}
        @test size(getFermiSurface(unitcell, 20)) == (20,3)

        # test if modifying the fermi energy is possible
        @test typeof(getFermiSurface(unitcell, 20, fermi_energy=1.0)) == Array{Float64,2}

        # test more options
        @test typeof(getFermiSurface(unitcell, 20, enforce_hermitian=true)) == Array{Float64,2}
        @test typeof(getFermiSurface(unitcell, 20, refold_to_first_BZ=false)) == Array{Float64,2}

    # end the testset here
    end;

    # begin the testset for plotting
    @testset "diamond (plotting)" begin

        # get a unitcell
        unitcell = getUnitcellDiamond()

        # get the fermi surface
        fermi_surface = getFermiSurface(unitcell, 20)

        # test if modifying the fermi energy is possible
        @test_nowarn close(plotFermiSurface(fermi_surface, showPlot=false))
        @test_nowarn close(plotFermiSurface(fermi_surface, plot_title="myplot", showPlot=false))
        @test_nowarn close(plotFermiSurface(fermi_surface, plot_color="r", showPlot=false))
        @test_nowarn close(plotFermiSurface(fermi_surface, figsize=(10,8), showPlot=false))
        @test_nowarn close(plotFermiSurface(fermi_surface, save_filename="$(FOLDER_SPECTRA)testplot.pdf", showPlot=false))

        # test calculation at the same time
        @test_nowarn close(plotFermiSurface(unitcell, 20, showPlot=false))

    # end the testset here
    end;


# end the testset here
end;


# end the testset here
end;
