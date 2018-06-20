# ALL TESTS THAT CONCERN UNITCELLS

# include
using LatticePhysics
using Base.Test

# begin the Unitcell testblock
SVG_plotting_testset = @testset "SVG Plotting tests" begin




################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT LATTICE PLOTTING FUNCTIONS FOR
#   PLOTTING TO SVG IMAGES
#
################################################################################


# 2D
@testset "2D Plotting" begin

    # get a test unitcell
    unitcell = getUnitcellHoneycomb(4)
    # get a lattice
    lattice = getLattice(unitcell, [-10, -5])

    # TODO naive method plotLattice

    # Test naive plotting
    @test typeof(plotLattice2D(lattice)) == String
    # Naive plotting and Inkscape pdf export
    @test typeof(plotLattice2D(lattice, inkscape_export_pdf=true)) == String
    # Naive plotting and printing options
    @test typeof(plotLattice2D(lattice, print_used_options=true)) == String

    # custom canvas border size and size
    @test typeof(plotLattice2D(lattice, border_percentage=0.02, size_long_side=2000)) == String

    # periodic connections
    @test typeof(plotLattice2D(lattice, visualize_periodic=true)) == String
    # certain bond thickness
    @test typeof(plotLattice2D(lattice, bond_thickness=-0.2)) == String
    @test typeof(plotLattice2D(lattice, bond_thickness=10.0)) == String
    # certain bond colors
    @test typeof(plotLattice2D(lattice, colorcode_bonds_automation="GREY")) == String
    @test typeof(plotLattice2D(lattice, colorcode_bonds_automation="COLOR")) == String
    @test typeof(plotLattice2D(lattice, colorcode_bonds=Dict("tx"=>[155,0,0], "ty"=>[0,155,0], "tz"=>[0,0,155]))) == String

    # custom site colors
    @test typeof(plotLattice2D(lattice, colorcode_sites=Dict(1=>[155,0,50], 2=>[0,155,50]))) == String
    # custom site labels
    @test typeof(plotLattice2D(lattice, site_labels="POSITION INDEX")) == String
    @test typeof(plotLattice2D(lattice, site_labels="LATTICE INDEX")) == String
    # custom site radius
    @test typeof(plotLattice2D(lattice, site_radius=40.0)) == String
    @test typeof(plotLattice2D(lattice, site_radius=-0.3)) == String
    # cuistom site border
    @test typeof(plotLattice2D(lattice, site_border_width_percentage=0.2)) == String


# end the testset here
end;



# 3D
@testset "3D Plotting" begin

    # get a test unitcell
    unitcell = getUnitcell_10_3_a(4)
    # get a lattice
    lattice = getLattice(unitcell, [-2, -2, -2])


# end the testset here
end;


# 2D Plaquettes
@testset "2D Plaquette Plotting" begin

    # get a test unitcell
    unitcell = getUnitcellHoneycomb(4)
    # get a lattice
    lattice = getLatticeInBox(unitcell, [10.0, 5.0])
    # get the plaquettes
    plaquettes = getPlaquettesOfSite(lattice, 3, 6)
    # get the plaquette values
    plaquette_values = [rand()>0.5 ? 1.0 : -1.0 for p in plaquettes]

    # PLOTTING RELATED (SAME AS 2D PLOTTING)

    # Test naive plotting
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values)) == String
    # Naive plotting and Inkscape pdf export
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, inkscape_export_pdf=true)) == String
    # Naive plotting and printing options
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, print_used_options=true)) == String

    # custom canvas border size and size
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, border_percentage=0.02, size_long_side=2000)) == String

    # periodic connections
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, visualize_periodic=true)) == String
    # certain bond thickness
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, bond_thickness=-0.2)) == String
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, bond_thickness=10.0)) == String
    # certain bond colors
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, colorcode_bonds_automation="GREY")) == String
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, colorcode_bonds_automation="COLOR")) == String
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, colorcode_bonds=Dict("tx"=>[155,0,0], "ty"=>[0,155,0], "tz"=>[0,0,155]))) == String

    # custom site colors
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, colorcode_sites=Dict(1=>[155,0,50], 2=>[0,155,50]))) == String
    # custom site labels
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, site_labels="POSITION INDEX")) == String
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, site_labels="LATTICE INDEX")) == String
    # custom site radius
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, site_radius=40.0)) == String
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, site_radius=-0.3)) == String
    # cuistom site border
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, site_border_width_percentage=0.2)) == String


    # PLAQUETTE RELATED

    # cuistom colorcode
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, colorcode_plaquettes=Dict(-1.0=>[0,125,0], 1.0=>[255,0,125]))) == String
    # cuistom opacity
    @test typeof(plotPlaquettes2D(lattice, plaquettes, plaquette_values, opacity_plaquettes=0.1)) == String



# end the testset here
end;



# end the SVG plotting testset here
end;
