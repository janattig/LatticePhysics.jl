# ALL TESTS THAT CONCERN UNITCELLS

# include
using LatticePhysics
using Base.Test

# begin the Unitcell testblock
plaquette_testset = @testset "Plaquette tests" begin


################################################################################
#
#    FINDING PLAQUETTES IN LATTICE OBJECTS
#
################################################################################

# begin the testset
@testset "Finding Plaquettes" begin

    # First, 2D
    @testset "Plaquettes in 2D (honeycomb)" begin

        # get a test unitcell
        unitcell = getUnitcellHoneycomb()
        # get a test lattice
        lattice = getLatticeByBondDistance(unitcell, 5)


        @testset "Finding plaquettes of one site" for testsite in [1,2,3]
            # does finding give the correct type
            @test typeof(getPlaquettesOfSite(lattice, testsite, 6)) == Array{Array{Int64, 1}, 1}
            # correct number of plaquettes found? 3 hexagons per site in honeycomb but nothing smaller or slightly bigger
            @test length(getPlaquettesOfSite(lattice, testsite, 6)) == 3
            @test sum([length(getPlaquettesOfSite(lattice, testsite, l)) for l in [3,4,5,7,8,9]]) == 0
        end;

        @testset "Finding plaquettes of complete lattice" begin
            # does finding give the correct type
            @test typeof(getPlaquettesOfLattice(lattice, 6)) == Array{Array{Int64, 1}, 1}
            # correct number of plaquettes found? 3 hexagons per site in honeycomb but nothing smaller or slightly bigger
            @test length(getPlaquettesOfLattice(lattice, 6)) > 0
            @test sum([length(getPlaquettesOfLattice(lattice, l)) for l in [3,4,5,7,8,9]]) == 0
        end;

    # end the 2D testset here
    end;



    # Then, 3D
    @testset "Plaquettes in 3D (hyperoctagon)" begin

        # get a test unitcell
        unitcell = getUnitcellHyperoctagon()
        # get a test lattice
        lattice = getLatticeByBondDistance(unitcell, 7)
        # set the testsite
        testsite = 1


        @testset "Finding plaquettes of one site" for testsite in [1,2,3]
            # does finding give the correct type
            @test typeof(getPlaquettesOfSite(lattice, testsite, 10)) == Array{Array{Int64, 1}, 1}
            # correct number of plaquettes found? 15 plaquettes of length 10 per site in hyperoctagon but nothing smaller or slightly bigger
            @test length(getPlaquettesOfSite(lattice, testsite, 10)) == 15
            @test sum([length(getPlaquettesOfSite(lattice, testsite, l)) for l in [3,4,5,6,7,8,9,11,12]]) == 0
        end;

        @testset "Finding plaquettes of complete lattice" begin
            # does finding give the correct type
            @test typeof(getPlaquettesOfLattice(lattice, 10)) == Array{Array{Int64, 1}, 1}
            # correct number of plaquettes found? 3 hexagons per site in honeycomb but nothing smaller or slightly bigger
            @test length(getPlaquettesOfLattice(lattice, 10)) > 0
            @test sum([length(getPlaquettesOfLattice(lattice, l)) for l in [3,4,5,6,7,8,9,11,12]]) == 0
        end;

    # end the 2D testset here
    end;


# end the testset here
end;





# end the plaquette test block here
end;
