# ALL TESTS THAT CONCERN MODFICIATION

# include
using LatticePhysics
using Base.Test

# begin the Modification testblock
path_testset = @testset "Interaction Matrices" begin

################################################################################
#
#   METHODS FOR CONSTRUCTION INTERACTION MATRICES FOR LATTICES
#
################################################################################

################################################################################
#
#   1) INTERACTION MATRICES IN REAL SPACE
#
################################################################################

# begin the testset
@testset "Interaction matrices in real space" begin

    # begin the testset
    @testset "2D (honeycomb)" begin

        # get a unitcell
        unitcell = getUnitcellHoneycomb()
        # get a lattice
        lattice = getLattice(unitcell, [-4,-2])

        # test the function for unitcell
        @test typeof(getInteractionMatrixRealSpace(unitcell)) == Array{Complex, 2}
        @test ishermitian(getInteractionMatrixRealSpace(unitcell, enforce_hermitian=true))
        @test size(getInteractionMatrixRealSpace(unitcell)) == (2,2)

        # test the function for lattice
        @test typeof(getInteractionMatrixRealSpace(lattice)) == Array{Complex, 2}
        @test ishermitian(getInteractionMatrixRealSpace(lattice, enforce_hermitian=true))
        @test size(getInteractionMatrixRealSpace(lattice)) == (16,16)

    # end the testset here
    end;


    # begin the testset
    @testset "2D (square)" begin

        # get a unitcell
        unitcell = getUnitcellSquare()
        # get a lattice
        lattice = getLattice(unitcell, [-4,-2])

        # test the function for unitcell
        @test typeof(getInteractionMatrixRealSpace(unitcell)) == Array{Complex, 2}
        @test ishermitian(getInteractionMatrixRealSpace(unitcell, enforce_hermitian=true))
        @test size(getInteractionMatrixRealSpace(unitcell)) == (1,1)

        # test the function for lattice
        @test typeof(getInteractionMatrixRealSpace(lattice)) == Array{Complex, 2}
        @test ishermitian(getInteractionMatrixRealSpace(lattice, enforce_hermitian=true))
        @test size(getInteractionMatrixRealSpace(lattice)) == (8,8)

    # end the testset here
    end;


    # begin the testset
    @testset "3D (diamond)" begin

        # get a unitcell
        unitcell = getUnitcellDiamond()
        # get a lattice
        lattice = getLattice(unitcell, [-4,-2,-2])

        # test the function for unitcell
        @test typeof(getInteractionMatrixRealSpace(unitcell)) == Array{Complex, 2}
        @test ishermitian(getInteractionMatrixRealSpace(unitcell, enforce_hermitian=true))
        @test size(getInteractionMatrixRealSpace(unitcell)) == (2,2)

        # test the function for lattice
        @test typeof(getInteractionMatrixRealSpace(lattice)) == Array{Complex, 2}
        @test ishermitian(getInteractionMatrixRealSpace(lattice, enforce_hermitian=true))
        @test size(getInteractionMatrixRealSpace(lattice)) == (32,32)

    # end the testset here
    end;

# end the testset here
end;


################################################################################
#
#   2) INTERACTION MATRICES IN MOMENTUM SPACE
#
################################################################################

# begin the testset
@testset "Interaction matrices in momentum space" begin

    # begin the testset
    @testset "2D (honeycomb)" begin

        # get a unitcell
        unitcell = getUnitcellHoneycomb()
        # get a lattice
        lattice = getLattice(unitcell, [-4,-2])

        # k vector (random)
        k = [0.1, 0.151235]

        # test the function for unitcell
        @test typeof(getInteractionMatrixKSpace(unitcell, k)) == Array{Complex, 2}
        @test ishermitian(getInteractionMatrixKSpace(unitcell, k, enforce_hermitian=true))
        @test size(getInteractionMatrixKSpace(unitcell, k)) == (2,2)

        # test the function for lattice
        @test typeof(getInteractionMatrixKSpace(lattice, k)) == Array{Complex, 2}
        @test ishermitian(getInteractionMatrixKSpace(lattice, k, enforce_hermitian=true))
        @test size(getInteractionMatrixKSpace(lattice, k)) == (16,16)

    # end the testset here
    end;


    # begin the testset
    @testset "2D (square)" begin

        # get a unitcell
        unitcell = getUnitcellSquare()
        # get a lattice
        lattice = getLattice(unitcell, [-4,-2])

        # k vector (random)
        k = [0.1, 0.151235]

        # test the function for unitcell
        @test typeof(getInteractionMatrixKSpace(unitcell, k)) == Array{Complex, 2}
        @test ishermitian(getInteractionMatrixKSpace(unitcell, k, enforce_hermitian=true))
        @test size(getInteractionMatrixKSpace(unitcell, k)) == (1,1)

        # test the function for lattice
        @test typeof(getInteractionMatrixKSpace(lattice, k)) == Array{Complex, 2}
        @test ishermitian(getInteractionMatrixKSpace(lattice, k, enforce_hermitian=true))
        @test size(getInteractionMatrixKSpace(lattice, k)) == (8,8)

    # end the testset here
    end;


    # begin the testset
    @testset "3D (diamond)" begin

        # get a unitcell
        unitcell = getUnitcellDiamond()
        # get a lattice
        lattice = getLattice(unitcell, [-4,-2,-2])

        # k vector (random)
        k = [0.1, 0.151235, -0.21]

        # test the function for unitcell
        @test typeof(getInteractionMatrixKSpace(unitcell, k)) == Array{Complex, 2}
        @test ishermitian(getInteractionMatrixKSpace(unitcell, k, enforce_hermitian=true))
        @test size(getInteractionMatrixKSpace(unitcell, k)) == (2,2)

        # test the function for lattice
        @test typeof(getInteractionMatrixKSpace(lattice, k)) == Array{Complex, 2}
        @test ishermitian(getInteractionMatrixKSpace(lattice, k, enforce_hermitian=true))
        @test size(getInteractionMatrixKSpace(lattice, k)) == (32,32)

    # end the testset here
    end;

# end the testset here
end;





# end the path test block here
end;
