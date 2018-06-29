# ALL TESTS THAT CONCERN UNITCELLS

# include
using LatticePhysics
using Base.Test

# begin the Unitcell testblock
unitcell_testset = @testset "Unitcells" begin


################################################################################
#
#   TYPE DEFINITIONS OF OBJECT CLASSES IN JULIA
#
################################################################################

# begin the testset
@testset "General Unitcell functions" begin

    # get a test unitcell
    unitcell = getUnitcellHoneycomb()

    @testset "File access with JLD" begin
        # test saving
        @test typeof(saveUnitcell(unitcell))==String
        # test loading
        @test typeof(loadUnitcell(unitcell))==String
        @test typeof(loadUnitcell(unitcell.filename))==Unitcell
    end;


    @testset "Test and Information functions" begin
        # test printing information
        @test printInfo(unitcell)==nothing
        @test printInfo(unitcell, detailed=true)==nothing
        # test testing method
        @test typeof(testUnitcell(unitcell, 2, 2))==Bool
    end;


    @testset "Connectivity information functions" begin
        # testing connectivity information
        @test typeof(getConnectivityList(unitcell)) == Array{Array,1}
        # testing connection information
        @test typeof(getConnectionList(unitcell)) == Array{Array,1}
        # testing connection strength information
        @test typeof(getConnectionStrengthList(unitcell)) == Array{Any,1}
    end;

# end the testset here
end;








################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT UNITCELLS in 2D and 3D
#
################################################################################

# begin the testset
@testset "Implemented Unitcells in 2D" begin

    #   2D UNITCELLS
    #   - SQUARE / RECTANGLE
    #       - getUnitcellSquare
    @testset "2D Square unitcell" begin
        @test testUnitcell(getUnitcellSquare(), 2, 2)
        @test testUnitcell(getUnitcellSquare(1), 2, 2)
        @test testUnitcell(getUnitcellSquare(2), 2, 2)
        @test testUnitcell(getUnitcellSquare(3), 2, 2)
    end;
    #       - getUnitcellExtendedSquare
    @testset "2D Extended Square unitcell" begin
        @test testUnitcell(getUnitcellExtendedSquare(), 2, 2)
        @test testUnitcell(getUnitcellExtendedSquare(1), 2, 2)
    end;
    #       - getUnitcellCheckerboard
    @testset "2D Checkerboard unitcell" begin
        @test testUnitcell(getUnitcellCheckerboard(), 2, 2)
        @test testUnitcell(getUnitcellCheckerboard(1), 2, 2)
    end;
    #       - getUnitcellShastrySutherland
    @testset "2D Shastry-Sutherland unitcell" begin
        @test testUnitcell(getUnitcellShastrySutherland(), 2, 2)
        @test testUnitcell(getUnitcellShastrySutherland(1), 2, 2)
        @test testUnitcell(getUnitcellShastrySutherland(4), 2, 2)
    end;
    #       - getUnitcellAdvancedSquare
    @testset "2D Advanced Square unitcell" begin
        @test testUnitcell(getUnitcellAdvancedSquare(), 2, 2)
        @test testUnitcell(getUnitcellAdvancedSquare(1), 2, 2)
    end;
    #       - getUnitcellSquareOctagon
    @testset "2D Square-octagon unitcell" begin
        @test testUnitcell(getUnitcellSquareOctagon(), 2, 2)
        @test testUnitcell(getUnitcellSquareOctagon(1), 2, 2)
        @test testUnitcell(getUnitcellSquareOctagon(4), 2, 2)
        @test testUnitcell(getUnitcellSquareOctagon(2), 2, 2)
        @test testUnitcell(getUnitcellSquareOctagon(5), 2, 2)
    end;
    #       - getUnitcellBCC2D
    @testset "2D BCC (square) unitcell" begin
        @test testUnitcell(getUnitcellBCC2D(), 2, 2)
        @test testUnitcell(getUnitcellBCC2D(1), 2, 2)
    end;
    #       - getUnitcellFullyConnectedSquare
    @testset "2D Fully connected square unitcell" begin
        @test testUnitcell(getUnitcellFullyConnectedSquare(), 2, 2)
        @test testUnitcell(getUnitcellFullyConnectedSquare(1), 2, 2)
    end;


    #   - TRIANGULAR
    #       - getUnitcellTriangular
    @testset "2D Triangular unitcell" begin
        @test testUnitcell(getUnitcellTriangular(), 2, 2)
        @test testUnitcell(getUnitcellTriangular(1), 2, 2)
        @test testUnitcell(getUnitcellTriangular(3), 2, 2)
    end;
    #       - getUnitcellHoneycomb
    @testset "2D Honeycomb unitcell" begin
        @test testUnitcell(getUnitcellHoneycomb(), 2, 2)
        @test testUnitcell(getUnitcellHoneycomb(1), 2, 2)
        @test testUnitcell(getUnitcellHoneycomb(4), 2, 2)
        @test testUnitcell(getUnitcellHoneycomb(2), 2, 2)
        @test testUnitcell(getUnitcellHoneycomb(5), 2, 2)
        @test testUnitcell(getUnitcellHoneycomb(3), 2, 2)
    end;
    #       - getUnitcellKagome
    @testset "2D Kagome unitcell" begin
        @test testUnitcell(getUnitcellKagome(), 2, 2)
        @test testUnitcell(getUnitcellKagome(1), 2, 2)
    end;
    #       - getUnitcellKagomeMinus
    @testset "2D Kagome Minus unitcell" begin
        @test testUnitcell(getUnitcellKagomeMinus(), 2, 2)
        @test testUnitcell(getUnitcellKagomeMinus(1), 2, 2)
        @test testUnitcell(getUnitcellKagomeMinus(4), 2, 2)
    end;
    #       - getUnitcellHoneycombXXX
    @testset "2D Honeycomb-XXX unitcell" begin
        @test testUnitcell(getUnitcellHoneycombXXX(), 2, 2)
        @test testUnitcell(getUnitcellHoneycombXXX(1), 2, 2)
        @test testUnitcell(getUnitcellHoneycombXXX(2), 2, 2)
    end;

# end the testset here
end;



# begin the testset
@testset "Implemented Unitcells in 3D" begin

    #   3D UNITCELLS
    #   - CUBIC / FCC
    #       - getUnitcellDiamond
    @testset "3D Diamond unitcell" begin
        @test testUnitcell(getUnitcellDiamond(), 3, 3)
        @test testUnitcell(getUnitcellDiamond(1), 3, 3)
        @test testUnitcell(getUnitcellDiamond(2), 3, 3)
        @test testUnitcell(getUnitcellDiamond(3), 3, 3)
    end;
    #       - getUnitcellBCC
    @testset "3D BCC unitcell" begin
        @test testUnitcell(getUnitcellBCC(), 3, 3)
        @test testUnitcell(getUnitcellBCC(1), 3, 3)
    end;
    #       - getUnitcellPyrochlore
    @testset "3D Pyrochlore unitcell" begin
        @test testUnitcell(getUnitcellPyrochlore(), 3, 3)
        @test testUnitcell(getUnitcellPyrochlore(1), 3, 3)
    end;

    #   - (X,3)y FAMILY
    #       - getUnitcell_8_3_a
    @testset "3D (8,3)a unitcell" begin
        @test testUnitcell(getUnitcell_8_3_a(), 3, 3)
        @test testUnitcell(getUnitcell_8_3_a(1), 3, 3)
        @test testUnitcell(getUnitcell_8_3_a(4), 3, 3)
    end;
    #       - getUnitcell_8_3_b
    @testset "3D (8,3)b unitcell" begin
        @test testUnitcell(getUnitcell_8_3_b(), 3, 3)
        @test testUnitcell(getUnitcell_8_3_b(1), 3, 3)
        @test testUnitcell(getUnitcell_8_3_b(4), 3, 3)
    end;
    #       - getUnitcell_8_3_c
    @testset "3D (8,3)c unitcell" begin
        @test testUnitcell(getUnitcell_8_3_c(), 3, 3)
        @test testUnitcell(getUnitcell_8_3_c(1), 3, 3)
        @test testUnitcell(getUnitcell_8_3_c(4), 3, 3)
    end;
    #       - getUnitcell_8_3_n
    @testset "3D (8,3)n unitcell" begin
        @test testUnitcell(getUnitcell_8_3_n(), 3, 3)
        @test testUnitcell(getUnitcell_8_3_n(1), 3, 3)
        @test testUnitcell(getUnitcell_8_3_n(4), 3, 3)
    end;
    #       - getUnitcell_9_3_a
    @testset "3D (9,3)a unitcell" begin
        @test testUnitcell(getUnitcell_9_3_a(), 3, 3)
        @test testUnitcell(getUnitcell_9_3_a(1), 3, 3)
        @test testUnitcell(getUnitcell_9_3_a(4), 3, 3)
        @test testUnitcell(getUnitcell_9_3_a(2), 3, 3)
        @test testUnitcell(getUnitcell_9_3_a(5), 3, 3)
    end;
    #       - getUnitcell_10_3_a / getUnitcellHyperoctagon
    @testset "3D (10,3)a / hyperoctagon unitcell" begin
        @test testUnitcell(getUnitcell_10_3_a(), 3, 3)
        @test testUnitcell(getUnitcell_10_3_a(1), 3, 3)
        @test testUnitcell(getUnitcell_10_3_a(4), 3, 3)
        @test testUnitcell(getUnitcellHyperoctagon(), 3, 3)
        @test testUnitcell(getUnitcellHyperoctagon(1), 3, 3)
        @test testUnitcell(getUnitcellHyperoctagon(4), 3, 3)
    end;
    #       - getUnitcell_10_3_b / getUnitcellHyperhoneycomb
    @testset "3D (10,3)b / hyperhoneycomb unitcell" begin
        @test testUnitcell(getUnitcell_10_3_b(), 3, 3)
        @test testUnitcell(getUnitcell_10_3_b(1), 3, 3)
        @test testUnitcell(getUnitcell_10_3_b(4), 3, 3)
        @test testUnitcell(getUnitcell_10_3_b(2), 3, 3)
        @test testUnitcell(getUnitcell_10_3_b(5), 3, 3)
        @test testUnitcell(getUnitcellHyperhoneycomb(), 3, 3)
        @test testUnitcell(getUnitcellHyperhoneycomb(1), 3, 3)
        @test testUnitcell(getUnitcellHyperhoneycomb(4), 3, 3)
        @test testUnitcell(getUnitcellHyperhoneycomb(2), 3, 3)
        @test testUnitcell(getUnitcellHyperhoneycomb(5), 3, 3)
    end;
    #       - getUnitcell_10_3_c
    @testset "3D (10,3)c unitcell" begin
        @test testUnitcell(getUnitcell_10_3_c(), 3, 3)
        @test testUnitcell(getUnitcell_10_3_c(1), 3, 3)
        @test testUnitcell(getUnitcell_10_3_c(4), 3, 3)
    end;
    #       - getUnitcell_10_3_d
    @testset "3D (10,3)d unitcell" begin
        @test testUnitcell(getUnitcell_10_3_d(), 3, 3)
        @test testUnitcell(getUnitcell_10_3_d(1), 3, 3)
        @test testUnitcell(getUnitcell_10_3_d(4), 3, 3)
    end;

# end the testset here
end;






################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT UNITCELL CONSTRUCTION RELATED STUFF
#
################################################################################

# begin the testset
@testset "Unitcell construction related" begin

    # get test unitcells
    unitcell_2d = getUnitcellHoneycomb()
    unitcell_3d = getUnitcellDiamond()

    # test getting unitcell from sites
    @testset "Getting unitcells from collection of sites" begin
        @test testUnitcell(getUnitcellFromSites(unitcell_2d.basis, unitcell_2d.lattice_vectors), 2, 2)
        @test testUnitcell(getUnitcellFromSites(unitcell_3d.basis, unitcell_3d.lattice_vectors), 3, 3)
    end;

    # print unitcell code
    @testset "Printing unitcell generating code" begin
        @test printUnitcellGeneratingCode(unitcell_2d)==nothing
        @test printUnitcellGeneratingCode(unitcell_3d)==nothing
    end;

# end the testset here
end;




# end the unitcell test block here
end;
