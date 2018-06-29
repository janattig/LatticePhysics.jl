# ALL TESTS THAT CONCERN MODFICIATION

# include
using LatticePhysics
using Base.Test

# begin the Modification testblock
path_testset = @testset "Path (in momentum space) tests" begin

################################################################################
#
#   METHODS FOR CONSTRUCTION OF PATHS INSIDE THE BZ OF A UNITCELL
#
################################################################################

################################################################################
#
#   1) TYPE PATH
#       - type definition
#       - printInfo function
#       - path String function (e.g. X-G-M-X)
#
################################################################################

# begin the testset
@testset "Type definition and information functions" begin

    # Constructors - how to construct a path
    @testset "Constructors" begin

        # construct some points
        points = Array{Float64,1}[
            [0.0, 0.0],
            [1.0, 1.0],
            [-pi, 0.0],
            [ pi,  pi]
        ]
        # construct some point names
        point_names = String[
            "A",
            "B",
            "A'",
            "something"
        ]
        # construct segment resolution
        segments = Int64[
            100,
            200,
            1200
        ]

        # Test the three possible constructors
        @test typeof(Path(points, point_names, segments)) == Path
        @test typeof(Path(points, point_names)) == Path
        @test typeof(Path()) == Path

    # end the testset here
    end;

    # Information functions printInfo and String
    @testset "Information functions" begin

        # construct some points
        points = Array{Float64,1}[
            [0.0, 0.0],
            [1.0, 1.0],
            [-pi, 0.0],
            [ pi,  pi]
        ]
        # construct some point names
        point_names = String[
            "A",
            "B",
            "A'",
            "something"
        ]
        # construct segment resolution
        segments = Int64[
            100,
            200,
            1200
        ]
        # construct the path
        path = Path(points, point_names, segments)

        # Test the path string
        @test typeof(getPathString(path)) == String
        # Test the printInfo function
        @test printInfo(path) == nothing
        @test printInfo(path, detailed=true) == nothing

    # end the testset here
    end;

# end the testset here
end;



################################################################################
#
#   2) CONSTRUCTION FUNCTIONS PATH
#       - function to add points
#       - TODO function to remove points
#       - function to scale the total resolution
#       - function to set the total resolution
#
################################################################################

# begin the testset
@testset "Construction functions for paths" begin

    # How to add points to a path
    @testset "Adding points" begin

        # construct an empty path
        path = Path()

        # Test three point addings
        @test addPointToPath!(path, [0.0, 0.0], "Gamma") == nothing
        @test addPointToPath!(path, [ pi,  pi], "K") == nothing
        @test addPointToPath!(path, [-pi, 1.0], "A", 200) == nothing

    # end the testset here
    end;

    # Modifying the segment resolution
    @testset "Changing resolution" begin

        # construct some points
        points = Array{Float64,1}[
            [0.0, 0.0],
            [1.0, 1.0],
            [-pi, 0.0],
            [ pi,  pi]
        ]
        # construct some point names
        point_names = String[
            "A",
            "B",
            "A'",
            "something"
        ]
        # construct the path (resolution up to hear: 300)
        path = Path(points, point_names)

        # Test to scale the total resolution by 2
        @test scaleResolution!(path, 2.0) == nothing
        @test sum(path.segment_resolution) == 600

        # Test to scale the total resolution by 0.25
        @test scaleResolution!(path, 0.25) == nothing
        @test sum(path.segment_resolution) == 150


        # construct the path again (resolution now: 300)
        path = Path(points, point_names)

        # Test to set the total resolution to 900
        @test setTotalResolution!(path, 900) == nothing
        @test sum(path.segment_resolution) == 900

        # Test to set the total resolution to 450
        @test setTotalResolution!(path, 450) == nothing
        @test sum(path.segment_resolution) == 450

    # end the testset here
    end;

# end the testset here
end;



################################################################################
#
#   3) DEFAULT PATHS
#       - getDefaultPathTriangular
#       - getDefaultPathSquare
#       - getDefaultPathSquareOctagon
#       - getDefaultPathFCC
#
################################################################################

# begin the testset
@testset "Default Paths" begin

    @testset "Triangular" begin
        # just get the path and test it
        @test typeof(getDefaultPathTriangular()) == Path
        # just get the path with a fixed resolution and test it
        @test typeof(getDefaultPathTriangular(resolution=1000)) == Path
    end;

    @testset "Square" begin
        # just get the path and test it
        @test typeof(getDefaultPathSquare()) == Path
        # get different versions of the path and test them
        @test typeof(getDefaultPathSquare(1)) == Path
        @test typeof(getDefaultPathSquare(2)) == Path
        # just get the path with a fixed resolution and test it
        @test typeof(getDefaultPathSquare(resolution=1000)) == Path
    end;

    @testset "Square Octagon" begin
        # just get the path and test it
        @test typeof(getDefaultPathSquareOctagon()) == Path
        # just get the path with a fixed resolution and test it
        @test typeof(getDefaultPathSquareOctagon(resolution=1000)) == Path
    end;

    @testset "FCC" begin
        # just get the path and test it
        @test typeof(getDefaultPathFCC()) == Path
        # just get the path with a fixed resolution and test it
        @test typeof(getDefaultPathFCC(resolution=1000)) == Path
    end;

# end the testset here
end;





# end the path test block here
end;
