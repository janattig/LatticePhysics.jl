################################################################################
#
#   METHODS FOR CONSTRUCTION OF PATHS INSIDE THE BZ OF A UNITCELL
#
#   STRUCTURE OF THE FILE
#
#   1) TYPE PATH
#
#   2) CONSTRUCTION FUNCTIONS PATH
#
#   3) DEFAULT PATHS
#
################################################################################



# The type Path
mutable struct Path

    # Array of Point names
    point_names::Array{String, 1}

    # Array of Point coordinates
    points::Array{Array{Float64,1}, 1}

    # Array of resolutions (for later calculations)
    segment_resolution::Array{Int64, 1}





    # The Default constructor
    function Path(point_names::Array{String,1}, points::Array{Array{Float64,1},1}, segment_resolution::Array{Int64,1})
        return new(point_names, points, segment_resolution)
    end

    # The constructor for a new path without segment resolution information
    function Path(point_names::Array{String,1}, points::Array{Array{Float64,1},1})
        return new(point_names, points, ones(Int64, length(point_names)-1).*100)
    end

    # The constructor for a new path without information
    function Path()
        return new(String[], Array{Float64,1}[], Int64[])
    end

end














# SOME DEFAULT PATHS
DEFAULT_PATH_FCC = Array[
    ["gamma"; [0,0,0]],
    ["X"; [2*pi, 0, 0]],
    ["W"; [2*pi, pi, 0]],
    ["L"; [pi, pi, pi]],
    ["gamma"; [0,0,0]],
    ["K"; [3*pi/2, 3*pi/2, 0]],
    ["X"; [2*pi, 0, 0]]
]
export DEFAULT_PATH_FCC

DEFAULT_PATH_TRIANGULAR = Array[
    ["gamma"; [0,0]],
    ["K"; [2*pi/sqrt(3.0), 2*pi/3]],
    ["M"; [2*pi/sqrt(3.0), 0]],
    ["gamma"; [0,0]]
]
export DEFAULT_PATH_TRIANGULAR

DEFAULT_PATH_SQUAREOCTAGON_2 = Array[
    ["gamma"; [0,0]],
    ["K"; [2,  0].*(pi / (1.0 + 1.0/sqrt(2.0)))],
    ["M"; [1, -1].*(pi / (1.0 + 1.0/sqrt(2.0)))],
    ["gamma"; [0,0]]
]
export DEFAULT_PATH_SQUAREOCTAGON_2


DEFAULT_PATH_SQUARE_LONG = Array[
    ["M";     [pi,   0.0]],
    ["Gamma"; [0.0,  0.0]],
    ["K'";    [pi,   -pi]],
    ["M";     [pi,   0.0]],
    ["M'";    [0.0,   pi]],
    ["K";     [pi,    pi]],
    ["Gamma"; [0.0,  0.0]]
]
export DEFAULT_PATH_SQUARE_LONG

DEFAULT_PATH_SQUARE_SHORT = Array[
    ["Gamma"; [0.0,  0.0]],
    ["M";     [pi,   0.0]],
    ["K";     [pi,    pi]],
    ["Gamma"; [0.0,  0.0]]
]
export DEFAULT_PATH_SQUARE_SHORT

DEFAULT_PATH_SQUARE = DEFAULT_PATH_SQUARE_SHORT
export DEFAULT_PATH_SQUARE
