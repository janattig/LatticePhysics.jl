################################################################################
#
#   METHODS FOR CONSTRUCTION OF PATHS INSIDE THE BZ OF A UNITCELL
#
#   STRUCTURE OF THE FILE
#
#   1) TYPE PATH
#       - type definition
#       - printInfo function
#
#   2) TODO CONSTRUCTION FUNCTIONS PATH
#       - TODO function to add points
#       - TODO function to remove points
#       - TODO function to set the total resolution
#
#   3) TODO DEFAULT PATHS
#
################################################################################


################################################################################
#
#   The type Path
#
################################################################################
"""
    mutable struct Path

The type that contains information on a path (in momentum space). Fields are

    points             :: Array{Array{Float64, 1}, 1}
    point_names        :: Array{String, 1}
    segment_resolution :: Array{Int64, 1}

Note that the length of the `segment_resolution` array should be smaller than
the length of the `points` array by exactly 1.


New `Path` objects can be created with

    Path(points::Array{Array{Float64,1},1}, point_names::Array{String,1}, segment_resolution::Array{Int64,1})
    Path(points::Array{Array{Float64,1},1}, point_names::Array{String,1})
    Path()

or by one of the several default functions to create a default path.





# Examples

```julia-repl
julia> path = Path()
```
"""
mutable struct Path

    # Array of Point coordinates
    points::Array{Array{Float64,1}, 1}

    # Array of Point names
    point_names::Array{String, 1}

    # Array of resolutions (for later calculations)
    segment_resolution::Array{Int64, 1}





    # The Default constructor
    function Path(points::Array{Array{Float64,1},1}, point_names::Array{String,1}, segment_resolution::Array{Int64,1})
        return new(points, point_names, segment_resolution)
    end

    # The constructor for a new path without segment resolution information
    function Path(points::Array{Array{Float64,1},1}, point_names::Array{String,1})
        return new(points, point_names, ones(Int64, length(point_names)-1).*100)
    end

    # The constructor for a new path without information
    function Path()
        return new(Array{Float64,1}[], String[], Int64[])
    end

end







################################################################################
#
#   Some functions to print information on the path
#
################################################################################

# Function to print some information on a path
function printInfo(path::Path, detailed::Bool=false)
    # distinguish detailed vs. non-detailed
    if detailed
        # print the complete path
        println("Path overview:")
        # alternate between points and segments
        for i in 1:length(path.points)-1
            # print the point
            println("  ($(i)) $(path.point_names[i]) at $(path.points[i])")
            # print the outgoing segment
            println("         |  (resolution: $(path.segment_resolution[i]))")
        end
        # print the last point
        println("  ($(length(points))) $(path.point_names[length(points)]) at $(path.points[length(points)])")
    else
        # not detailed, just give the number of segments and the total resolution
        println("Path contains $(length(path.points)) points ($(length(path.segment_resolution)) segments) with a total resolution of $(sum(path.segment_resolution)).")
    end
end







################################################################################
#
#   CONSTRUCTION OF PATHS
#
################################################################################


# Add a point
function addPointToPath!(path::Path, point::Array{Float64,1}, point_name::String, resolution::Int64=100)
    # push the values into the lists
    push!(path.points, point)
    push!(path.point_names, point_name)
    # maybe push resolution, if there were already some points
    if length(path.points) > 1
        push!(path.segment_resolution, resolution)
    end
    # return nothing
    return nothing
end



################################################################################
#
#   DEFAULT PATHS
#
################################################################################

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
