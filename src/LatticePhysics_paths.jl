################################################################################
#
#   METHODS FOR CONSTRUCTION OF PATHS INSIDE THE BZ OF A UNITCELL
#
#   STRUCTURE OF THE FILE
#
#   1) TYPE PATH
#       - type definition
#       - printInfo function
#       - path String function (e.g. X-G-M-X)
#
#   2) CONSTRUCTION FUNCTIONS PATH
#       - function to add points
#       - TODO function to remove points
#       - function to scale the total resolution
#       - function to set the total resolution
#
#   3) DEFAULT PATHS
#       - getDefaultPathTriangular
#       - getDefaultPathSquare
#       - getDefaultPathSquareOctagon
#       - getDefaultPathFCC
#
#   4) TODO FUNCTION TO CREATE A DEFAULT PATH BASED ON A UNITCELL OBJECT
#
#   5) FUNCTION TO PLOT A PATH (in both 2D and 3D)
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


# export the type
export Path






################################################################################
#
#   Some functions to print information on the path
#
################################################################################

# Function to print some information on a path
"""
    printInfo(path::Path [; detailed::Bool])

prints (detailed) information on a `Path` object `path`. If detailed output is desired, the complete path will be printed.


# Examples

```julia-repl
julia> printInfo(path)
...

julia> printInfo(path, detailed=true)
...
```
"""
function printInfo(path::Path; detailed::Bool=false)
    # distinguish detailed vs. non-detailed
    if detailed
        # print the complete path
        println("Path overview:")
        # maybe already abort if no points or only one point in path
        if length(path.points) == 0
            println("   (no points defined)")
            return nothing
        elseif length(path.points) == 1
            println("  ($(1)) \"$(path.point_names[1])\" at $(path.points[1])  (only point in path)")
            return nothing
        end
        # alternate between points and segments
        for i in 1:length(path.points)-1
            # print the point
            println("  ($(i)) \"$(path.point_names[i])\" at $(path.points[i])")
            # print the outgoing segment
            println("         |  (resolution: $(path.segment_resolution[i]))")
        end
        # print the last point
        println("  ($(length(path.points))) \"$(path.point_names[length(path.points)])\" at $(path.points[length(path.points)])")
    else
        # not detailed, just give the number of segments and the total resolution
        println("Path contains $(length(path.points)) points ($(length(path.segment_resolution)) segments) with a total resolution of $(sum(path.segment_resolution)).")
        println("($(getPathString(path)))")
    end
end

# export the function
export printInfo



"""
    getPathString(path::Path)

compiles a string that represents the given path by chaining together all point names.


# Examples

```julia-repl
julia> getPathString(path)
"X--Gamma--K--M--X"

```
"""
function getPathString(path::Path)
    # create a new string
    path_string = ""
    # distinguish by point number
    if length(path.points) == 0
        # no points in path
        path_string = "(no points defined)"
    elseif length(path.points) == 1
        # one point in path
        path_string = "$(path.point_names[1])"
    else
        # alternate between points and segments
        for i in 1:length(path.points)-1
            # add the point and a segment
            path_string = "$(path_string)$(path.point_names[i])--"
        end
        # add the last point
        path_string = "$(path_string)$(path.point_names[end])"
    end
    # return the string
    return path_string
end

# export the function
export getPathString


################################################################################
#
#   CONSTRUCTION OF PATHS
#
################################################################################


# Add a point
"""
    addPointToPath!(path::Path, point::Array{Float64,1}, point_name::String [, resolution::Int64=100])

adds a new point to an existing `Path` object. The new point has to be specified by a name and a location.
Optionally, a resolution of the segment to the preceeding point can be given as well.
Note, that only a resolution will be added if there are already points in the path.
In this function, the path object will be changed and no new object will be created.


# Examples

```julia-repl
julia> addPointToPath!(path, [0.0, 0.0], "Gamma")

julia> addPointToPath!(path, [0.0, pi], "M", 150)
```
"""
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
# Function for when points including pi are added (not of type Float64)
function addPointToPath!(path::Path, point::Array, point_name::String, resolution::Int64=100)
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

# export the function
export addPointToPath!







# scale the resolution by some factor
"""
    scaleResolution!(path::Path, factor::Float64)

scales all segment resolutions of the path by a factor and converts them back to `Int64`.
The path object will be changed and no new object will be created.


# Examples

```julia-repl
julia> scaleResolution!(path, 1.5)
```
"""
function scaleResolution!(path::Path, factor::Float64)
    # multiply all segments
    for s in 1:length(path.segment_resolution)
        path.segment_resolution[s] = round(Int64, path.segment_resolution[s]*factor)
    end
end

# export the function
export scaleResolution!


# set the total resolution
"""
    setTotalResolution!(path::Path, resolution::Int64)

scales all segment resolutions of the path by a factor to match the total resolution `resolution`
and converts them back to `Int64`. The new sum over all segments will give approximately `resolution` (up to float/int conversion).
The path object will be changed and no new object will be created.


# Examples

```julia-repl
julia> setTotalResolution!(path, 1500)
```
"""
function setTotalResolution!(path::Path, resolution::Int64)
    # determine the factor
    factor = resolution / sum(path.segment_resolution)
    # apply the factor
    scaleResolution!(path, factor)
end

# export the function
export setTotalResolution!

















################################################################################
#
#   DEFAULT PATHS
#
################################################################################



"""
    getDefaultPathTriangular( [; resolution::Int64=900])

creates the default path for the triangular lattice. Points in this path are

    [0.0, 0.0]               (Gamma)
    [2*pi/sqrt(3.0), 2*pi/3] (K)
    [2*pi/sqrt(3.0), 0.0]    (M)
    [0.0, 0.0]               (Gamma)

Additionally, a resolution can be set so that the entire path in total has this resolution.

# Examples
```julia-repl
julia> path = getDefaultPathTriangular()
LatticePhysics.Path(...)

julia> path = getDefaultPathTriangular(resolution=1200)
LatticePhysics.Path(...)
```
"""
function getDefaultPathTriangular( ; resolution::Int64=900)
    # create a new path object
    path = Path()
    # insert points
    addPointToPath!(path, [0.0, 0.0], "Gamma")
    addPointToPath!(path, [2*pi/sqrt(3.0), 2*pi/3], "K")
    addPointToPath!(path, [2*pi/sqrt(3.0), 0.0], "M")
    addPointToPath!(path, [0.0, 0.0], "Gamma")
    # set the total resolution
    setTotalResolution!(path, resolution)
    # return the path
    return path
end

# export the function
export getDefaultPathTriangular




"""
    getDefaultPathSquare(version::Int64=1 [; resolution::Int64])

creates the default path for the square lattice. Points in this path are dependent on the version.

Additionally, a resolution can be set so that the entire path in total has this resolution.



# Versions

#### Version 1 (DEFAULT) - short

Points are given by

    [0.0, 0.0]  (Gamma)
    [ pi, 0.0]  (M)
    [ pi,  pi]  (K)
    [0.0, 0.0]  (Gamma)


#### Version 2 - long / extended

Points are given by

    [ pi, 0.0]  (M)
    [0.0, 0.0]  (Gamma)
    [ pi, -pi]  (K')
    [ pi, 0.0]  (M)
    [0.0,  pi]  (M')
    [ pi,  pi]  (K)
    [0.0, 0.0]  (Gamma)


# Examples
```julia-repl
julia> path = getDefaultPathSquare()
LatticePhysics.Path(...)

julia> path = getDefaultPathSquare(2)
LatticePhysics.Path(...)

julia> path = getDefaultPathSquare(resolution=1200)
LatticePhysics.Path(...)
```
"""
function getDefaultPathSquare(version::Int64=1; resolution::Int64=1000)
    # create a new path object
    path = Path()
    # distinguish between version
    if version == 1
        # insert points
        addPointToPath!(path, [0.0, 0.0], "Gamma")
        addPointToPath!(path, [ pi, 0.0], "M")
        addPointToPath!(path, [ pi,  pi], "K")
        addPointToPath!(path, [0.0, 0.0], "Gamma")
    elseif version == 2
        # insert points
        addPointToPath!(path, [ pi, 0.0], "M")
        addPointToPath!(path, [0.0, 0.0], "Gamma")
        addPointToPath!(path, [ pi, -pi], "K'")
        addPointToPath!(path, [ pi, 0.0], "M")
        addPointToPath!(path, [0.0,  pi], "M'")
        addPointToPath!(path, [ pi,  pi], "K")
        addPointToPath!(path, [0.0, 0.0], "Gamma")
    else
        println("version $(version) unknown")
    end
    # set the total resolution
    setTotalResolution!(path, resolution)
    # return the path
    return path
end

# export the function
export getDefaultPathSquare








"""
    getDefaultPathSquareOctagon( [; resolution::Int64=900])

creates the default path for the square octagon lattice (version 2). Points in this path are

    [ 0.0,  0.0]                                 (Gamma)
    [ 2.0,  0.0] * (pi / (1.0 + 1.0/sqrt(2.0)))  (K)
    [-1.0, -1.0] * (pi / (1.0 + 1.0/sqrt(2.0)))  (M)
    [ 0.0,  0.0]                                 (Gamma)

Additionally, a resolution can be set so that the entire path in total has this resolution.

# Examples
```julia-repl
julia> path = getDefaultPathSquareOctagon()
LatticePhysics.Path(...)

julia> path = getDefaultPathSquareOctagon(resolution=1200)
LatticePhysics.Path(...)
```
"""
function getDefaultPathSquareOctagon( ; resolution::Int64=900)
    # create a new path object
    path = Path()
    # insert points
    addPointToPath!(path, [0.0, 0.0], "Gamma")
    addPointToPath!(path, [2,  0].*(pi / (1.0 + 1.0/sqrt(2.0))), "K")
    addPointToPath!(path, [1, -1].*(pi / (1.0 + 1.0/sqrt(2.0))), "M")
    addPointToPath!(path, [0.0, 0.0], "Gamma")
    # set the total resolution
    setTotalResolution!(path, resolution)
    # return the path
    return path
end

# export the function
export getDefaultPathSquareOctagon
















"""
    getDefaultPathFCC( [; resolution::Int64=1200])

creates the default path for the FCC lattice (3D). Points in this path are

    [   0.0,    0.0,  0.0]   (Gamma)
    [  2*pi,    0.0,  0.0]   (X)
    [  2*pi,     pi,  0.0]   (W)
    [    pi,     pi,   pi]   (L)
    [   0.0,    0.0,  0.0]   (Gamma)
    [3*pi/2, 3*pi/2,  0.0]   (K)
    [  2*pi,     pi,  0.0]   (W)

Additionally, a resolution can be set so that the entire path in total has this resolution.

# Examples
```julia-repl
julia> path = getDefaultPathFCC()
LatticePhysics.Path(...)

julia> path = getDefaultPathFCC(resolution=600)
LatticePhysics.Path(...)
```
"""
function getDefaultPathFCC( ; resolution::Int64=1200)
    # create a new path object
    path = Path()
    # insert points
    addPointToPath!(path, [   0.0,    0.0,  0.0], "Gamma")
    addPointToPath!(path, [  2*pi,    0.0,  0.0], "X")
    addPointToPath!(path, [  2*pi,     pi,  0.0], "W")
    addPointToPath!(path, [    pi,     pi,   pi], "L")
    addPointToPath!(path, [   0.0,    0.0,  0.0], "Gamma")
    addPointToPath!(path, [3*pi/2, 3*pi/2,  0.0], "K")
    addPointToPath!(path, [  2*pi,     pi,  0.0], "W")
    # set the total resolution
    setTotalResolution!(path, resolution)
    # return the path
    return path
end

# export the function
export getDefaultPathFCC

















# PLOTTING IN 2D (not exported)
function plotPath2D(
            path::Path;
            plot_color="r",
            new_figure::Bool=true,
            showPlot::Bool=true
        )


    ###########################
    #   INITIAL SETTINGS
    ###########################


    # only if new figure is desired
    if new_figure

        # configure plot environment
        rc("font", family="serif")

        # create a new figure
        fig = figure()

    end




    ###########################
    #   PLOT BRILLOUIN ZONE
    ###########################

    # compile lists of x and y values
    x_values = Float64[p[1] for p in path.points]
    y_values = Float64[p[2] for p in path.points]

    # STEP 1 - scatter the points
    scatter(x_values, y_values, color=plot_color)

    # STEP 2 - draw all lines of edges
    plot(
            x_values, y_values, "-",
            color=plot_color
        )

    # STEP 3 - set all labels
    for i in 1:length(path.points)
        text(
            x_values[i], y_values[i], path.point_names[i],
            ha="left", va="top", size=15, color=plot_color
        )
    end



    ###########################
    #   CONFIGURE AXIS & TITLE
    ###########################

    # equal axis
    gca()[:set_aspect]("equal")



    ###########################
    #   FINISH THE PLOT
    ###########################

    # tighten the layout
    tight_layout()

    # maybe show the plot
    if showPlot
        show()
    end

    # return the figure object
    return gcf()

end

# PLOTTING IN 3D (not exported)
function plotPath3D(
            path::Path;
            plot_color="r",
            new_figure::Bool=true,
            showPlot::Bool=true
        )


    ###########################
    #   INITIAL SETTINGS
    ###########################


    # only if new figure is desired
    if new_figure

        # configure plot environment
        rc("font", family="serif")

        # create a new figure
        fig = figure()

    end




    ###########################
    #   PLOT BRILLOUIN ZONE
    ###########################

    # compile lists of x and y values
    x_values = Float64[p[1] for p in path.points]
    y_values = Float64[p[2] for p in path.points]
    z_values = Float64[p[3] for p in path.points]

    # STEP 1 - scatter the points
    scatter3D(x_values, y_values, z_values, color=plot_color)

    # STEP 2 - draw all lines of edges
    plot3D(
            x_values, y_values, z_values, "-",
            color=plot_color
        )

    # STEP 3 - set all labels
    for i in 1:length(path.points)
        text3D(
            x_values[i], y_values[i], z_values[i], path.point_names[i],
            ha="left", va="top", size=15, color=plot_color
        )
    end




    ###########################
    #   CONFIGURE AXIS & TITLE
    ###########################

    # equal axis
    gca()[:set_aspect]("equal")



    ###########################
    #   FINISH THE PLOT
    ###########################

    # tighten the layout
    tight_layout()

    # maybe show the plot
    if showPlot
        show()
    end

    # return the figure object
    return gcf()

end


# GENERAL PLOTTING (exported)

"""
    plotPath(
                path::Path
             [; plot_color="r",
                new_figure::Bool=true,
                showPlot::Bool=true ]
            )

Plots the Path using `PyPlot`. First, all corner points are scattered using `scatter3d`.
Second, all vertices are plotted using `plot3d`.
Third, all labels are put at the location of the points.

For specifying details while plotting, optional keywords can be passed:
- `plot_color` is the color of the BZ in the plot.
- `new_figure` specifies if a new figure should be created before plotting (or if the current open figure should be used)
- `showPlot` specifies if the plot should be opened after plotting.




# Examples

```julia-repl
julia> plotPath(path)
PyPlot.Figure(...)

julia> plotPath(path, plot_color="m")
PyPlot.Figure(...)

```
"""
function plotPath(
            path::Path;
            plot_color="r",
            new_figure::Bool=true,
            showPlot::Bool=true
        )
    # check if points are 3D or 2D
    if length(path.points[1]) == 2
        return plotPath2D(
            path,
            plot_color=plot_color,
            new_figure=new_figure,
            showPlot=showPlot
        )
    elseif length(path.points[1]) == 3
        return plotPath3D(
            path,
            plot_color=plot_color,
            new_figure=new_figure,
            showPlot=showPlot
        )
    else
        println("Bad dimension of k space: $(length(path.points[1]))")
    end

end
export plotPath
