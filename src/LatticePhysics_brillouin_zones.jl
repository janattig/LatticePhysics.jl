################################################################################
#
#   METHODS FOR CONSTRUCTION AND PLOTTING OF BRILLOUIN ZONES (BZ) OF A UNITCELL
#
#   STRUCTURE OF THE FILE
#
#   1) TYPE BRILLOUINZONE
#       - type definition
#       - TODO printInfo function
#
#   2) TODO CONSTRUCTION FUNCTIONS BZ
#
#   3) DEFAULT BZ
#       - TODO getDefaultBZSquare
#       - getDefaultBZFCC
#
#   4) TODO FUNCTION TO CREATE A DEFAULT BZ BASED ON A UNITCELL OBJECT
#
#   5) PLOTTING OF A BZ
#
################################################################################



################################################################################
#
#   The type BrillouinZone
#
################################################################################
"""
    struct BrillouinZone

The type that contains information on a Brillouin zone (of a `Unitcell` object). Fields are

    points :: Array{Array{Float64, 1}, 1}
    edges  :: Array{Array{Int64, 1}, 1}
    faces  :: Array{Array{Int64, 1}, 1}

whereas `edges` and `faces` contain lists of indices of points given in the `points` field.

New `BrillouinZone` object
s can be created only with the default constructor of giving all field values explicitly.



# Examples

```julia-repl
julia> bz = BrillouinZone(points, edges, faces)
LatticePhysics.BrillouinZone(...)
```
"""
struct BrillouinZone

    # FIELDS
    # Note: 2D has no faces, ONLY edges

    # the contained points (corners of the BZ)
    points::Array{Array{Float64, 1}, 1}

    # edges of the BZ (lists with combination of point indices)
    edges::Array{Array{Int64,1}, 1}

    # faces of the BZ (lists with combination of point indices)
    faces::Array{Array{Int64,1}, 1}

end

# export the type
export BrillouinZone









################################################################################
#
#   DEFAULT BRILLOUIN ZONES
#
################################################################################



"""
    getDefaultBZFCC()

creates the default Brillouin zone for the FCC lattice (3D) as a `BrillouinZone` object.

# Examples
```julia-repl
julia> bz = getDefaultBZFCC()
LatticePhysics.BrillouinZone(...)
```
"""
function getDefaultBZFCC()
    # create a list of all points
    points = Array{Float64,1}[
        [-2,  0, -1], #0
        [-2,  1,  0], #1
        [-2, -1,  0],
        [-2,  0,  1],
        [-4,  2,  1],
        [-3,  2,  0],
        [-4,  1,  2], #6
        [-3,  0,  2],
        [-2,  4, -1],
        [-2,  3,  0], #9
        [-1,  2,  0],
        [ 0,  2,  1],
        [ 1,  2,  0],
        [ 0,  2, -1], #13
        [-1,  4, -2],
        [ 0,  3, -2],
        [-1,  0, -2],
        [ 0,  1, -2], #17
        [ 0, -1, -2],
        [ 1,  0, -2],
        [ 2,  1, -4], #20
        [ 2,  0, -3],
        [ 1,  2, -4],
        [ 0,  2, -3],
        [-2, -1,  4],
        [-2,  0,  3],
        [-1,  0,  2], #26
        [ 0,  1,  2],
        [ 1,  0,  2],
        [ 0, -1,  2],
        [-1, -2,  4],
        [ 0, -2,  3],
        [-1, -2,  0], #32
        [ 0, -2, -1],
        [ 0, -2,  1],
        [ 1, -2,  0],
        [ 2, -4,  1],
        [ 2, -3,  0],
        [ 1, -4,  2], #38
        [ 0, -3,  2],
        [ 2,  1,  0],
        [ 2,  0, -1],
        [ 2, -1,  0],
        [ 2,  0,  1],
        [ 4, -1, -2],
        [ 3,  0, -2],
        [ 4, -2, -1],
        [ 3, -2,  0]  #47
    ] .* pi
    # create a list of all edges
    edges = Array{Int64,1}[
        [19, 34, 36, 43, 42, 20, 19],
        [1, 17, 18, 14, 11, 2, 1],
        [41, 44, 43, 42, 41],
        [11, 12, 13, 14, 11],
        [3, 33, 35, 30, 27, 4, 3],
        [29, 44, 43, 36, 35, 30, 29],
        [12, 28, 29, 44, 41, 13, 12],
        [1, 3, 4, 2, 1],
        [2, 11, 12, 28, 27, 4, 2],
        [27, 28, 29, 30, 27],
        [13, 41, 42, 20, 18, 14, 13],
        [33, 34, 36, 35, 33],
        [17, 19, 20, 18, 17],
        [1, 17, 19, 34, 33, 3, 1]
    ]
    # create a list of all faces
    faces = deepcopy(edges)
    # create the object
    bz = BrillouinZone(points, edges, faces)
    # return the BZ object
    return bz
end

# export the function
export getDefaultBZFCC















# PLOTTING IN 2D (not exported)
function plotBrillouinZone2D(
            brillouin_zone::BrillouinZone;
            filter_points::Bool=true,
            plot_color="k",
            new_figure=true,
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
        fig = figure(figsize=figsize)

    end




    ###########################
    #   PLOT BRILLOUIN ZONE
    ###########################

    # compile lists of x and y values
    x_values = Float64[p[1] for p in brillouin_zone.points]
    y_values = Float64[p[2] for p in brillouin_zone.points]

    # STEP 1 - scatter the points (maybe filter)
    if filter_points
        # check all indices for usage
        used = Bool[false for x in x_values]
        # iterate over all faces and edges
        for e in brillouin_zone.edges
        for i in e
            used[i] = true
        end
        end
        for f in brillouin_zone.faces
        for i in f
            used[i] = true
        end
        end
        # scatter only filtered points
        scatter(x_values[used], y_values[used], color=plot_color)
    else
        # no filter, scatter everything
        scatter(x_values, y_values, color=plot_color)
    end

    # STEP 2 - draw all lines of edges
    for l in brillouin_zone.edges
        # plot the edge
        plot(
                [x_values[i] for i in l],
                [y_values[i] for i in l],
                color=plot_color
            )
    end

    # STEP 3 - plot all surfaces (does not apply for 3D)




    ###########################
    #   CONFIGURE AXIS & TITLE
    ###########################

    # equal axis and no labels
    gca()[:set_aspect]("equal")
    gca()[:axis]("off")
    axis("off")



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
function plotBrillouinZone3D(
            brillouin_zone::BrillouinZone;
            filter_points::Bool=true,
            plot_color="k",
            new_figure=true,
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
    x_values = Float64[p[1] for p in brillouin_zone.points]
    y_values = Float64[p[2] for p in brillouin_zone.points]
    z_values = Float64[p[3] for p in brillouin_zone.points]

    # STEP 1 - scatter the points (maybe filter)
    if filter_points
        # check all indices for usage
        used = Bool[false for x in x_values]
        # iterate over all faces and edges
        for e in brillouin_zone.edges
        for i in e
            used[i] = true
        end
        end
        for f in brillouin_zone.faces
        for i in f
            used[i] = true
        end
        end
        # scatter only filtered points
        scatter3D(x_values[used], y_values[used], z_values[used], color=plot_color)
    else
        # no filter, scatter everything
        scatter3D(x_values, y_values, z_values, color=plot_color)
    end

    # STEP 2 - draw all lines of edges
    for l in brillouin_zone.edges
        # plot the edge
        plot3D(
                [x_values[i] for i in l],
                [y_values[i] for i in l],
                [z_values[i] for i in l],
                color=plot_color
            )
    end

    # TODO STEP 3 - plot all surfaces




    ###########################
    #   CONFIGURE AXIS & TITLE
    ###########################

    # equal axis and no labels
    gca()[:set_aspect]("equal")
    gca()[:axis]("off")
    axis("off")



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
    plotBrillouinZone(
                brillouin_zone::BrillouinZone
             [; filter_points::Bool=true,
                plot_color="k",
                new_figure=true,
                showPlot::Bool=true ]
            )

Plots the Brillouin zone given by the `BrillouinZone` object `brillouin_zone` using `PyPlot` in three steps.
First, all corner points are scattered using `scatter3d`.
Second, all edges are plotted using `plot3d`.
TODO: Third, the faces of the BZ are plotted.

For specifying details while plotting, optional keywords can be passed:
- `filter_points` determines if only points should be scattered which are used in edges or faces.
- `plot_color` is the color of the BZ in the plot.
- `new_figure` specifies if a new figure should be created before plotting (or if the current open figure should be used)
- `showPlot` specifies if the plot should be opened after plotting.




# Examples

```julia-repl
julia> plotBrillouinZone(bz)
PyPlot.Figure(...)

julia> plotBrillouinZone(bz, plot_color="b", filter_points=false)
PyPlot.Figure(...)

```
"""
function plotBrillouinZone(
            brillouin_zone::BrillouinZone;
            filter_points::Bool=true,
            plot_color="k",
            new_figure=true,
            showPlot::Bool=true
        )
    # check if points are 3D or 2D
    if length(brillouin_zone.points[1]) == 2
        return plotBrillouinZone2D(
            brillouin_zone,
            filter_points=filter_points,
            plot_color=plot_color,
            new_figure=new_figure,
            showPlot=showPlot
        )
    elseif length(brillouin_zone.points[1]) == 3
        return plotBrillouinZone3D(
            brillouin_zone,
            filter_points=filter_points,
            plot_color=plot_color,
            new_figure=new_figure,
            showPlot=showPlot
        )
    else
        println("Bad dimension of k space: $(length(brillouin_zone.points[1]))")
    end

end
export plotBrillouinZone
