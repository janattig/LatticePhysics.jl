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
#   2) DEFAULT BZ
#       - TODO getDefaultBZSquare
#       - getDefaultBZFCC
#
#   3) FUNCTION TO CREATE A DEFAULT BZ BASED ON A UNITCELL OBJECT
#       - for 2D unitcells
#       - TODO for 3D unitcells
#
#   4) PLOTTING OF A BZ
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

New `BrillouinZone` objects can be created only with the default constructor of giving all field values explicitly.



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


    # DEFAULT BRILLOUIN ZONE
    function BrillouinZone(points::Array{Array{Float64, 1}, 1}, edges::Array{Array{Int64,1}, 1}, faces::Array{Array{Int64,1}, 1})
        return new(points, edges, faces)
    end

    # EMPTY BRILLOUIN ZONE
    function BrillouinZone()
        return new(Array{Float64, 1}[], Array{Int64, 1}[], Array{Int64, 1}[])
    end

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







################################################################################
#
#   CONSTRUCTION OF BRILLOUIN ZONES
#
################################################################################

# CONSTRUCT 2D (not exported)
# line intersection https://rosettacode.org/wiki/Find_the_intersection_of_two_lines#Julia
function createBrillouinZone2D(unitcell::Unitcell; max_ij::Int64=5)

    ##########
    # STEP 1 - Construct the reciprocal lattice vectors b1 and b2
    ##########

    # get the lattice vectors
    a1 = unitcell.lattice_vectors[1]
    a2 = unitcell.lattice_vectors[2]
    # get the reciprocal lattice vectors
    b1 = [a2[2], -a2[1]]
    b2 = [a1[2], -a1[1]]
    # normalize the vectors
    b1 .*= 2*pi/sum(b1.*a1)
    b2 .*= 2*pi/sum(b2.*a2)


    ##########
    # STEP 2 - Construct the reciprocal lattice points that are relevant
    ##########

    # list of reciprocal points
    k_points = Array{Float64, 1}[]

    # build list of points
    for i in -max_ij:max_ij
    for j in -max_ij:max_ij
        # add to the list
        push!(k_points, b1.*i .+ b2.*j)
    end
    end


    ##########
    # STEP 3 - Construct all mid points as well as directions of normals
    ##########

    # list of mid points
    mid_points = Array{Float64,1}[
        k.*0.5  for k in k_points if sum(k.*k) > 1e-5
    ]

    # list of normals
    normals = Array{Float64,1}[
        [k[2], -k[1]] ./ (k[1]*k[1] + k[2]*k[2]) for k in mid_points
    ]


    ##########
    # STEP 4 - Find all intersection points of lines
    ##########

    # list of intersection points
    intersections = Array{Float64, 1}[]

    # find the intersections
    for i in 1:length(mid_points)
    for j in 1:length(mid_points)
        # continue if same index
        if i==j
            continue
        end
        # get data for point 1
        p1 = mid_points[i]
        d1 = normals[i]
        # get data for point 2
        p2 = mid_points[j]
        d2 = normals[j]
        # check the interesection point
        delta_1 =  d1[2]*p1[1] - d1[1]*p1[2]
        delta_2 =  d2[2]*p2[1] - d2[1]*p2[2]
        delta   = -d1[2]*d2[1] + d2[2]*d1[1]
        # push to the intersections (if not devided by zero)
        if abs(delta) > 1e-8
            push!(intersections,[
                (d1[1]*delta_2 - d2[1]*delta_1) / delta,
                (d1[2]*delta_2 - d2[2]*delta_1) / delta
            ])
        end
    end
    end

    ##########
    # STEP 5 - Find closest intersection points --> points of BZ
    ##########

    # list of distances
    intersection_distances = Float64[
        sum(isec.*isec) for isec in intersections
    ]

    # find the minimal intersection distance
    min_distance = minimum(intersection_distances)

    # list of points which are only the minimal intersection distance away
    points = Array{Float64, 1}[]

    # check all points
    for i in 1:length(intersections)
        # check the relative deviation to the calculated minimal distance
        if intersection_distances[i] < min_distance * (1.0 + 1e-8)
            # check if not inserted already
            inserted_already = false
            for p in points
                if sum((p.-intersections[i]).*(p.-intersections[i])) < 1e-8
                    inserted_already = true
                end
            end
            # found a new point
            if !inserted_already
                push!(points, intersections[i])
            end
        end
    end



    ##########
    # STEP 6 - Connect the points of the BZ --> edges of BZ
    ##########

    # list of points if they are contained in the loop
    in_loop = Bool[false for p in points]

    # the list of the loop (to bring the points in correct order)
    loop = Int64[]

    # start with the first point
    push!(loop, 1)
    in_loop[1] = true

    # look for next point (which is closest to the given point in all remaing ones)
    while length(loop) < length(points)
        # check for the closest
        closest_distance = sum(points[1].*points[1]) * 1000
        closest_point = -1
        # check all other points
        for i in 1:length(points)
            # if already in loop, continue
            if in_loop[i]
                continue
            end
            # compare to the nearest one
            if sum((points[i].-points[loop[end]]).*(points[i].-points[loop[end]])) < closest_distance
                closest_point    = i
                closest_distance = sum((points[i].-points[loop[end]]).*(points[i].-points[loop[end]]))
            end
        end
        # insert the closest point as the next in line
        push!(loop, closest_point)
        in_loop[closest_point] = true
    end

    # push the starting point into the loop to wrap up and make a closed line
    push!(loop, loop[1])

    # list of edges (only contains this one loop)
    edges = Array{Int64, 1}[loop]



    ##########
    # STEP 7 - Finish up
    ##########

    # return the finished BZ
    return BrillouinZone(
        points,
        edges,
        Array{Int64, 1}[]
    )
end

# TODO CONSTRUCT 3D (not exported)
# plane intersection: http://www.ambrsoft.com/TrigoCalc/Plan3D/3PlanesIntersection_.htm
function createBrillouinZone3D(unitcell::Unitcell; max_ij::Int64=5)

    ##########
    # STEP 1 - Construct the reciprocal lattice vectors b1, b2 and b3
    ##########

    # get the lattice vectors
    a1 = unitcell.lattice_vectors[1]
    a2 = unitcell.lattice_vectors[2]
    a3 = unitcell.lattice_vectors[3]
    # get the reciprocal lattice vectors
    b1 = cross(a2, a3)
    b2 = cross(a3, a1)
    b3 = cross(a1, a2)
    # normalize the vectors
    b1 .*= 2*pi/sum(b1.*a1)
    b2 .*= 2*pi/sum(b2.*a2)
    b3 .*= 2*pi/sum(b3.*a3)


    ##########
    # STEP 2 - Construct the reciprocal lattice points that are relevant
    ##########

    # list of reciprocal points
    k_points = Array{Float64, 1}[]

    # build list of points
    for i in -max_ij:max_ij
    for j in -max_ij:max_ij
    for l in -max_ij:max_ij
        # add to the list
        push!(k_points, b1.*i .+ b2.*j .+ b3.*l)
    end
    end
    end


    ##########
    # STEP 3 - Construct all mid points as well as directions of normals
    ##########

    # list of mid points
    mid_points = Array{Float64,1}[
        k.*0.5  for k in k_points if sum(k.*k) > 1e-5
    ]

    # list of normals
    normals = deepcopy(mid_points)


    ##########
    # STEP 4 - Find all intersection points of lines
    ##########

    # list of intersection points
    intersections = Array{Float64, 1}[]

    # find the intersections
    for i in 1:length(mid_points)
    for j in 1:length(mid_points)
        # continue if same index
        if i==j
            continue
        end
        # get data for point 1
        p1 = mid_points[i]
        d1 = normals[i]
        # get data for point 2
        p2 = mid_points[j]
        d2 = normals[j]
        # check the interesection point
        delta_1 =  d1[2]*p1[1] - d1[1]*p1[2]
        delta_2 =  d2[2]*p2[1] - d2[1]*p2[2]
        delta   = -d1[2]*d2[1] + d2[2]*d1[1]
        # push to the intersections (if not devided by zero)
        if abs(delta) > 1e-8
            push!(intersections,[
                (d1[1]*delta_2 - d2[1]*delta_1) / delta,
                (d1[2]*delta_2 - d2[2]*delta_1) / delta
            ])
        end
    end
    end

    ##########
    # STEP 5 - Find closest intersection points --> points of BZ
    ##########

    # list of distances
    intersection_distances = Float64[
        sum(isec.*isec) for isec in intersections
    ]

    # find the minimal intersection distance
    min_distance = minimum(intersection_distances)

    # list of points which are only the minimal intersection distance away
    points = Array{Float64, 1}[]

    # check all points
    for i in 1:length(intersections)
        # check the relative deviation to the calculated minimal distance
        if intersection_distances[i] < min_distance * (1.0 + 1e-8)
            # check if not inserted already
            inserted_already = false
            for p in points
                if sum((p.-intersections[i]).*(p.-intersections[i])) < 1e-8
                    inserted_already = true
                end
            end
            # found a new point
            if !inserted_already
                push!(points, intersections[i])
            end
        end
    end



    ##########
    # STEP 6 - Connect the points of the BZ --> edges of BZ
    ##########

    # list of points if they are contained in the loop
    in_loop = Bool[false for p in points]

    # the list of the loop (to bring the points in correct order)
    loop = Int64[]

    # start with the first point
    push!(loop, 1)
    in_loop[1] = true

    # look for next point (which is closest to the given point in all remaing ones)
    while length(loop) < length(points)
        # check for the closest
        closest_distance = sum(points[1].*points[1]) * 1000
        closest_point = -1
        # check all other points
        for i in 1:length(points)
            # if already in loop, continue
            if in_loop[i]
                continue
            end
            # compare to the nearest one
            if sum((points[i].-points[loop[end]]).*(points[i].-points[loop[end]])) < closest_distance
                closest_point    = i
                closest_distance = sum((points[i].-points[loop[end]]).*(points[i].-points[loop[end]]))
            end
        end
        # insert the closest point as the next in line
        push!(loop, closest_point)
        in_loop[closest_point] = true
    end

    # push the starting point into the loop to wrap up and make a closed line
    push!(loop, loop[1])

    # list of edges (only contains this one loop)
    edges = Array{Int64, 1}[loop]



    ##########
    # STEP 7 - Finish up
    ##########

    # return the finished BZ
    return BrillouinZone(
        points,
        edges,
        Array{Int64, 1}[]
    )
    # TODO so far only an empty BZ is returned
    return BrillouinZone()
end




# CONSTRUCT IN ANY DIMENSION
function createBrillouinZone(unitcell::Unitcell; max_ij::Int64=5)
    # distinguish the number of dimensions
    if length(unitcell.lattice_vectors) == 2 && length(unitcell.basis[1]) == 2
        # return the 2D case
        return createBrillouinZone2D(unitcell, max_ij=max_ij)
    elseif length(unitcell.lattice_vectors) == 3 && length(unitcell.basis[1]) == 3
        # return the 3D case
        return createBrillouinZone3D(unitcell, max_ij=max_ij)
    else
        # print an error and return an empty BZ
        println("dimensions not fitting for BZ calculation")
        return BrillouinZone(
            Array{Float64, 1}[],
            Array{Int64, 1}[],
            Array{Int64, 1}[]
        )
    end
end
export createBrillouinZone





################################################################################
#
#   PLOTTING OF BRILLOUIN ZONES
#
################################################################################



# PLOTTING IN 2D (not exported)
function plotBrillouinZone2D(
            brillouin_zone::BrillouinZone;
            filter_points::Bool=true,
            plot_color="k",
            new_figure::Bool=true,
            zoom_to_BZ::Bool=true,
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

    # STEP 1 - scatter the points (maybe filter)
    if filter_points && (length(brillouin_zone.edges)>0 && length(brillouin_zone.faces)>0)
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
        # zoom to the used points as well
        max_dim = 0
        max_dim = max(max_dim, maximum(x_values[used]))
        max_dim = max(max_dim, maximum(y_values[used]))
        max_dim = max(max_dim, maximum(x_values[used].*-1))
        max_dim = max(max_dim, maximum(y_values[used].*-1))
    else
        # no filter, scatter everything
        scatter(x_values, y_values, color=plot_color)
        # zoom to all points
        max_dim = 0
        max_dim = max(max_dim, maximum(x_values))
        max_dim = max(max_dim, maximum(y_values))
        max_dim = max(max_dim, maximum(x_values.*-1))
        max_dim = max(max_dim, maximum(y_values.*-1))
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

    # STEP 3 - plot all surfaces (does not apply for 2D)




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

    # maybe zoom
    if zoom_to_BZ
        xlim(-max_dim*1.05, max_dim*1.05)
        ylim(-max_dim*1.05, max_dim*1.05)
    end

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
            new_figure::Bool=true,
            zoom_to_BZ::Bool=true,
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
    if filter_points && (length(brillouin_zone.edges)>0 && length(brillouin_zone.faces)>0)
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
        # find out maximum dimension
        max_dim = 0
        max_dim = max(max_dim, maximum(x_values[used]))
        max_dim = max(max_dim, maximum(y_values[used]))
        max_dim = max(max_dim, maximum(z_values[used]))
        max_dim = max(max_dim, maximum(x_values[used].*-1))
        max_dim = max(max_dim, maximum(y_values[used].*-1))
        max_dim = max(max_dim, maximum(z_values[used].*-1))
    else
        # no filter, scatter everything
        scatter3D(x_values, y_values, z_values, color=plot_color)
        # find out maximum dimension
        max_dim = 0
        max_dim = max(max_dim, maximum(x_values))
        max_dim = max(max_dim, maximum(y_values))
        max_dim = max(max_dim, maximum(z_values))
        max_dim = max(max_dim, maximum(x_values.*-1))
        max_dim = max(max_dim, maximum(y_values.*-1))
        max_dim = max(max_dim, maximum(z_values.*-1))
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

    # maybe zoom
    if zoom_to_BZ
        xlim(-max_dim*1.05, max_dim*1.05)
        ylim(-max_dim*1.05, max_dim*1.05)
        zlim(-max_dim*1.05, max_dim*1.05)
    end

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
                new_figure::Bool=true,
                zoom_to_BZ::Bool=true,
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
            new_figure::Bool=true,
            zoom_to_BZ::Bool=true,
            showPlot::Bool=true
        )
    # check if it contains points at all
    if length(brillouin_zone.points) == 0
        println("BZ does not contain points")
        return
    end
    # check if points are 3D or 2D
    if length(brillouin_zone.points[1]) == 2
        return plotBrillouinZone2D(
            brillouin_zone,
            filter_points=filter_points,
            plot_color=plot_color,
            new_figure=new_figure,
            zoom_to_BZ=zoom_to_BZ,
            showPlot=showPlot
        )
    elseif length(brillouin_zone.points[1]) == 3
        return plotBrillouinZone3D(
            brillouin_zone,
            filter_points=filter_points,
            plot_color=plot_color,
            new_figure=new_figure,
            zoom_to_BZ=zoom_to_BZ,
            showPlot=showPlot
        )
    else
        println("Bad dimension of k space: $(length(brillouin_zone.points[1]))")
    end

end
export plotBrillouinZone
