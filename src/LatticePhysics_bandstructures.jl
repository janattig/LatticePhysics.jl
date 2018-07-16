################################################################################
#
#   METHODS FOR CONSTRUCTION OF BANDSTRUCTURES (ALONG PATHS)
#
#   STRUCTURE OF THE FILE
#
#   1) TYPE BANDSTRUCTURE
#       - type definition
#       - TODO printInfo function
#
#   2) CALCULATION OF BAND STRUCTURES OF UNTICELL OBJECTS
#
#   3) PLOTTING OF BAND STRUCTURES
#       - plotting of Bandstructure objects
#       - plotting of bandstructures of unitcells along paths
#
################################################################################




################################################################################
#
#   TYPE BANDSTRUCTURE
#       - type definition
#       - printInfo function
#
################################################################################
struct Bandstructure

    # the path along which the band structure is calcualted
    path::Path

    # bands for each segment
    # bands[i] gives all bands of segment i
    # bands[i][j] gives all energy values for band j of segment i
    # bands[i][j][k] gives the energy value at kpoint index k of band j in segment i
    bands::Array{Array{Array{Float64, 1}, 1}, 1}

    # ONLY DEFAULT CONSTRUCTOR

end

# export the type
export Bandstructure









################################################################################
#
#   BAND STRUCTURE CALCULATION
#
################################################################################
"""
    getBandStructure(
                unitcell::Unitcell,
                path::Path
             [; resolution::Int64=-1,
                enforce_hermitian::Bool=false ]
            )

    getBandStructure(
                matrixFunction::Function,
                path::Path
             [; resolution::Int64=-1 ]
            )

Calculates the band struture of a `Unitcell` object (or a matrix given by a function `matrixFunction`)
along some path given by a `Path` object and returns a `Bandstructure` object.


# Examples

```julia-repl
julia> bandstructure = getBandStructure(unitcell, path)
LatticePhysics.Bandstructure(...)

julia> bandstructure = getBandStructure(unitcell, path, resolution=1000)
LatticePhysics.Bandstructure(...)
```
"""
function getBandStructure(
                unitcell::Unitcell,
                path::Path;
                resolution::Int64=-1,
                enforce_hermitian::Bool=false
            )

    # maybe modify the path resolution
    if resolution > 0
        setTotalResolution!(path, resolution)
    end

    # build up the list of all bands of all segments
    segments_total = Array{Array{Float64,1},1}[]
    # iterate over all path segments and push empty lists into the segments list
    for i in 1:length(path.segment_resolution)
        # build an empty band structure for this segment
        segment = Array{Float64, 1}[]
        for b in 1:length(unitcell.basis)
            push!(segment, zeros(Float64, path.segment_resolution[i]))
        end
        # push the segment band structure into the complete segment list
        push!(segments_total, segment)
    end


    # iterate over all segments
    for s in 1:length(path.segment_resolution)
        # get the grid in between two points
        segment_resolution = path.segment_resolution[s]
        # get all multipliers of k vectors (i.e. all alpha in (1-alpha)*k_1 + alpha*k_2)
        multipliers = linspace(0, 1, segment_resolution)
        # get the local start and end point of the segment
        k1 = path.points[s]
        k2 = path.points[s+1]
        # calculate all energies
        for i in 1:segment_resolution
            # get the current k
            k = (k2 .* multipliers[i]) .+ (k1 .* (1-multipliers[i]))
            # get the interaction matrix for this k
            matrix = getInteractionMatrixKSpace(unitcell, k, enforce_hermitian=enforce_hermitian)
            # diagonalize the matrix
            eigenvalues = eigvals(matrix)
            # save all the eigenvalues to their lists
            for b in 1:length(unitcell.basis)
                segments_total[s][b][i] = eigenvalues[b]
            end
        end
    end

    # generate a new band structure object
    bandstructure = Bandstructure(path, segments_total)

    # return the band structure
    return bandstructure
end
function getBandStructure(
                matrixFunction::Function,
                path::Path;
                resolution::Int64=-1
            )

    # maybe modify the path resolution
    if resolution > 0
        setTotalResolution!(path, resolution)
    end

    # build up the list of all bands of all segments
    segments_total = Array{Array{Float64,1},1}[]
    # get the number of bands
    N_bands = size(matrixFunction(path.points[1]),1)
    # iterate over all path segments and push empty lists into the segments list
    for i in 1:length(path.segment_resolution)
        # build an empty band structure for this segment
        segment = Array{Float64, 1}[]
        for b in 1:N_bands
            push!(segment, zeros(Float64, path.segment_resolution[i]))
        end
        # push the segment band structure into the complete segment list
        push!(segments_total, segment)
    end


    # iterate over all segments
    for s in 1:length(path.segment_resolution)
        # get the grid in between two points
        segment_resolution = path.segment_resolution[s]
        # get all multipliers of k vectors (i.e. all alpha in (1-alpha)*k_1 + alpha*k_2)
        multipliers = linspace(0, 1, segment_resolution)
        # get the local start and end point of the segment
        k1 = path.points[s]
        k2 = path.points[s+1]
        # calculate all energies
        for i in 1:segment_resolution
            # get the current k
            k = (k2 .* multipliers[i]) .+ (k1 .* (1-multipliers[i]))
            # get the interaction matrix for this k
            matrix = matrixFunction(k)
            # diagonalize the matrix
            eigenvalues = eigvals(matrix)
            # save all the eigenvalues to their lists
            for b in 1:length(eigenvalues)
                segments_total[s][b][i] = eigenvalues[b]
            end
        end
    end

    # generate a new band structure object
    bandstructure = Bandstructure(path, segments_total)

    # return the band structure
    return bandstructure
end
export getBandStructure





################################################################################
#
#   BAND STRUCTURE PLOTTING
#
################################################################################


"""
    plotBandstructure(
            bandstructure::Bandstructure
         [; limits_energy="AUTO",
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,4),
            showPlot::Bool=true,
            save_filename::String="NONE" ]
            )

    plotBandstructure(
            unitcell::Unitcell,
            path::Path
         [; resolution::Int64=-1,
            enforce_hermitian::Bool=false,
            ... ]
        )

    plotBandstructure(
            matrixFunction::Function,
            path::Path
         [; resolution::Int64=-1,
            ... ]
        )

Plots the band struture of a passed `Bandstructure` object along some its path and returns the plot as a `PyPlot.Figure` object.
Alternatively, one can pass a `Unitcell` (or a matrix function) and `Path` to calculate the band structure which is plotted.

Additional options include setting plotting related options of `PyPlot` as well as determining if the plot is saved or shown.


# Examples

```julia-repl
julia> plotBandstructure(unitcell, path)
PyPlot.Figure(...)

julia> plotBandstructure(unitcell, path, showPlot=false)
PyPlot.Figure(...)

julia> plotBandstructure(unitcell, save_filename="myplot.pdf")
PyPlot.Figure(...)

julia> plotBandstructure(bandstructure)
PyPlot.Figure(...)
```
"""
function plotBandstructure(
            bandstructure::Bandstructure;
            limits_energy="AUTO",
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,4),
            showPlot::Bool=true,
            save_filename::String="NONE"
        )

    ###########################
    #   INITIAL SETTINGS
    ###########################

    # get the path from the bandstructure
    path = bandstructure.path

    # configure plot environment
    rc("font", family="serif")

    # create a new figure
    fig = figure(figsize=figsize)




    ###########################
    #   PLOT BANDS
    ###########################

    # plot the band structure
    for s in 1:length(bandstructure.bands)
        # plot the segment
        for b in bandstructure.bands[s]
            # plot the band
            plot(
                collect(1:path.segment_resolution[s]) .+ sum(path.segment_resolution[1:s-1]), # xvalues
                b,
                "-$(plot_color)"
            )
        end
    end



    ###########################
    #   SET ALL TICKS (POINTS)
    ###########################

    # get the current axis
    ax = gca()
    axx = ax[:get_xaxis]()
    # compile tick positions and labels
    point_pos = Int64[]
    push!(point_pos, 1)
    for l in 1:length(path.segment_resolution)
        push!(point_pos, sum(path.segment_resolution[1:l]))
    end
    point_labels = String[path.point_names[i] for i in 1:length(path.points)]
    # configure tick labels
    xticks(point_pos, point_labels)
    # configure ticks
    axx[:set_tick_params](which="both", direction="out")
    axx[:set_tick_params](which="top", color="none")
    axy = ax[:get_yaxis]()
    axy[:set_tick_params](which="both", direction="out")

    # plot vertical lines for each point
    for p in point_pos
        axvline(p,color=[0.6, 0.6, 0.6], linestyle="--")
    end


    ###########################
    #   CONFIGURE AXIS & TITLE
    ###########################

    # label the axis
    xlabel("momentum")
    ylabel("energy")

    # energy limits
    # check if specific boundaries are desired
    if !(limits_energy == "AUTO")
        ylim(limits_energy[1], limits_energy[2])
    end

    # momentum limits (x axis)
    xlim(0, maximum(point_pos)+1)

    # set the title
    if plot_title == "AUTO"
        # set the title to an automatically generated title
        title("energy spectrum along path $(getPathString(path))")
    elseif plot_title == ""
        # do nothing title related
    else
        # set the title to the given title
        title(plot_title)
    end





    ###########################
    #   FINISH THE PLOT
    ###########################

    # tighten the layout
    tight_layout()

    # save the plot
    if save_filename != "NONE"
        # make sure the directory exists
        if contains(save_filename, "/")
    		# get the containing folder
    		folder = save_filename[1:findlast(save_filename, '/')]
    		# build the path to that folder
    		mkpath(folder)
    	end
        # save the plot
        savefig(save_filename)
    end

    # maybe show the plot
    if showPlot
        show()
    end

    # return the figure object
    return fig
end
function plotBandstructure(
            unitcell::Unitcell,
            path::Path;
            resolution::Int64=-1,
            enforce_hermitian::Bool=false,
            limits_energy="AUTO",
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,4),
            showPlot::Bool=true,
            save_filename::String="NONE"
        )
    # calculate the bandstructure
    bandstructure = getBandStructureAlongPath(unitcell, path, resolution=resolution, enforce_hermitian=enforce_hermitian)
    # call the respective function
    return plotBandstructure(
                bandstructure;
                limits_energy=limits_energy,
                plot_title=plot_title,
                plot_color=plot_color,
                figsize=figsize,
                showPlot=showPlot,
                save_filename=save_filename
            )
end
function plotBandstructure(
            matrixFunction::Function,
            path::Path;
            resolution::Int64=-1,
            limits_energy="AUTO",
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,4),
            showPlot::Bool=true,
            save_filename::String="NONE"
        )
    # calculate the bandstructure
    bandstructure = getBandStructureAlongPath(matrixFunction, path, resolution=resolution)
    # call the respective function
    return plotBandstructure(
                bandstructure;
                limits_energy=limits_energy,
                plot_title=plot_title,
                plot_color=plot_color,
                figsize=figsize,
                showPlot=showPlot,
                save_filename=save_filename
            )
end
export plotBandstructure
