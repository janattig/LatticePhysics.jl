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










################################################################################
#
#   BAND STRUCTURE CALCULATION
#
################################################################################
"""
    getBandStructureAlongPath(
                unitcell::Unitcell,
                path::Path
             [; resolution::Int64=-1,
                enforce_hermitian::Bool=false ]
            )

Calculates the band struture of a `Unitcell` object along some path given by a `Path` object
and returns a `Bandstructure` object.


# Examples

```julia-repl
julia> bandstructure = getBandStructureAlongPath(unitcell, path)
LatticePhysics.Bandstructure(...)

julia> bandstructure = getBandStructureAlongPath(unitcell, path, resolution=1000)
LatticePhysics.Bandstructure(...)
```
"""
function getBandStructureAlongPath(
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






################################################################################
#
#   BAND STRUCTURE PLOTTING
#
################################################################################

# uses PyPlot
using PyPlot


function plotBandstructure(
            bandstructure::Bandstructure;
            limits_energy="AUTO",
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,4),
            showPlot::Bool=true
        )

    # configure plot environment
    rc("font", family="serif")

    # create a new figure
    fig = figure(figsize=figsize)

    # set the title
    if plot_title == "AUTO"
        # set the title to an automatically generated title
        title("energy spectrum along path $(getPathString())")
    elseif plot_title == ""
        # do nothing title related
    else
        # set the title to the given title
        title(plot_title)
    end
    for l in hlines[1:end-1]
    axvline(l,color=[0.6, 0.6, 0.6], linestyle="--")
    end
    xlabel("momentum")
    ylabel("energy")
    for b in bandstructure
    plot(collect(1:resolution_actual), b, "-$(plot_color)")
    end
    ax = gca()
    axx = ax[:get_xaxis]()
    xtpos = []
    push!(xtpos, 0)
    for h in hlines
    push!(xtpos, h)
    end
    xtlabs = [p[1] for p in path]
    xticks(xtpos, xtlabs)
    #axx[:set_ticks]([])
    axx[:set_tick_params](which="both", direction="out")
    axx[:set_tick_params](which="top", color="none")
    axy = ax[:get_yaxis]()
    axy[:set_tick_params](which="both", direction="out")
    # check if specific boundaries are desired
    if !(limits_energy == "AUTO")
    ylim(limits_energy[1], limits_energy[2])
    end
    # tighten the layout
    tight_layout()
    # save the plot
    figurename = split(lattice.filename, FOLDER_SPECTRA[end])[end]
    figurename1 = "$(FOLDER_SPECTRA)bandstructure_path_$(figurename[1:end-4]).pdf"
    figurename2 = "$(FOLDER_SPECTRA)bandstructure_path_$(figurename[1:end-4]).png"
    #buildFolderSpectra()
    savefig(figurename1)
    savefig(figurename2)
    # maybe show the plot
    if showPlot
        show()
    end

    # return the figure object
    return fig
end
