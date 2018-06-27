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
