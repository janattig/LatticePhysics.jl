################################################################################
#
#   METHODS FOR LUTTINGER TISZA CALCULATION
#
#   STRUCTURE OF THE FILE
#
#   1) TYPE LTBANDSTRUCTURE
#       - type definition
#       - printInfo function
#
#   2) CALCULATION OF LT BAND STRUCTURES OF UNTICELL OBJECTS
#       - LT constraint and deviation functions
#       - spin interaction matrices
#       - calculation of band structures
#
#   3) TODO PLOTTING OF LT BAND STRUCTURES
#       - TODO plotting of Bandstructure objects
#       - TODO plotting of bandstructures of unitcells along paths
#
#   4) TODO CALCULATION OF LT GROUND STATES (k space manifold)
#
#   5) TODO PLOTTING OF LT GROUND STATES (k space manifold)
#       - TODO plotting from points
#       - TODO plotting from unitcell
#
################################################################################










################################################################################
#
#   TYPE LTBANDSTRUCTURE
#       - type definition
#       - printInfo function
#
################################################################################
struct LTBandstructure

    # the path along which the band structure is calcualted
    path::Path

    # bands for each segment
    # bands[i] gives all bands of segment i
    # bands[i][j] gives all energy values for band j of segment i
    # bands[i][j][k] gives the energy value at kpoint index k of band j in segment i
    bands::Array{Array{Array{Float64, 1}, 1}, 1}

    # constraint value for all bands
    # value is the minimum of all sum(|s_i - 1.0|^2) for constructed s_i
    constraint_value::Array{Array{Array{Float64, 1}, 1}, 1}

    # ONLY DEFAULT CONSTRUCTOR

end


# export the type
export LTBandstructure


# INFORMATION FUNCTION
function printInfo(bandstructure::LTBandstructure; constraint::Float64=1e-6)
    # print the header
    println("Bandstructure (LT), constraint satisfied for var(length) < $(constraint)")
    print("\t$(bandstructure.path.point_names[1])")
    for p in bandstructure.path.point_names[2:end]
        print("\t->-\t$(p)")
    end
    println("")
    # print each band
    for b in 1:length(bandstructure.bands[1])
        # first the point
        print("$(b))\t|")
        # then all segments with closing points
        for s in 1:length(bandstructure.bands)
            # segment s starts
            print("\t")
            # calculate the number of constraint fullfilling values
            print("$(round(100.0*sum([c>constraint ? 0 : 1 for c in bandstructure.constraint_value[s][b]])/length(bandstructure.constraint_value[s][b]),2))%")
            # print the next point
            print("\t|")
        end
        # print a new line
        println("")
    end
end















# function to calculate the deviation from best result
function deviation(spin_eigenvectors::Array{Array{Complex{Float64},1},1}, spin_dimension::Int64, alpha::Array{Float64,1})
    # build up the global spin vector
    spin_vector = zeros(length(spin_eigenvectors[1]))
    for s in 1:length(spin_eigenvectors)
        spin_vector = spin_vector .+ alpha[s].*spin_eigenvectors[s]
    end
    # find out the individual lengths of spins
    spin_lengths = spin_vector .* conj.(spin_vector)
    spin_lengths = [sum(spin_lengths[s:s+spin_dimension-1]) for s in 1:spin_dimension:length(spin_lengths)]
    # find out the deviation from unity
    dl = sum(abs.((spin_lengths .- 1)))
    return dl
end


# Definition of constraint
function getLTConstraint(spin_eigenvectors::Array{Array{Complex{Float64},1},1}, spin_dimension::Int64)
    # detrmine what to do based on the number of eigenvectors
    if length(spin_eigenvectors) == 1
        # find out the individual lengths of spins
        spin_lengths = spin_eigenvectors[1] .* conj.(spin_eigenvectors[1])
        spin_lengths = [sum(spin_lengths[s:s+spin_dimension-1]) for s in 1:spin_dimension:length(spin_lengths)]
        # find out the deviation from unity
        dl = sum(abs.((spin_lengths .- mean(spin_lengths))))
    else
        # optimize the function
        return Optim.minimum(Optim.optimize(x -> deviation(spin_eigenvectors, spin_dimension, x), ones(length(spin_eigenvectors))))
    end
end





# Function to create a bond strength matrix
function getBondInteractionMatrixHeisenbergKitaev(connection::Array{Any,1})
    # new 3x3 matrix
    bond_matrix = zeros(3,3)
    # get the bond strength
    strength = connection[3]
    # check what type the bond is
    if typeof(strength) == String
        if strength == "J1"
            bond_matrix[1,1] = 1.0
            bond_matrix[2,2] = 1.0
            bond_matrix[3,3] = 1.0
        elseif strength == "J2"
            bond_matrix[1,1] = 1.0
            bond_matrix[2,2] = 1.0
            bond_matrix[3,3] = 1.0
        elseif strength == "Jx" || strength == "tx"
            bond_matrix[1,1] = 1.0
        elseif strength == "Jy" || strength == "ty"
            bond_matrix[2,2] = 1.0
        elseif strength == "Jz" || strength == "tz"
            bond_matrix[3,3] = 1.0
        end
    else
        bond_matrix[1,1] = strength
        bond_matrix[2,2] = strength
        bond_matrix[3,3] = strength
    end
    # return the matrix
    return bond_matrix
end
function getBondInteractionMatrixHeisenberg(connection::Array{Any,1})
    # new 3x3 matrix
    bond_matrix = zeros(1,1)
    # get the bond strength
    strength = connection[3]
    # check what type the bond is
    if typeof(strength) == String
        if strength == "J1"
            bond_matrix[1,1] = 1.0
        elseif strength == "J2"
            bond_matrix[1,1] = 1.0
        end
    else
        bond_matrix[1,1] = strength
    end
    # return the matrix
    return bond_matrix
end

# Function to produce interaction matrices
function getSpinInteractionMatrixKSpace(unitcell::Unitcell, k_vector::Array{Float64,1}, bondInteractionMatrix::Function)
    # get the spin dimension
    spin_dimension = size(bondInteractionMatrix(unitcell.connections[1]), 1)
    # create a new matrix
    matrix = zeros(Complex, spin_dimension*length(unitcell.basis), spin_dimension*length(unitcell.basis))
    # iterate over all connections
    for c in unitcell.connections
        # get the indices
        index_from  = Int(c[1])
        index_to    = Int(c[2])
        strength    = c[3]
        wrap        = c[4]
        # get the difference vector
        pos_delta   = unitcell.basis[index_to] .- unitcell.basis[index_from]
        if size(unitcell.lattice_vectors,1) > 0
            for pair in zip(wrap, unitcell.lattice_vectors)
                pos_delta .+= pair[1].*pair[2]
            end
        end
        # get the bond interaction matrix
        bond_interaction_matrix = bondInteractionMatrix(c)
        # add to the complete matrix twice
        for i in 1:spin_dimension
        for j in 1:spin_dimension
            matrix[(index_from-1)*spin_dimension + i, (index_to-1)*spin_dimension + j] += 0.5 * bond_interaction_matrix[i,j] * exp(-im * sum(pos_delta .* k_vector))
            matrix[(index_to-1)*spin_dimension + j, (index_from-1)*spin_dimension + i] += 0.5 * bond_interaction_matrix[i,j] * exp( im * sum(pos_delta .* k_vector))
        end
        end
    end
    # return the matrix
    return matrix
end


################################################################################
#
#   BAND STRUCTURE CALCULATION
#
################################################################################
"""
    getLTBandStructure(
                unitcell::Unitcell,
                path::Path,
                bondInteractionMatrix::Function
             [; resolution::Int64=-1,
                enforce_hermitian::Bool=false ]
            )

Calculates the Luttinger Tisza band struture of a `Unitcell` object
along some path given by a `Path` object and returns a `LTBandstructure` object.

Note 1: Optionally a bond interaction matrix can be given as a function
that constructs a matrix for a passed bond argument. If this function is not passed,
the default function will be used.

Note 2: The bond interaction matrix specifies the dimension of interacting spins.


# Examples

```julia-repl
julia> bandstructure = getLTBandStructure(unitcell, path)
LatticePhysics.LTBandstructure(...)

julia> bandstructure = getLTBandStructure(unitcell, path, resolution=1000)
LatticePhysics.LTBandstructure(...)
```
"""
function getLTBandStructure(
                unitcell::Unitcell,
                path::Path,
                bondInteractionMatrix::Function = getBondInteractionMatrixHeisenberg;
                resolution::Int64=-1,
                enforce_hermitian::Bool=false,
                epsilon_degenerate::Float64=1e-6
            )

    # maybe modify the path resolution
    if resolution > 0
        setTotalResolution!(path, resolution)
    end

    # build up the list of all bands of all segments (and all constraints)
    segments_total    = Array{Array{Float64,1},1}[]
    constraints_total = Array{Array{Float64,1},1}[]

    # get the spin dimension
    spin_dimension = size(bondInteractionMatrix(unitcell.connections[1]), 1)
    # iterate over all path segments and push empty lists into the segments list
    for i in 1:length(path.segment_resolution)
        # build an empty band structure for this segment
        segment     = Array{Float64, 1}[]
        constraints = Array{Float64, 1}[]
        for b in 1:length(unitcell.basis)*spin_dimension
            push!(segment,     zeros(Float64, path.segment_resolution[i]))
            push!(constraints, zeros(Float64, path.segment_resolution[i]))
        end
        # push the segment band structure into the complete segment list
        push!(segments_total,    segment)
        push!(constraints_total, constraints)
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
            matrix = getSpinInteractionMatrixKSpace(unitcell, k, bondInteractionMatrix)
            # diagonalize the matrix
            eigenfactorization = eigfact(matrix)
            eigenvalues  = eigenfactorization[:values]
            eigenvectors = eigenfactorization[:vectors]
            # save all the eigenvalues to their lists
            for b in 1:length(eigenvalues)
                segments_total[s][b][i] = eigenvalues[b]
            end
            # compute the constraint, first find out what bands are degenerate
            # list of bands they are degenerate with
            degenerate = zeros(Int64,   length(eigenvalues)) .- 1
            treat      = zeros(Int64,   length(eigenvalues))
            for b in 2:length(eigenvalues)
                # treat the current band
                treat[b] = 1
                # check if degenerate
                if eigenvalues[b] - epsilon_degenerate <= eigenvalues[b-1]
                    # band b is degenerate with band b-1
                    degenerate[b-1] = b
                    # dont treat the band as it is treated before
                    treat[b] = 0
                end
            end
            # iterate over all bands
            for b in 1:length(eigenvalues)
                # if not treated, continue
                if treat[b] == 0
                    continue
                end
                # if treated, compile list of all bands
                degenerate_bands = Int64[b]
                while degenerate[degenerate_bands[end]] != -1
                    push!(degenerate_bands, degenerate[degenerate_bands[end]])
                end
                # get the LT constraint
                LT_constraint = getLTConstraint([eigenvectors[:,j] for j in degenerate_bands], spin_dimension)
                # save the constraint
                for d in degenerate_bands
                    constraints_total[s][d][i] = LT_constraint
                end
            end
        end
    end

    # generate a new LT band structure object
    bandstructure = LTBandstructure(path, segments_total, constraints_total)

    # return the LT band structure
    return bandstructure
end
