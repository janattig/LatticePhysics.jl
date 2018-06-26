################################################################################
#
#   METHODS FOR CONSTRUCTION INTERACTION MATRICES FOR LATTICES
#
#   STRUCTURE OF THE FILE
#
#   1) INTERACTION MATRICES IN REAL SPACE
#
#   2) INTERACTION MATRICES IN MOMENTUM SPACE
#
#   NOTE: No Majorana fermions so far on the level of matrices, because of some
#   problems with the various gauges that can be tuned. Want: sublattice
#   idenfitication in bipartite systems.
#
################################################################################



#-----------------------------------------------------------------------------------------------------------------------------
#
#   Interaction matrix in REAL space
#
#   Parameters are
#   - lattice: The complete lattice of which the interaction matrix should be constructed
#   - enforce_hermitian (optional): if the matrix should be made hermitian by 0.5*(A + A_dagger)
#
#-----------------------------------------------------------------------------------------------------------------------------
"""
    getInteractionMatrixRealSpace(unitcell::Unitcell  [ ; enforce_hermitian::Bool=false ])
    getInteractionMatrixRealSpace(lattice::Lattice    [ ; enforce_hermitian::Bool=false ])

Constructs the interaction matrix (in *position* space) of a given `Unitcell` or `Lattice` object.
The matrix is a NxN matrix where N is the number of sites in the given object.
Entries (i,j) contain the interactions between sites i and j.

The matrix is of type `Array{Complex,2}`, i.e. it has complex entries to account for eventual complex hopping amplitudes.

Furthermore, the matrix can be constructed hermitian by using `enforce_hermitian=true` in which case
the interactions between i and j are not only added to element (i,j) of the matrix but also to element (j,i) as conjugated.



# Examples

```julia-repl
julia> getInteractionMatrixRealSpace(lattice)
2×2 Array{Complex,2}:
...

julia> getInteractionMatrixRealSpace(unitcell, enforce_hermitian=true)
2×2 Array{Complex,2}:
...

```
"""
function getInteractionMatrixRealSpace(lattice::Lattice; enforce_hermitian::Bool=false)
    # generate a new matrix
    matrix = zeros(Complex, length(lattice.positions),length(lattice.positions))
    # decide if hermitian or not hermitian enforced
    if enforce_hermitian
        # iterate over all connections and add strength to two entries
        for c in lattice.connections
            # get the indices
            index_from  = Int(c[1])
            index_to    = Int(c[2])
            strength    = c[3]
            # add to the matrix twice
            matrix[index_from, index_to] += 0.5 * strength
            matrix[index_to, index_from] += 0.5 * conj(strength)
        end
    else
        # iterate over all connections and add strength to one entry
        for c in lattice.connections
            # get the indices
            index_from  = Int(c[1])
            index_to    = Int(c[2])
            strength    = c[3]
            # add to the matrix once
            matrix[index_from, index_to] += strength
        end
    end
    # return the matrix
    return matrix
end
function getInteractionMatrixRealSpace(unitcell::Unitcell; enforce_hermitian::Bool=false)
    # generate a new matrix
    matrix = zeros(Complex, length(unitcell.basis),length(unitcell.basis))
    # decide if hermitian or not hermitian enforced
    if enforce_hermitian
        # iterate over all connections and add strength to two entries
        for c in unitcell.connections
            # get the indices
            index_from  = Int(c[1])
            index_to    = Int(c[2])
            strength    = c[3]
            # add to the matrix twice
            matrix[index_from, index_to] += 0.5 * strength
            matrix[index_to, index_from] += 0.5 * conj(strength)
        end
    else
        # iterate over all connections and add strength to one entry
        for c in unitcell.connections
            # get the indices
            index_from  = Int(c[1])
            index_to    = Int(c[2])
            strength    = c[3]
            # add to the matrix once
            matrix[index_from, index_to] += strength
        end
    end
    # return the matrix
    return matrix
end
export getInteractionMatrixRealSpace



#-----------------------------------------------------------------------------------------------------------------------------
#
#   Interaction matrix in RECIPROCAL (K) space
#
#   Parameters are
#   - lattice: The complete lattice of which the interaction matrix should be constructed
#   - k_vector: The reciprocal point k at which the matrix is constructed
#   - enforce_hermitian (optional): if the matrix should be made hermitian by 0.5*(A + A_dagger)
#
#-----------------------------------------------------------------------------------------------------------------------------
"""
    getInteractionMatrixKSpace(unitcell::Unitcell, k_vector::Array{Float64,1} [ ; enforce_hermitian::Bool=false ])
    getInteractionMatrixKSpace(lattice::Lattice,   k_vector::Array{Float64,1} [ ; enforce_hermitian::Bool=false ])

Constructs the interaction matrix (in *momentum* space at point `k_vector`) of a given `Unitcell` or `Lattice` object.
The matrix is a NxN matrix where N is the number of sites in the given object.
Entries (i,j) contain the interactions between sites i and j as well as the phase factor exp(i k*delta).

The matrix is of type `Array{Complex,2}`, i.e. it has complex entries to account for the phase factor.

Furthermore, the matrix can be constructed hermitian by using `enforce_hermitian=true` in which case
the interactions between i and j are not only added to element (i,j) of the matrix but also to element (j,i) as conjugated.



# Examples

```julia-repl
julia> getInteractionMatrixKSpace(lattice, k)
2×2 Array{Complex,2}:
...

julia> getInteractionMatrixKSpace(lattice, [pi/2.0, 0.0], enforce_hermitian=true)
2×2 Array{Complex,2}:
...

```
"""
function getInteractionMatrixKSpace(lattice::Lattice, k_vector::Array{Float64,1}; enforce_hermitian::Bool=false)
    # generate a new matrix
    matrix = zeros(Complex, size(lattice.positions,1),size(lattice.positions,1))
    # decide if hermitian or not hermitian enforced
    if enforce_hermitian
        # iterate over all connections
        for c in lattice.connections
            # get the indices
            index_from  = Int(c[1])
            index_to    = Int(c[2])
            strength    = c[3]
            wrap        = c[4]
            # get the difference vector
            pos_delta   = lattice.positions[index_to] .- lattice.positions[index_from]
            if size(lattice.lattice_vectors,1) > 0
                for pair in zip(wrap, lattice.lattice_vectors)
                    pos_delta .+= pair[1].*pair[2]
                end
            end
            # add to the matrix twice
            matrix[index_from, index_to] += 0.5 * strength       * exp(-im * sum(pos_delta .* k_vector))
            matrix[index_to, index_from] += 0.5 * conj(strength) * exp(-im * sum(pos_delta .* k_vector))
        end
    else
        # iterate over all connections
        for c in lattice.connections
            # get the indices
            index_from  = Int(c[1])
            index_to    = Int(c[2])
            strength    = c[3]
            wrap        = c[4]
            # get the difference vector
            pos_delta   = lattice.positions[index_to] .- lattice.positions[index_from]
            if size(lattice.lattice_vectors,1) > 0
                for pair in zip(wrap, lattice.lattice_vectors)
                    pos_delta .+= pair[1].*pair[2]
                end
            end
            # add to the matrix once
            matrix[index_from, index_to] += strength * exp(-im * sum(pos_delta .* k_vector))
        end
    end
    # return the matrix
    return matrix
end
function getInteractionMatrixKSpace(unitcell::Unitcell, k_vector::Array{Float64,1}; enforce_hermitian::Bool=false)
    # generate a new matrix
    matrix = zeros(Complex, size(unitcell.basis,1),size(unitcell.basis,1))
    # decide if hermitian or not hermitian enforced
    if enforce_hermitian
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
            # add to the matrix twice
            matrix[index_from, index_to] += 0.5 * strength       * exp(-im * sum(pos_delta .* k_vector))
            matrix[index_to, index_from] += 0.5 * conj(strength) * exp(-im * sum(pos_delta .* k_vector))
        end
    else
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
            # add to the matrix once
            matrix[index_from, index_to] += strength * exp(-im * sum(pos_delta .* k_vector))
        end
    end
    # return the matrix
    return matrix
end
export getInteractionMatrixKSpace
