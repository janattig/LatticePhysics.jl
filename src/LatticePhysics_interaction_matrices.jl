#-----------------------------------------------------------------------------------------------------------------------------
#
#   METHODS FOR CONSTRUCTION INTERACTION MATRICES FOR LATTICES
#
#-----------------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------------
#
#   Interaction matrix in REAL space
#
#   Parameters are
#   - lattice: The complete lattice of which the interaction matrix should be constructed
#   - enforce_hermitian (optional): if the matrix should be made hermitian by 0.5*(A + A_dagger)
#
#-----------------------------------------------------------------------------------------------------------------------------
function getInteractionMatrixRealSpace(lattice::Lattice; enforce_hermitian=false)
    # generate a new matrix
    matrix = zeros(size(lattice.positions,1),size(lattice.positions,1))
    # iterate over all connections
    for c in lattice.connections
        # get the indices
        index_from  = Int(c[1])
        index_to    = Int(c[2])
        strength    = c[3]
        # just add to the matrix
        matrix[index_from, index_to] += strength
    end
    # eventually ensure the hermitian nature of the matrix
    if enforce_hermitian
        matrix = 0.5*(matrix .+ transpose(conj(matrix)))
    end
    # return the matrix
    return matrix
end
function getInteractionMatrixRealSpace(unitcell::Unitcell; enforce_hermitian=false)
    # generate a new matrix
    matrix = zeros(size(unitcell.basis,1),size(unitcell.basis,1))
    # iterate over all connections
    for c in unitcell.connections
        # get the indices
        index_from  = Int(c[1])
        index_to    = Int(c[2])
        strength    = c[3]
        # just add to the matrix
        matrix[index_from, index_to] += strength
    end
    # eventually ensure the hermitian nature of the matrix
    if enforce_hermitian
        matrix = 0.5*(matrix .+ transpose(conj(matrix)))
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
function getInteractionMatrixKSpace(lattice::Lattice, k_vector::Array{Float64,1}; enforce_hermitian=false, majorana=false)
    # generate a new matrix
    matrix = zeros(Complex, size(lattice.positions,1),size(lattice.positions,1))
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
        # take majorana fermionic statistics into account
        majorana_factor = 1
        if majorana && index_from < index_to
            majorana_factor = 1.0im
        elseif majorana && index_from > index_to
            majorana_factor = -1.0im
        end
        # just add to the matrix
        matrix[index_from, index_to] += strength * exp(-sum(pos_delta.*k_vector) * im) * majorana_factor
    end
    # eventually ensure the hermitian nature of the matrix
    if enforce_hermitian
        matrix = 0.5*(matrix .+ transpose(conj(matrix)))
    end
    # return the matrix
    return matrix
end
function getInteractionMatrixKSpace(unitcell::Unitcell, k_vector::Array{Float64,1}; enforce_hermitian=false, majorana=false)
    # generate a new matrix
    matrix = zeros(Complex, size(unitcell.basis,1),size(unitcell.basis,1))
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
        # take majorana fermionic statistics into account
        majorana_factor = 1
        if majorana && index_from < index_to
            majorana_factor = 1.0im
        elseif majorana && index_from > index_to
            majorana_factor = -1.0im
        end
        # just add to the matrix
        matrix[index_from, index_to] += strength * exp(-sum(pos_delta.*k_vector) * im) * majorana_factor
    end
    # eventually ensure the hermitian nature of the matrix
    if enforce_hermitian
        matrix = 0.5*(matrix .+ transpose(conj(matrix)))
    end
    # return the matrix
    return matrix
end
export getInteractionMatrixKSpace

#-----------------------------------------------------------------------------------------------------------------------------
#
#   Interaction matrix INFORMATION
#   (returns two matrices that indicate the number of periodic and non-periodic connections between respective sites)
#
#   Parameters are
#   - lattice: The complete lattice of which the interaction matrix should be constructed
#   or
#   - unitcell: The unitcell of which the interaction matrix should be constructed
#
#-----------------------------------------------------------------------------------------------------------------------------
function getInteractionMatrixInformation(lattice::Lattice)
    # matrices indicating the number of connections
    con_periodic =  zeros(Int64, size(lattice.positions,1),size(lattice.positions,1))
    con_inside   =  zeros(Int64, size(lattice.positions,1),size(lattice.positions,1))
    # go through all connections of the lattice
    for c in lattice.connections
        # get the indices
        index_from  = Int(c[1])
        index_to    = Int(c[2])
        strength    = c[3]
        wrap        = c[4]
        # check if inside or periodic and add to the counter
        if sum([abs(el) for el in wrap]) == 0
            con_inside[index_from, index_to] = con_inside[index_from, index_to] + 1
        else
            con_periodic[index_from, index_to] = con_periodic[index_from, index_to] + 1
        end
    end
    # return the matrices
    return con_inside, con_periodic
end
function getInteractionMatrixInformation(unitcell::Unitcell)
    # matrices indicating the number of connections
    con_periodic =  zeros(Int64, size(unitcell.basis,1),size(unitcell.basis,1))
    con_inside   =  zeros(Int64, size(unitcell.basis,1),size(unitcell.basis,1))
    # go through all connections of the lattice
    for c in unitcell.connections
        # get the indices
        index_from  = Int(c[1])
        index_to    = Int(c[2])
        strength    = c[3]
        wrap        = c[4]
        # check if inside or periodic and add to the counter
        if sum([abs(el) for el in wrap]) == 0
            con_inside[index_from, index_to] = con_inside[index_from, index_to] + 1
        else
            con_periodic[index_from, index_to] = con_periodic[index_from, index_to] + 1
        end
    end
    # return the matrices
    return con_inside, con_periodic
end


# print the interaction matrix information in compressed LaTeX format
function printInteractionMatrixInformation(connections_inside, connections_periodic)
    # get the dimension of the matrices
    len = size(connections_periodic,1)
    # starting lines
    println("\\begin{equation*}")
    println("   \\mathbf{H} = \\bordermatrix{")
    print("       ")
    for i in 1:len
        print("& $(i) ")
    end
    println("\\cr")
    # iterate over all lines of the matrix
    for j in 1:len
        # start line by printing the current index
        print("       $(j) ")
        # go through all respective numbers
        for i in 1:len
            # check if periodic or internal connections
            if connections_inside[i,j] == 0 && connections_periodic[i,j] == 0
                print("& 0 ")
            elseif connections_inside[i,j] > 0 && connections_periodic[i,j] == 0
                #print("& {\\color{blue} 0},{\\color{red} $(connections_inside[i,j])} ")
                print("& {\\color{blue} $(connections_inside[i,j])} ")
            elseif connections_inside[i,j] == 0 && connections_periodic[i,j] > 0
                #print("& {\\color{blue} $(connections_periodic[i,j])},{\\color{red} 0} ")
                print("& {\\color{red} $(connections_periodic[i,j])} ")
            else
                print("& {\\color{red} $(connections_periodic[i,j])},{\\color{blue} $(connections_inside[i,j])} ")
            end
        end
        # end the line by printing a linebreak
        println("\\cr")
    end
    # end lines
    println("   }")
    println("\\end{equation*}")
end
function printInteractionMatrixInformation(lattice::Lattice)
    println("Interaction matrix for lattice \"$(lattice.filename)\"")
    connections_inside, connections_periodic = getInteractionMatrixInformation(lattice)
    printInteractionMatrixInformation(connections_inside, connections_periodic)
end
function printInteractionMatrixInformation(unitcell::Unitcell)
    println("Interaction matrix for unitcell \"$(unitcell.filename)\"")
    connections_inside, connections_periodic = getInteractionMatrixInformation(unitcell)
    printInteractionMatrixInformation(connections_inside, connections_periodic)
end
export printInteractionMatrixInformation
