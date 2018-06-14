################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT LATTICE CONSTRUCTION RELATED STUFF
#   AS WELL AS DIFFERNT LATTICE MODIFICATION STUFF
#
#   STRUCTURE OF THE FILE
#
################################################################################





#-----------------------------------------------------------------------------------------------------------------------------
#
#   INDIVIDUAL FUNCTIONS FOR LATTICES COMPOSED FROM UNITCELLS
#   BUILDING FUNCTIONS FOR NORMAL PLACEMENT OF CELLS
#
#   Parameters that have to be passed:
#   - Unitcell object from which lattice is build
#   - Array with integers indicating the extent of the lattice along this lattice vector
#
#   Functions will start with most specific and will be generalized later on
#
#-----------------------------------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------------------------------
#
#   Construction of PERIODIC Lattices (i.e. lattices, where all directions have periodic boundary conditions)
#
#-----------------------------------------------------------------------------------------------------------------------------

# FOR 2D AND 3D
function getLatticePeriodic2D(unitcell::Unitcell, repetition_array::Array{Int64}; save::Bool=false, load::Bool=false)

    # extract the cardinal directions of the lattice from the array
    N_a1 = abs(repetition_array[1])
    N_a2 = abs(repetition_array[2])

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_periodic_$(N_a1)_$(N_a2).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "FOLDER_LATTICES$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
        filename = replace(filename, ".jld", "_periodic_$(N_a1)_$(N_a2).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    end

    # check if can be loaded
    if isfile(filename) && load
        return Lattice(filename)
    end

    # build the lattice

    # load the data from the unit cell
    uc_basis            = unitcell.basis
    uc_connections      = unitcell.connections
    uc_lattice_vectors  = unitcell.lattice_vectors

    # GENERATE NEW POSITIONS
	positions = Array[]
    positions_indices = Int64[]

	# define the index function to get the index of an element in the array
	function index(i,j,alpha)
		return size(uc_basis,1)*((i-1)*N_a2+j-1) + alpha
	end

    # define positions array to be filled
	for element in 1:N_a1*N_a2*size(uc_basis, 1)
		push!(positions, uc_basis[1])
        push!(positions_indices, 1)
	end
	# set all positions to their correct values
	for i in 1:N_a1
	for j in 1:N_a2
    for a in 1:size(uc_basis, 1)
		positions[index(i,j,a)] = uc_basis[a] + i*uc_lattice_vectors[1] + j*uc_lattice_vectors[2]
        positions_indices[index(i,j,a)] = a
	end
	end
    end

    # GENERATE NEW CONNECTIONS
	connections = Array[]

	# iterate over all unit cells
	for i in 1:N_a1
	for j in 1:N_a2
		# add all connections for unitcell (i,j)
		for connection in uc_connections
			# calculate the index from where the connection goes
			index_from = index(i,j,connection[1])
			# calculate the aimed unitcell
			i_to = i + connection[4][1]
			j_to = j + connection[4][2]
			# check if the connection goes around in a1 direction
            offset_a1 = 0
			while i_to < 1
                i_to += N_a1
                offset_a1 += -1
            end
            while i_to > N_a1
                i_to -= N_a1
                offset_a1 += 1
			end
			# check if the connection goes around in a2 direction
            offset_a2 = 0
			while j_to < 1
                j_to += N_a2
                offset_a2 += -1
            end
            while j_to > N_a2
                j_to -= N_a2
                offset_a2 += 1
			end
			# get the index to where the connection goes
			index_to = index(i_to, j_to, connection[2])
            # generate a new connection
            # format [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            connection_new = [index_from; index_to; connection[3]; (offset_a1, offset_a2)]
			# register as connection
			push!(connections, connection_new)
		end
	end
	end

    # generate new lattice vectors, now spanning the WHOLE lattice
    lattice_vectors = Array[]
    a1 = uc_lattice_vectors[1] .* N_a1
    a2 = uc_lattice_vectors[2] .* N_a2
    push!(lattice_vectors, a1)
    push!(lattice_vectors, a2)

    # save everything to a Lattice object
    lattice = Lattice(
        unitcell,
        [N_a1, N_a2],
        lattice_vectors,
        positions,
        positions_indices,
        connections,
        filename
        )
    # save the lattice object
    if save
        saveLattice(lattice)
    end

    # return the lattice
    return lattice

end
function getLatticePeriodic3D(unitcell::Unitcell, repetition_array::Array{Int64}; save::Bool=false, load::Bool=false)

    # extract the cardinal directions of the lattice from the array
    N_a1 = abs(repetition_array[1])
    N_a2 = abs(repetition_array[2])
    N_a3 = abs(repetition_array[3])

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_periodic_$(N_a1)_$(N_a2)_$(N_a3).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "FOLDER_LATTICES$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
        filename = replace(filename, ".jld", "_periodic_$(N_a1)_$(N_a2)_$(N_a3).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    end

    # check if can be loaded
    if isfile(filename) && load
        return Lattice(filename)
    end

    # build the lattice

    # load the data from the unit cell
    uc_basis            = unitcell.basis
    uc_connections      = unitcell.connections
    uc_lattice_vectors  = unitcell.lattice_vectors


    # GENERATE NEW POSITIONS
	positions = Array[]
    positions_indices = Int64[]

	# define the index function to get the index of an element in the array
	function index(i,j,k,alpha)
		return size(uc_basis,1)*((i-1)*N_a2*N_a3 + (j-1)*N_a3 + (k-1)) + alpha
	end

    # define positions array to be filled
	for element in 1:N_a1*N_a2*N_a3*size(uc_basis, 1)
		push!(positions, uc_basis[1])
        push!(positions_indices, 1)
	end
	# set all positions to their correct values
	for i in 1:N_a1
	for j in 1:N_a2
	for k in 1:N_a3
    for a in 1:size(uc_basis, 1)
		positions[index(i,j,k,a)] = uc_basis[a] .+ i.*uc_lattice_vectors[1] .+ j.*uc_lattice_vectors[2] .+ k.*uc_lattice_vectors[3]
        positions_indices[index(i,j,k,a)] = a
	end
	end
    end
    end

    # GENERATE NEW CONNECTIONS
	connections = Array[]

	# iterate over all unit cells
	for i in 1:N_a1
	for j in 1:N_a2
	for k in 1:N_a3
		# add all connections for unitcell (i,j)
		for connection in uc_connections
			# calculate the index from where the connection goes
			index_from = index(i,j,k,Int(connection[1]))
			# calculate the aimed unitcell
			i_to = Int(i + connection[4][1])
			j_to = Int(j + connection[4][2])
			k_to = Int(k + connection[4][3])
			# check if the connection goes around in a1 direction
            offset_a1 = 0
			while i_to < 1
                i_to += N_a1
                offset_a1 += -1
            end
            while i_to > N_a1
                i_to -= N_a1
                offset_a1 += 1
			end
			# check if the connection goes around in a2 direction
            offset_a2 = 0
			while j_to < 1
                j_to += N_a2
                offset_a2 += -1
            end
            while j_to > N_a2
                j_to -= N_a2
                offset_a2 += 1
			end
			# check if the connection goes around in a3 direction
            offset_a3 = 0
			while k_to < 1
                k_to += N_a3
                offset_a3 += -1
            end
            while k_to > N_a3
                k_to -= N_a3
                offset_a3 += 1
			end
			# get the index to where the connection goes
			index_to = index(i_to, j_to, k_to, Int(connection[2]))
            # generate a new connection
            # format [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            connection_new = [index_from; index_to; connection[3]; (offset_a1, offset_a2, offset_a3)]
			# register as connection
			push!(connections, connection_new)
		end
	end
	end
    end


    # generate new lattice vectors, now spanning the WHOLE lattice
    lattice_vectors = Array[]
    a1 = uc_lattice_vectors[1] .* N_a1
    a2 = uc_lattice_vectors[2] .* N_a2
    a3 = uc_lattice_vectors[3] .* N_a3
    push!(lattice_vectors, a1)
    push!(lattice_vectors, a2)
    push!(lattice_vectors, a3)

    # save everything to a Lattice object
    lattice = Lattice(
        unitcell,
        [N_a1, N_a2, N_a3],
        lattice_vectors,
        positions,
        positions_indices,
        connections,
        filename
        )
    # save the lattice object
    if save
        saveLattice(lattice)
    end

    # return the lattice object
    return lattice

end


"""
    getLatticePeriodic(unitcell::Unitcell, repetition_array::Array{Int64} [; save::Bool, load::Bool])

Function to construct a finite lattice with periodic boundary conditions on all sides out of a `Unitcell` object.
The number of unitcells that are put together in each elementery direction is passed in the
`repetition_array`.

Additionally, the newly created `Lattice` object can directly be saved. If this has been done before, passing a `load=true`
will allow to load the object instead of creating it again.

Note that this function works for both 2D and 3D unitcells.



# Examples

```julia-repl
julia> getLatticePeriodic(unitcell, [10, 20])
LatticePhysics.Lattice(...)

julia> getLatticePeriodic(unitcell, [10, 20], load=true)
LatticePhysics.Lattice(...)
```
"""
function getLatticePeriodic(unitcell::Unitcell, repetition_array::Array{Int64}; save::Bool=false, load::Bool=false)

    # check how many lattice vectors the unitcell has
    number_lv = size(unitcell.lattice_vectors,1)

    # determine which method to call
    if number_lv == 1
        println("Periodic lattices in 1D not implemented yet")
        return
    elseif number_lv == 2
        return getLatticePeriodic2D(unitcell, repetition_array, save=save, load=load)
    elseif number_lv == 3
        return getLatticePeriodic3D(unitcell, repetition_array, save=save, load=load)
    else
        println("Periodic lattices in dimensions larger 3D not implemented yet")
        return
    end

end

export getLatticePeriodic



#-----------------------------------------------------------------------------------------------------------------------------
#
#   Construction of OPEN Lattices (i.e. lattices, where all directions have open boundary conditions)
#
#-----------------------------------------------------------------------------------------------------------------------------

# FOR 2D AND 3D
function getLatticeOpen2D(unitcell::Unitcell, repetition_array::Array{Int64}; save=true, load=false)

    # extract the cardinal directions of the lattice from the array
    N_a1 = repetition_array[1]
    N_a2 = repetition_array[2]

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_open_$(N_a1)_$(N_a2).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "FOLDER_LATTICES$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
        filename = replace(filename, ".jld", "_open_$(N_a1)_$(N_a2).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    end

    # check if can be loaded
    if isfile(filename) && load
        return Lattice(filename)
    end

    # build the lattice

    # load the data from the unit cell
    uc_basis            = unitcell.basis
    uc_connections      = unitcell.connections
    uc_lattice_vectors  = unitcell.lattice_vectors


    # GENERATE NEW POSITIONS
	positions = Array[]
    positions_indices = []

	# define the index function to get the index of an element in the array
	function index(i,j,alpha)
		return size(uc_basis,1)*((i-1)*N_a2+j-1) + alpha
	end

    # define positions array to be filled
	for element in 1:N_a1*N_a2*size(uc_basis, 1)
		push!(positions, uc_basis[1])
		push!(positions_indices, 1)
	end
	# set all positions to their correct values
	for i in 1:N_a1
	for j in 1:N_a2
    for a in 1:size(uc_basis, 1)
		positions[index(i,j,a)] = uc_basis[a] + i*uc_lattice_vectors[1] + j*uc_lattice_vectors[2]
        positions_indices[index(i,j,a)] = a
	end
	end
    end

    # GENERATE NEW CONNECTIONS
	connections = Array[]

	# iterate over all unit cells
	for i in 1:N_a1
	for j in 1:N_a2
		# add all connections for unitcell (i,j)
		for connection in uc_connections
			# calculate the index from where the connection goes
			index_from = index(i,j,connection[1])
			# calculate the aimed unitcell
			i_to = i + connection[4][1]
			j_to = j + connection[4][2]
			# check if the connection goes around in a1 direction
            offset_a1 = 0
			while i_to < 1
                i_to += N_a1
                offset_a1 += -1
            end
            while i_to > N_a1
                i_to -= N_a1
                offset_a1 += 1
			end
			# check if the connection goes around in a2 direction
            offset_a2 = 0
			while j_to < 1
                j_to += N_a2
                offset_a2 += -1
            end
            while j_to > N_a2
                j_to -= N_a2
                offset_a2 += 1
			end
            # check if it is periodic, if yes: ignore
            if (offset_a1, offset_a2) != (0,0)
                continue
            end
			# get the index to where the connection goes
			index_to = index(i_to, j_to, connection[2])
            # generate a new connection
            # format [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            connection_new = [index_from; index_to; connection[3]; (offset_a1, offset_a2)]
			# register as connection
			push!(connections, connection_new)
		end
	end
	end


    # generate new lattice vectors, now spanning the WHOLE lattice (zero lattice vectors indicating no periodicity = finite)
    lattice_vectors = Array[]

    # save everything to a Lattice object
    lattice = Lattice(
        unitcell,
        [N_a1, N_a2],
        lattice_vectors,
        positions,
        positions_indices,
        connections,
        filename
        )
    # save the lattice object
    if save
        saveLattice(lattice)
    end

    # return the lattice
    return lattice

end
function getLatticeOpen3D(unitcell::Unitcell, repetition_array::Array{Int64}; save=true, load=false)

    # extract the cardinal directions of the lattice from the array
    N_a1 = repetition_array[1]
    N_a2 = repetition_array[2]
    N_a3 = repetition_array[3]

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_open_$(N_a1)_$(N_a2)_$(N_a3).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "FOLDER_LATTICES$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
        filename = replace(filename, ".jld", "_open_$(N_a1)_$(N_a2)_$(N_a3).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    end

    # check if can be loaded
    if isfile(filename) && load
        return Lattice(filename)
    end

    # build the lattice

    # load the data from the unit cell
    uc_basis            = unitcell.basis
    uc_connections      = unitcell.connections
    uc_lattice_vectors  = unitcell.lattice_vectors


    # GENERATE NEW POSITIONS
	positions = Array[]
    positions_indices = []

	# define the index function to get the index of an element in the array
	function index(i,j,k,alpha)
		return size(uc_basis,1)*((i-1)*N_a2*N_a3 + (j-1)*N_a3 + (k-1)) + alpha
	end

    # define positions array to be filled
	for element in 1:N_a1*N_a2*N_a3*size(uc_basis, 1)
		push!(positions, uc_basis[1])
		push!(positions_indices, 1)
	end
	# set all positions to their correct values
	for i in 1:N_a1
	for j in 1:N_a2
	for k in 1:N_a3
    for a in 1:size(uc_basis, 1)
		positions[index(i,j,k,a)] = uc_basis[a] + i*uc_lattice_vectors[1] + j*uc_lattice_vectors[2] + k*uc_lattice_vectors[3]
        positions_indices[index(i,j,k,a)] = a
	end
	end
    end
    end

    # GENERATE NEW CONNECTIONS
	connections = Array[]

	# iterate over all unit cells
	for i in 1:N_a1
	for j in 1:N_a2
	for k in 1:N_a3
		# add all connections for unitcell (i,j)
		for connection in uc_connections
			# calculate the index from where the connection goes
			index_from = index(i,j,k,connection[1])
			# calculate the aimed unitcell
			i_to = i + connection[4][1]
			j_to = j + connection[4][2]
			k_to = k + connection[4][3]
			# check if the connection goes around in a1 direction
            offset_a1 = 0
			while i_to < 1
                i_to += N_a1
                offset_a1 += -1
            end
            while i_to > N_a1
                i_to -= N_a1
                offset_a1 += 1
			end
			# check if the connection goes around in a2 direction
            offset_a2 = 0
			while j_to < 1
                j_to += N_a2
                offset_a2 += -1
            end
            while j_to > N_a2
                j_to -= N_a2
                offset_a2 += 1
			end
			# check if the connection goes around in a3 direction
            offset_a3 = 0
			while k_to < 1
                k_to += N_a3
                offset_a3 += -1
            end
            while k_to > N_a3
                k_to -= N_a3
                offset_a3 += 1
			end
            # check if it is periodic, if yes: ignore
            if (offset_a1, offset_a2, offset_a3) != (0,0,0)
                continue
            end
			# get the index to where the connection goes
			index_to = index(i_to, j_to, k_to, connection[2])
            # generate a new connection
            # format [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            connection_new = [index_from; index_to; connection[3]; (offset_a1, offset_a2, offset_a3)]
			# register as connection
			push!(connections, connection_new)
		end
	end
	end
    end


    # generate new lattice vectors, now spanning the WHOLE lattice (zero lattice vectors indicating no periodicity = finite)
    lattice_vectors = Array[]

    # save everything to a Lattice object
    lattice = Lattice(
        unitcell,
        [N_a1, N_a2, N_a3],
        lattice_vectors,
        positions,
        positions_indices,
        connections,
        filename
        )
    # save the lattice object
    if save
        saveLattice(lattice)
    end

    # return the lattice object
    return lattice

end

export getLatticeOpen2D
export getLatticeOpen3D


# FOR UNKNOWN DIMENSION
function getLatticeOpen(unitcell::Unitcell, repetition_array::Array{Int64}; save=true, load=false)

    # check how many lattice vectors the unitcell has
    number_lv = size(unitcell.lattice_vectors,1)

    # determine which method to call
    if number_lv == 1
        println("Open lattices in 1D not implemented yet")
        return
    elseif number_lv == 2
        return getLatticeOpen2D(unitcell, repetition_array, save=save, load=load)
    elseif number_lv == 3
        return getLatticeOpen3D(unitcell, repetition_array, save=save, load=load)
    else
        println("Open lattices in dimensions larger 3D not implemented yet")
        return
    end

end

export getLatticeOpen



#-----------------------------------------------------------------------------------------------------------------------------
#
#   Construction of SEMIPERIODIC Lattices
#   (i.e. lattices, where some directions have open and some have periodic boundary conditions)
#
#   negative numbers denote periodic directions
#   positive numbers denote open directions
#
#-----------------------------------------------------------------------------------------------------------------------------

# FOR 2D AND 3D
function getLatticeSemiperiodic2D(unitcell::Unitcell, repetition_array::Array{Int64}; save=true, load=false)

    # extract the cardinal directions of the lattice from the array
    N_a1 = repetition_array[1]
    N_a2 = repetition_array[2]

    # check if correct method is called
    if N_a1*N_a2 > 0
        # something is wrong, either flake or totally periodic
        println("wrong method called, either no periodic direction or all directions periodic desired, but can only deliver semiperiodic!")
        return
    end

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_semiperiodic_$(N_a1)_$(N_a2).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "FOLDER_LATTICES$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
        filename = replace(filename, ".jld", "_semiperiodic_$(N_a1)_$(N_a2).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    end

    # check if can be loaded
    if isfile(filename) && load
        return Lattice(filename)
    end

    # build the lattice

    # load the data from the unit cell
    uc_basis            = unitcell.basis
    uc_connections      = unitcell.connections
    uc_lattice_vectors  = unitcell.lattice_vectors

    # GENERATE NEW POSITIONS
	positions = Array[]
    positions_indices = []

    # check which direction periodic
    if N_a1 < 0
        N_a1 *= -1
        periodic_direction = 1
    else
        N_a2 *= -1
        periodic_direction = 2
    end

	# define the index function to get the index of an element in the array
	function index(i,j,alpha)
		return size(uc_basis,1)*((i-1)*N_a2+j-1) + alpha
	end

    # define positions array to be filled
	for element in 1:N_a1*N_a2*size(uc_basis, 1)
		push!(positions, uc_basis[1])
        push!(positions_indices, 1)
	end
	# set all positions to their correct values
	for i in 1:N_a1
	for j in 1:N_a2
    for a in 1:size(uc_basis, 1)
		positions[index(i,j,a)] = uc_basis[a] + i*uc_lattice_vectors[1] + j*uc_lattice_vectors[2]
        positions_indices[index(i,j,a)] = a
	end
	end
    end

    # GENERATE NEW CONNECTIONS
	connections = Array[]

	# iterate over all unit cells
	for i in 1:N_a1
	for j in 1:N_a2
		# add all connections for unitcell (i,j)
		for connection in uc_connections
			# calculate the index from where the connection goes
			index_from = index(i,j,connection[1])
			# calculate the aimed unitcell
			i_to = i + connection[4][1]
			j_to = j + connection[4][2]
			# check if the connection goes around in a1 direction
            offset_a1 = 0
			while i_to < 1
                i_to += N_a1
                offset_a1 += -1
            end
            while i_to > N_a1
                i_to -= N_a1
                offset_a1 += 1
			end
			# check if the connection goes around in a2 direction
            offset_a2 = 0
			while j_to < 1
                j_to += N_a2
                offset_a2 += -1
            end
            while j_to > N_a2
                j_to -= N_a2
                offset_a2 += 1
			end
            # check if it outside the open boundary
            if periodic_direction == 1 && offset_a2 != 0
                continue
            elseif periodic_direction == 2 && offset_a1 != 0
                continue
            end
			# get the index to where the connection goes
			index_to = index(i_to, j_to, connection[2])
            # generate a new connection
            # format [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            if periodic_direction == 1
                connection_new = [index_from; index_to; connection[3]; (offset_a1)]
            else
                connection_new = [index_from; index_to; connection[3]; (offset_a2)]
            end
			# register as connection
			push!(connections, connection_new)
		end
	end
	end

    # generate new lattice vectors, now spanning the WHOLE lattice (zero lattice vectors indicating no periodicity = finite)
    lattice_vectors = Array[]
    if periodic_direction == 1
        a1 = uc_lattice_vectors[1] .* N_a1
        push!(lattice_vectors, a1)
    else
        a2 = uc_lattice_vectors[2] .* N_a2
        push!(lattice_vectors, a2)
    end

    # save everything to a Lattice object
    lattice = Lattice(
        unitcell,
        [N_a1, N_a2],
        lattice_vectors,
        positions,
        positions_indices,
        connections,
        filename
        )
    # save the lattice object
    if save
        saveLattice(lattice)
    end

    # return the lattice
    return lattice
end
function getLatticeSemiperiodic3D(unitcell::Unitcell, repetition_array::Array{Int64}; save=true, load=false)

    # extract the cardinal directions of the lattice from the array
    N_a1 = repetition_array[1]
    N_a2 = repetition_array[2]
    N_a3 = repetition_array[3]

    # check which directions are periodic
    periodic_directions = []
    if N_a1 < 0
        push!(periodic_directions, 1)
    end
    if N_a2 < 0
        push!(periodic_directions, 2)
    end
    if N_a3 < 0
        push!(periodic_directions, 3)
    end

    # check if correct method called
    if length(periodic_directions) in [0,3]
        # something is wrong, either flake or totally periodic
        println("wrong method called, either no periodic direction or all directions periodic desired, but can only deliver semiperiodic!")
        return
    end

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_semiperiodic_$(N_a1)_$(N_a2)_$(N_a3).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "FOLDER_LATTICES$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
        filename = replace(filename, ".jld", "_semiperiodic_$(N_a1)_$(N_a2)_$(N_a3).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    end

    # check if can be loaded
    if isfile(filename) && load
        return Lattice(filename)
    end

    # build the lattice

    # load the data from the unit cell
    uc_basis            = unitcell.basis
    uc_connections      = unitcell.connections
    uc_lattice_vectors  = unitcell.lattice_vectors


    # GENERATE NEW POSITIONS
	positions = Array[]
	positions_indices = []

    # turn all directions to positive numbers
    N_a1 = abs(N_a1)
    N_a2 = abs(N_a2)
    N_a3 = abs(N_a3)

	# define the index function to get the index of an element in the array
	function index(i,j,k,alpha)
		return size(uc_basis,1)*((i-1)*N_a2*N_a3 + (j-1)*N_a3 + (k-1)) + alpha
	end

    # define positions array to be filled
	for element in 1:N_a1*N_a2*N_a3*size(uc_basis, 1)
		push!(positions, uc_basis[1])
        push!(positions_indices, 1)
	end
	# set all positions to their correct values
	for i in 1:N_a1
	for j in 1:N_a2
	for k in 1:N_a3
    for a in 1:size(uc_basis, 1)
		positions[index(i,j,k,a)] = uc_basis[a] + i*uc_lattice_vectors[1] + j*uc_lattice_vectors[2] + k*uc_lattice_vectors[3]
        positions_indices[index(i,j,k,a)] = a
	end
	end
    end
    end

    # GENERATE NEW CONNECTIONS
	connections = Array[]

	# iterate over all unit cells
	for i in 1:N_a1
	for j in 1:N_a2
	for k in 1:N_a3
		# add all connections for unitcell (i,j)
		for connection in uc_connections
			# calculate the index from where the connection goes
			index_from = index(i,j,k,connection[1])
			# calculate the aimed unitcell
			i_to = i + connection[4][1]
			j_to = j + connection[4][2]
			k_to = k + connection[4][3]
			# check if the connection goes around in a1 direction
            offset_a1 = 0
			while i_to < 1
                i_to += N_a1
                offset_a1 += -1
            end
            while i_to > N_a1
                i_to -= N_a1
                offset_a1 += 1
			end
			# check if the connection goes around in a2 direction
            offset_a2 = 0
			while j_to < 1
                j_to += N_a2
                offset_a2 += -1
            end
            while j_to > N_a2
                j_to -= N_a2
                offset_a2 += 1
			end
			# check if the connection goes around in a3 direction
            offset_a3 = 0
			while k_to < 1
                k_to += N_a3
                offset_a3 += -1
            end
            while k_to > N_a3
                k_to -= N_a3
                offset_a3 += 1
			end
            # offsets
            offsets = []
            # check if the offsets are correct
            if offset_a1 != 0 && !(1 in periodic_directions)
                continue
            elseif 1 in periodic_directions
                push!(offsets, offset_a1)
            end
            if offset_a2 != 0 && !(2 in periodic_directions)
                continue
            elseif 2 in periodic_directions
                push!(offsets, offset_a2)
            end
            if offset_a3 != 0 && !(3 in periodic_directions)
                continue
            elseif 3 in periodic_directions
                push!(offsets, offset_a3)
            end
			# get the index to where the connection goes
			index_to = index(i_to, j_to, k_to, connection[2])
            # generate a new connection
            # format [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            if length(offsets) == 1
                connection_new = [index_from; index_to; connection[3]; offsets[1]]
            else # offset has to have 2 entries
                connection_new = [index_from; index_to; connection[3]; (offsets[1], offsets[2])]
            end
			# register as connection
			push!(connections, connection_new)
		end
	end
	end
    end


    # generate new lattice vectors, now spanning the WHOLE lattice
    lattice_vectors = Array[]
    a1 = uc_lattice_vectors[1] .* N_a1
    a2 = uc_lattice_vectors[2] .* N_a2
    a3 = uc_lattice_vectors[3] .* N_a3
    if 1 in periodic_directions
        push!(lattice_vectors, a1)
    end
    if 2 in periodic_directions
        push!(lattice_vectors, a2)
    end
    if 3 in periodic_directions
        push!(lattice_vectors, a3)
    end

    # save everything to a Lattice object
    lattice = Lattice(
        unitcell,
        [N_a1, N_a2, N_a3],
        lattice_vectors,
        positions,
        positions_indices,
        connections,
        filename
        )
    # save the lattice object
    if save
        saveLattice(lattice)
    end

    # return the lattice object
    return lattice
end

export getLatticeSemiperiodic2D
export getLatticeSemiperiodic3D


# FOR UNKNOWN DIMENSION
function getLatticeSemiperiodic(unitcell::Unitcell, repetition_array::Array{Int64}; save=true, load=false)

    # check how many lattice vectors the unitcell has
    number_lv = size(unitcell.lattice_vectors,1)

    # determine which method to call
    if number_lv == 1
        println("Semiperiodic lattices in 1D not implemented yet")
        return
    elseif number_lv == 2
        return getLatticeSemiperiodic2D(unitcell, repetition_array, save=save, load=load)
    elseif number_lv == 3
        return getLatticeSemiperiodic3D(unitcell, repetition_array, save=save, load=load)
    else
        println("Semiperiodic lattices in dimensions larger 3D not implemented yet")
        return
    end

end

export getLatticeSemiperiodic




#-----------------------------------------------------------------------------------------------------------------------------
#
#   Construction of ALL Lattices (that are build from stacking/repeating unitcells)
#
#   negative numbers denote periodic directions
#   positive numbers denote open directions
#
#-----------------------------------------------------------------------------------------------------------------------------
function getLattice(unitcell::Unitcell, repetition_array::Array{Int64}; save=true, load=false)

    # check how many periodic directions
    number_pd = 0
    for N in repetition_array
        if N < 0
            number_pd = number_pd + 1
        end
    end

    # check how many open directions
    number_od = 0
    for N in repetition_array
        if N > 0
            number_od = number_od + 1
        end
    end

    # check if all directions are either periodic or open
    if number_od + number_pd != length(repetition_array)
        println("Some directions are not specified correctly! Abort...")
        return
    end

    # check which method to call
    if number_pd == 0
        # no periodc directions present, OPEN
        return getLatticeOpen(unitcell, repetition_array, save=save, load=load)
    elseif number_od == 0
        # no open directions present, PERIODIC (but multiply all ranges by -1)
        return getLatticePeriodic(unitcell, repetition_array.*-1, save=save, load=load)
    else
        # both periodic and open directions present, SEMIPERIODIC
        return getLatticeSemiperiodic(unitcell, repetition_array, save=save, load=load)
    end

end
function getLattice(unitcell::Unitcell, repetitions::Int64; save=true, load=false)
    # just parse through to more general method
    return getLattice(unitcell, repetitions.*ones(size(unitcell.lattice_vectors,1)), save=save, load=load)
end


export getLattice










#-----------------------------------------------------------------------------------------------------------------------------
#
#   BUILDING FUNCTIONS FOR LATTICES THAT ARE BUILD BY BOND DISTANCE FROM AN ORIGIN SITE
#   These lattices have of course open boundaries
#
#   Parameters that have to be passed:
#   - Unitcell object from which lattice is build
#   - Integer with extent of the lattice from the origin site outwards
#   - OPTIONAL: origin site index
#
#-----------------------------------------------------------------------------------------------------------------------------


# generate a lattice in 2D by bond distance (open boundary conditions)
function getLatticeByBondDistance2D(unitcell::Unitcell, bonddistance::Int64; origin::Int64=1, load=false, save=true)

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_by_bonddistance_$(bonddistance)_from_$(origin).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "FOLDER_LATTICES$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
        filename = replace(filename, ".jld", "_by_bonddistance_$(bonddistance)_from_$(origin).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    end

    # check if the filename already exists
    if isfile(filename) && load
        # return the lattice
        return Lattice(filename)
    end

    # load the data from the unit cell
    uc_basis            = unitcell.basis
    uc_connections      = unitcell.connections
    uc_lattice_vectors  = unitcell.lattice_vectors

    # arrays for new positions and connections
    positions   = []
    connections = Array[]

    # checklist for all sites that are checked if added etc
    checklist   = []
    # format for sites on checklist: ([pos_x, pos_y], [uc_i, uc_j], index_in_UC, bd_current)

    # push the origin site to the checklist
    push!(
        checklist, (uc_basis[origin], [0,0], origin, 0)
    )

    # iterate while the checklist is not empty
    while size(checklist,1) > 0

        # get the item that was on the checklist for longest
        item_to_handle = shift!(checklist)

        # check if the item is already in the positions list
        found = false
        for p in positions
            if item_to_handle[2] == p[2] && item_to_handle[3] == p[3]
                # if yes, continue
                found = true
                break
            end
        end
        if found
            continue
        end

        # if not, push it into
        push!(positions, item_to_handle)
        # for all connections
        index_from = size(positions, 1)

        # insert all connections to sites that are already inside the positions list and all other sites into the checklist
        for c in uc_connections
            # check if correct connections
            if c[1] != item_to_handle[3]
                continue
            end
            # search for the other element of the connection
            i_to = item_to_handle[2][1] + c[4][1]
            j_to = item_to_handle[2][2] + c[4][2]
            a_to = c[2]
            index_to = -1
            for (index,item_handled) in enumerate(positions)
                if item_handled[2] == [i_to, j_to] && item_handled[3] == a_to
                    # make a new connection
                    index_to = index
                    # break the loop
                    break
                end
            end
            # determine whether the element is already inside the list
            if index_to > 0
                # element is at index index_to an can be linked
                connection_new = [index_from; index_to; c[3]; (0, 0)]
                # check if the connetion is already added
                if connection_new in connections
                    continue
                end
                # register as connection
                push!(connections, connection_new)
            else
                # element not in list yet, maybe should be added
                if item_to_handle[4] < bonddistance
                    # format for sites on checklist: ([pos_x, pos_y], [uc_i, uc_j], index_in_UC, bd_current)
                    push!(checklist, (
                        (i_to * uc_lattice_vectors[1]) .+ (j_to * uc_lattice_vectors[2]) .+ uc_basis[c[2]],
                        [i_to, j_to],
                        c[2],
                        item_to_handle[4]+1
                    ))
                end
            end
        end


    end

    # change the format of positions
    positions_TMP = positions

    # erase positions
    positions = Array[]
    positions_indices = []

    # insert the real positions
    for p in positions_TMP
        push!(positions, p[1])
        push!(positions_indices, p[3])
    end

    # insert missing connections (if (i to j) is present, insert (j to i))
    for c in connections
        c_proposed = [c[2]; c[1]; c[3]; (0, 0)]
        if !(c_proposed in connections)
            push!(connections, c_proposed)
        else
            #println("connection already there")
        end
    end

    # save everything to a Lattice object
    lattice = Lattice(
        unitcell,
        [],
        [],
        positions,
        positions_indices,
        connections,
        filename
        )
    # save the lattice object
    if save
        saveLattice(lattice)
    end

    # return the lattice object
    return lattice

end
# generate a lattice in 3D by bond distance (open boundary conditions)
function getLatticeByBondDistance3D(unitcell::Unitcell, bonddistance::Int64; origin::Int64=1, load=false, save=true)

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_by_bonddistance_$(bonddistance)_from_$(origin).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "FOLDER_LATTICES$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
        filename = replace(filename, ".jld", "_by_bonddistance_$(bonddistance)_from_$(origin).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    end

    # check if the filename already exists
    if isfile(filename) && load
        # return the lattice
        return Lattice(filename)
    end

    # load the data from the unit cell
    uc_basis            = unitcell.basis
    uc_connections      = unitcell.connections
    uc_lattice_vectors  = unitcell.lattice_vectors

    # arrays for new positions and connections
    positions   = []
    connections = Array[]

    # checklist for all sites that are checked if added etc
    checklist   = []
    # format for sites on checklist: ([pos_x, pos_y, pos_z], [uc_i, uc_j, uc_k], index_in_UC, bd_current)

    # push the origin site to the checklist
    push!(
        checklist, (uc_basis[origin], [0,0,0], origin, 0)
    )

    # iterate while the checklist is not empty
    while size(checklist,1) > 0

        # get the item that was on the checklist for longest
        item_to_handle = shift!(checklist)

        # check if the item is already in the positions list
        found = false
        for p in positions
            if item_to_handle[2] == p[2] && item_to_handle[3] == p[3]
                # if yes, continue
                found = true
                break
            end
        end
        if found
            continue
        end

        # if not, push it into
        push!(positions, item_to_handle)
        # for all connections
        index_from = size(positions, 1)

        # insert all connections to sites that are already inside the positions list and all other sites into the checklist
        for c in uc_connections
            # check if correct connections
            if c[1] != item_to_handle[3]
                continue
            end
            # search for the other element of the connection
            i_to = item_to_handle[2][1] + c[4][1]
            j_to = item_to_handle[2][2] + c[4][2]
            k_to = item_to_handle[2][3] + c[4][3]
            a_to = c[2]
            index_to = -1
            for (index,item_handled) in enumerate(positions)
                if item_handled[2] == [i_to, j_to, k_to] && item_handled[3] == a_to
                    # make a new connection
                    index_to = index
                    # break the loop
                    break
                end
            end
            # determine whether the element is already inside the list
            if index_to > 0
                # element is at index index_to an can be linked
                connection_new = [index_from; index_to; c[3]; (0, 0, 0)]
                # check if the connetion is already added
                if connection_new in connections
                    continue
                end
                # register as connection
                push!(connections, connection_new)
            else
                # element not in list yet, maybe should be added
                if item_to_handle[4] < bonddistance
                    # format for sites on checklist: ([pos_x, pos_y, pos_z], [uc_i, uc_j, uc_k], index_in_UC, bd_current)
                    push!(checklist, (
                        (i_to * uc_lattice_vectors[1]) .+ (j_to * uc_lattice_vectors[2]) .+ (k_to * uc_lattice_vectors[3]) .+ uc_basis[c[2]],
                        [i_to, j_to, k_to],
                        c[2],
                        item_to_handle[4]+1
                    ))
                end
            end
        end
    end

    # change the format of positions
    positions_TMP = positions

    # erase positions
    positions = Array[]
    positions_indices = []

    # insert the real positions
    for p in positions_TMP
        push!(positions, p[1])
        push!(positions_indices, p[3])
    end

    # insert missing connections (if (i to j) is present, insert (j to i))
    for c in connections
        c_proposed = [c[2]; c[1]; c[3]; (0, 0, 0)]
        if !(c_proposed in connections)
            push!(connections, c_proposed)
        else
            #println("connection already there")
        end
    end

    # save everything to a Lattice object
    lattice = Lattice(
        unitcell,
        [],
        [],
        positions,
        positions_indices,
        connections,
        filename
        )
    # save the lattice object
    if save
        saveLattice(lattice)
    end

    # return the lattice object
    return lattice
end

# generate a flake of any dimension by bond distance
function getLatticeByBondDistance(unitcell::Unitcell, bonddistance::Int64; origin::Int64=1, load=false, save=true)

    # check how many periodic dimensions the unitcell has
    N_dims = size(unitcell.lattice_vectors, 1)

    # depends on dimensions what to call
    if N_dims == 0
        # not possible
        println("cannot build a lattice without periodic lattice vectors!")
        return
    elseif N_dims == 1
        # not implemented yet
        println("building from bond distance not implemented yet for d=1")
    elseif N_dims == 2
        # just call the fitting routine with the first two entries
        return getLatticeByBondDistance2D(unitcell, bonddistance, origin=origin, load=load, save=save)
    elseif N_dims == 3
        # just call the fitting routine with the first two entries
        return getLatticeByBondDistance3D(unitcell, bonddistance, origin=origin, load=load, save=save)
    else
        # not implemented for any other lattice!
        println("building from bond distance not implemented for d=$(N_dims)")
        return
    end

end


export getLatticeByBondDistance2D
export getLatticeByBondDistance3D

export getLatticeByBondDistance






#-----------------------------------------------------------------------------------------------------------------------------
#
#   BUILDING FUNCTIONS FOR LATTICES THAT ARE BUILD BY BOND DISTANCE FROM AN ORIGIN SITE AND RESIDE IN A GIVEN SHAPE
#   These lattices have of course open boundaries
#
#   Parameters that have to be passed:
#   - Unitcell object from which lattice is build
#   - shape to determine the extent of the lattice from the origin site outwards
#     (a julia function that gives true / false if a site is inside the shape)
#   - the name of the shape for filename purposes
#   - OPTIONAL: origin site index
#
#-----------------------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------------------
#
#   General building in shape for 2D and 3D and general lattices
#
#-----------------------------------------------------------------------------------------------------------------------------
function getLatticeInShape2D(unitcell::Unitcell, shape::Function, shapename::String; origin::Int64=1, load=false, save=true)

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_in_shape_$(shapename).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "FOLDER_LATTICES$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
        filename = replace(filename, ".jld", "_in_shape_$(shapename).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    end

    # check if can be loaded
    if isfile(filename) && load
        return Lattice(filename)
    end

    # load the data from the unit cell
    uc_basis            = unitcell.basis
    uc_connections      = unitcell.connections
    uc_lattice_vectors  = unitcell.lattice_vectors

    # arrays for new positions and connections
    positions   = []
    connections = Array[]

    # checklist for all sites that are checked if added etc
    checklist   = []
    # format for sites on checklist: ([pos_x, pos_y], [uc_i, uc_j], index_in_UC, bd_current)

    # push the origin site to the checklist
    push!(
        checklist, (uc_basis[origin], [0,0], origin, 0)
    )

    # iterate while the checklist is not empty
    while size(checklist,1) > 0

        # get the item that was on the checklist for longest
        item_to_handle = shift!(checklist)

        # check if the item is already in the positions list
        found = false
        for p in positions
            if item_to_handle[2] == p[2] && item_to_handle[3] == p[3]
                # if yes, continue
                found = true
                break
            end
        end
        if found
            continue
        end

        # if not, push it into
        push!(positions, item_to_handle)
        # for all connections
        index_from = size(positions, 1)

        # insert all connections to sites that are already inside the positions list and all other sites into the checklist
        for c in uc_connections
            # check if correct connections
            if c[1] != item_to_handle[3]
                continue
            end
            # search for the other element of the connection
            i_to = item_to_handle[2][1] + c[4][1]
            j_to = item_to_handle[2][2] + c[4][2]
            a_to = c[2]
            index_to = -1
            for (index,item_handled) in enumerate(positions)
                if item_handled[2] == [i_to, j_to] && item_handled[3] == a_to
                    # make a new connection
                    index_to = index
                    # break the loop
                    break
                end
            end
            # determine whether the element is already inside the list
            if index_to > 0
                # element is at index index_to an can be linked
                connection_new = [index_from; index_to; c[3]; (0, 0)]
                # check if the connetion is already added
                if connection_new in connections
                    continue
                end
                # register as connection
                push!(connections, connection_new)
            else
                # element not in list yet, maybe should be added
                if shape((i_to * uc_lattice_vectors[1]) .+ (j_to * uc_lattice_vectors[2]) .+ uc_basis[c[2]] .- uc_basis[origin])
                    # format for sites on checklist: ([pos_x, pos_y], [uc_i, uc_j], index_in_UC, bd_current)
                    push!(checklist, (
                        (i_to * uc_lattice_vectors[1]) .+ (j_to * uc_lattice_vectors[2]) .+ uc_basis[c[2]],
                        [i_to, j_to],
                        c[2],
                        item_to_handle[4]+1
                    ))
                end
            end
        end


    end

    # change the format of positions
    positions_TMP = positions

    # erase positions
    positions = Array[]
    positions_indices = []

    # insert the real positions
    for p in positions_TMP
        push!(positions, p[1])
        push!(positions_indices, p[3])
    end

    # insert missing connections (if (i to j) is present, insert (j to i))
    for c in connections
        c_proposed = [c[2]; c[1]; c[3]; (0, 0)]
        if !(c_proposed in connections)
            push!(connections, c_proposed)
        else
            #println("connection already there")
        end
    end

    # save everything to a Lattice object
    lattice = Lattice(
        unitcell,
        [],
        [],
        positions,
        positions_indices,
        connections,
        filename
        )
    # save the lattice object
    if save
        saveLattice(lattice)
    end

    # return the lattice object
    return lattice

end
function getLatticeInShape3D(unitcell::Unitcell, shape::Function, shapename::String; origin::Int64=1, load=false, save=true)

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_in_shape_$(shapename).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "FOLDER_LATTICES$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
        filename = replace(filename, ".jld", "_in_shape_$(shapename).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    end

    # check if can be loaded
    if isfile(filename) && load
        return Lattice(filename)
    end

    # load the data from the unit cell
    uc_basis            = unitcell.basis
    uc_connections      = unitcell.connections
    uc_lattice_vectors  = unitcell.lattice_vectors

    # arrays for new positions and connections
    positions   = []
    connections = Array[]

    # checklist for all sites that are checked if added etc
    checklist   = []
    # format for sites on checklist: ([pos_x, pos_y, pos_z], [uc_i, uc_j, uc_k], index_in_UC, bd_current)

    # push the origin site to the checklist
    push!(
        checklist, (uc_basis[origin], [0,0,0], origin, 0)
    )

    # iterate while the checklist is not empty
    while size(checklist,1) > 0

        # get the item that was on the checklist for longest
        item_to_handle = shift!(checklist)

        # check if the item is already in the positions list
        found = false
        for p in positions
            if item_to_handle[2] == p[2] && item_to_handle[3] == p[3]
                # if yes, continue
                found = true
                break
            end
        end
        if found
            continue
        end

        # if not, push it into
        push!(positions, item_to_handle)
        # for all connections
        index_from = size(positions, 1)

        # insert all connections to sites that are already inside the positions list and all other sites into the checklist
        for c in uc_connections
            # check if correct connections
            if c[1] != item_to_handle[3]
                continue
            end
            # search for the other element of the connection
            i_to = item_to_handle[2][1] + c[4][1]
            j_to = item_to_handle[2][2] + c[4][2]
            k_to = item_to_handle[2][3] + c[4][3]
            a_to = c[2]
            index_to = -1
            for (index,item_handled) in enumerate(positions)
                if item_handled[2] == [i_to, j_to, k_to] && item_handled[3] == a_to
                    # make a new connection
                    index_to = index
                    # break the loop
                    break
                end
            end
            # determine whether the element is already inside the list
            if index_to > 0
                # element is at index index_to an can be linked
                connection_new = [index_from; index_to; c[3]; (0, 0, 0)]
                # check if the connetion is already added
                if connection_new in connections
                    continue
                end
                # register as connection
                push!(connections, connection_new)
            else
                # element not in list yet, maybe should be added
                if shape((i_to * uc_lattice_vectors[1]) .+ (j_to * uc_lattice_vectors[2]) .+ (k_to * uc_lattice_vectors[3]) .+ uc_basis[c[2]] .- uc_basis[origin])
                    # format for sites on checklist: ([pos_x, pos_y, pos_z], [uc_i, uc_j, uc_k], index_in_UC, bd_current)
                    push!(checklist, (
                        (i_to * uc_lattice_vectors[1]) .+ (j_to * uc_lattice_vectors[2]) .+ (k_to * uc_lattice_vectors[3]) .+ uc_basis[c[2]],
                        [i_to, j_to, k_to],
                        c[2],
                        item_to_handle[4]+1
                    ))
                end
            end
        end
    end

    # change the format of positions
    positions_TMP = positions

    # erase positions
    positions = Array[]
    positions_indices = []

    # insert the real positions
    for p in positions_TMP
        push!(positions, p[1])
        push!(positions_indices, p[3])
    end

    # insert missing connections (if (i to j) is present, insert (j to i))
    for c in connections
        c_proposed = [c[2]; c[1]; c[3]; (0, 0, 0)]
        if !(c_proposed in connections)
            push!(connections, c_proposed)
        else
            #println("connection already there")
        end
    end

    # save everything to a Lattice object
    lattice = Lattice(
        unitcell,
        [],
        [],
        positions,
        positions_indices,
        connections,
        filename
        )
    # save the lattice object
    if save
        saveLattice(lattice)
    end

    # return the lattice object
    return lattice

end

function getLatticeInShape(unitcell::Unitcell, shape::Function, shapename::String; origin::Int64=1, load=false, save=true)

    # check how many periodic dimensions the unitcell has
    N_dims = size(unitcell.lattice_vectors, 1)

    # depends on dimensions what to call
    if N_dims == 0
        # not possible
        println("cannot build a lattice without periodic lattice vectors!")
        return
    elseif N_dims == 1
        # not implemented yet
        println("building in shape not implemented yet for d=1")
    elseif N_dims == 2
        # just call the fitting routine with the first two entries
        return getLatticeInShape2D(unitcell, shape, shapename, origin=origin, load=load, save=save)
    elseif N_dims == 3
        # just call the fitting routine with the first two entries
        return getLatticeInShape3D(unitcell, shape, shapename, origin=origin, load=load, save=save)
    else
        # not implemented for any other lattice!
        println("building in shape not implemented for d=$(N_dims)")
        return
    end
end


export getLatticeInShape2D
export getLatticeInShape3D

export getLatticeInShape



#-----------------------------------------------------------------------------------------------------------------------------
#
#   Special cases for building in shape for 2D and 3D and general lattices
#
#-----------------------------------------------------------------------------------------------------------------------------

# SPECIAL CASE: SPHERE / CIRCLE with radius
function getLatticeInSphere(unitcell::Unitcell, radius::Float64; origin::Int64=1, load=false, save=true)

    # check how many periodic dimensions the unitcell has
    N_dims = size(unitcell.lattice_vectors, 1)

    # depends on dimensions what to call
    if N_dims == 0
        # not possible
        println("cannot build a lattice without periodic lattice vectors!")
        return
    elseif N_dims == 1
        # not implemented yet
        println("building in a sphere not implemented yet for d=1")
        return
    elseif N_dims == 2
        # determine the shape name
        shapename = "sphere_$(radius)"
        # determine the shape function
        shape_circle(point) = sum(point.*point) < radius
        # call the general shape function
        return getLatticeInShape2D(unitcell, shape_circle, shapename, origin=origin, load=load, save=save)
    elseif N_dims == 3
        # determine the shape name
        shapename = "sphere_$(radius)"
        # determine the shape function
        shape_sphere(point) = sum(point.*point) < radius
        # call the general shape function
        return getLatticeInShape3D(unitcell, shape_sphere, shapename, origin=origin, load=load, save=save)
    else
        # not implemented for any other lattice!
        println("building in a sphere not implemented for d=$(N_dims)")
        return
    end
end
export getLatticeInSphere


# SPECIAL CASE: BOX / RECTANGLE with extent_array denoting the length of the different sides of the box (centered around the origin)
function getLatticeInBox(unitcell::Unitcell, extent_array::Array{Float64}; origin::Int64=1, load=false, save=true)

    # check how many periodic dimensions the unitcell has
    N_dims = size(unitcell.lattice_vectors, 1)

    # depends on dimensions what to call
    if N_dims == 0
        # not possible
        println("cannot build a lattice without periodic lattice vectors!")
        return
    elseif N_dims == 1
        # not implemented yet
        println("building in box not implemented yet for d=1")
        return
    elseif N_dims == 2
        # get the dimensions
        length_x = extent_array[1]
        length_y = extent_array[2]
        # determine the shape name
        shapename = "box_$(length_x)x$(length_y)_around_$(origin)"
        # determine the shape function
        shape_box2d(point) = (abs(point[1])<=length_x/2.0) && (abs(point[2])<=length_y/2.0)
        # call the general shape function
        return getLatticeInShape2D(unitcell, shape_box2d, shapename, origin=origin, load=load, save=save)
    elseif N_dims == 3
        # get the dimensions
        length_x = extent_array[1]
        length_y = extent_array[2]
        length_z = extent_array[3]
        # determine the shape name
        shapename = "box_$(length_x)x$(length_y)x$(length_z)_around_$(origin)"
        # determine the shape function
        shape_box3d(point) = (abs(point[1])<=length_x/2.0) && (abs(point[2])<=length_y/2.0) && (abs(point[3])<=length_z/2.0)
        # call the general shape function
        return getLatticeInShape3D(unitcell, shape_box3d, shapename, origin=origin, load=load, save=save)
    else
        # not implemented for any other lattice!
        println("building in box not implemented for d=$(N_dims)")
        return
    end
end
export getLatticeInBox







#
# TODO
#
# getLatticeGrapheneFlake
# getLatticeTriangularFlake
#























#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#
#   LATTICE OPERATIONS (MODIFY LATTICE CONNECTIONS AND POSITIONS)
#
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------------------
#
#   Build a lattice-X / unitcell-X version of the lattice / unitcell
#   Procedure replaces all bonds by sites
#
#-----------------------------------------------------------------------------------------------------------------------------
function getTransformedLatticeBondToSite(lattice::Lattice; connection_strength="AUTO")
    # new positions and connections
    positions = copy(lattice.positions)
    connections = Array[]
    connectionsTreated = Array[]
    # generate a new filename
    filename_new = replace(lattice.filename, ".jld", "_X.jld")
    # neutral wrap
    neutral_wrap = (0,0,0)
    if size(lattice.lattice_vectors,1) == 2
   	    neutral_wrap = (0,0)
    elseif size(lattice.lattice_vectors,1) == 1
   	    neutral_wrap = (0)
    end
    # iterate over all connections
    for c in lattice.connections
        # only treat one of the two
        treated = false
        for ct in connectionsTreated
            if ct[1] == c[2] && ct[2] == c[1] && ct[3] == c[3]
                treated = true
                break
            end
        end
        if treated
            continue
        end
        # add a new position
        pnew = positions[Int(c[1])].+positions[Int(c[2])]
        for i in 1:length(c[4])
            pnew .+= c[4][i].*lattice.lattice_vectors[i]
        end
        pnew = pnew .* 0.5
        push!(positions, pnew)
        # add new connections
        if typeof(connection_strength) == String && connection_strength == "AUTO"
            con_strength = "con$(size(connectionsTreated,1))"
        elseif typeof(connection_strength) == String && connection_strength == "SAME"
            con_strength = c[3]
        elseif typeof(connection_strength) == String && connection_strength == "SQRT"
            if typeof(c[3]) == String
                con_strength = "sqrt($(c[3]))"
            else
                con_strength = sqrt(c[3])
            end
        else
            con_strength = connection_strength
        end
        push!(connections, [c[1], size(positions,1), con_strength, neutral_wrap])
        push!(connections, [size(positions,1), c[1], con_strength, neutral_wrap])
        push!(connections, [size(positions,1), c[2], con_strength, c[4]])
        if length(c[4]) == 3
            push!(connections, [c[2], size(positions,1), con_strength, (-c[4][1],-c[4][2],-c[4][3])])
        elseif length(c[4]) == 2
            push!(connections, [c[2], size(positions,1), con_strength, (-c[4][1],-c[4][2])])
        else
            push!(connections, [c[2], size(positions,1), con_strength, (-c[4])])
        end
        # add to treated
        push!(connectionsTreated, c)
    end
    # build a new lattice
    lattice_X = Lattice(lattice.lattice_vectors, positions, connections, filename_new)
    # return the new lattice
    return lattice_X
end
export getTransformedLatticeBondToSite

function getTransformedUnitcellBondToSite(unitcell::Unitcell; connection_strength="AUTO")
    # new positions and connections
    positions = copy(unitcell.basis)
    connections = Array[]
    connectionsTreated = Array[]
    # generate a new filename
    filename_new = replace(unitcell.filename, ".jld", "_X.jld")
    # neutral wrap
    neutral_wrap = (0,0,0)
    if size(unitcell.lattice_vectors,1) == 2
   	    neutral_wrap = (0,0)
    elseif size(unitcell.lattice_vectors,1) == 1
   	    neutral_wrap = (0)
    end
    # iterate over all connections
    for c in unitcell.connections
        # only treat one of the two
        treated = false
        for ct in connectionsTreated
            if ct[1] == c[2] && ct[2] == c[1] && ct[3] == c[3]
                treated = true
                break
            end
        end
        if treated
            continue
        end
        # add a new position
        pnew = positions[Int(c[1])].+positions[Int(c[2])]
        for i in 1:length(c[4])
            pnew .+= c[4][i].*unitcell.lattice_vectors[i]
        end
        pnew = pnew .* 0.5
        push!(positions, pnew)
        # add new connections
        if typeof(connection_strength) == String && connection_strength == "AUTO"
            con_strength = "con$(size(connectionsTreated,1))"
        elseif typeof(connection_strength) == String && connection_strength == "SAME"
            con_strength = c[3]
        elseif typeof(connection_strength) == String && connection_strength == "SQRT"
            if typeof(c[3]) == String
                con_strength = "sqrt($(c[3]))"
            else
                con_strength = sqrt(c[3])
            end
        else
            con_strength = connection_strength
        end
        push!(connections, [c[1], size(positions,1), con_strength, neutral_wrap])
        push!(connections, [size(positions,1), c[1], con_strength, neutral_wrap])
        push!(connections, [size(positions,1), c[2], con_strength, c[4]])
        if length(c[4]) == 3
            push!(connections, [c[2], size(positions,1), con_strength, (-c[4][1],-c[4][2],-c[4][3])])
        elseif length(c[4]) == 2
            push!(connections, [c[2], size(positions,1), con_strength, (-c[4][1],-c[4][2])])
        else
            push!(connections, [c[2], size(positions,1), con_strength, (-c[4])])
        end
        # add to treated
        push!(connectionsTreated, c)
    end
    # build a new unitcell
    uc_X = Unitcell(unitcell.lattice_vectors, positions, connections, filename_new)
    # return the new unitcell
    return uc_X
end
export getTransformedUnitcellBondToSite




#-----------------------------------------------------------------------------------------------------------------------------
#
#   Build a the squareroot version of the unitcell
#   Procedure replaces all fully connected plaquettes by sites and then sets the interaction strengths (TODO)
#
#-----------------------------------------------------------------------------------------------------------------------------
function getTransformedUnitcellFCPToSite(unitcell::Unitcell; printFCP=false, connection_strength="AUTO")
    # look for all fully connected plaquettes
    fcp_list = []
    connectionlist = getConnectionList(unitcell)
    neutral_wrap = (0,0,0)
    if length(connectionlist[1][1][4]) == 2
        neutral_wrap = (0,0)
    elseif length(connectionlist[1][1][4]) == 1
        neutral_wrap = (0)
    end
    # define a recursive looking method
    function getFCPs(fcp_current)
        # a site with index current_index is asked if it is in a fully connected plaquette with more sits than given by fcp
        # returns the fcps that the site is in (minimum: the fcp given)
        fcp_list_new = Array[]
        # iterate over all neighbors of the current site
        for nc in connectionlist[Int(fcp_current[end][1])]
            # current test site and test wrap
            test_site = Int(nc[2])
            test_wrap = collect(nc[4]) .+ fcp_current[end][2]
            # check, if "step back"
            if test_site == fcp_current[end-1][1] && test_wrap == fcp_current[end-1][2]
                continue
            end
            # no step back, check if this site is connected to all fcp sites that are already in the fcp
            completely_connected = true
            for site_in_fcp in fcp_current
                # check if new site is connected to this site
                connected = false
                for c in connectionlist[test_site]
                    # check if the destination and wrap agree
                    if Int(c[2]) == Int(site_in_fcp[1]) && collect(c[4]) .+ test_wrap == site_in_fcp[2]
                        # found a match
                        connected = true
                        break
                    end
                end
                # if not connected: complete connection broken
                if !connected
                    completely_connected = false
                    break
                end
            end
            # if it is connected to everything: make a new fcp out of it and go deper into the recursion
            if completely_connected
                fcp_new = []
                for fcp_current_element in fcp_current
                    push!(fcp_new, fcp_current_element)
                end
                push!(fcp_new, [test_site, test_wrap])
                # recursion
                fcp_list_recursion = getFCPs(fcp_new)
                # add those to the list
                for fcp in fcp_list_recursion
                    push!(fcp_list_new, fcp)
                end
            end
        end
        # if empty, insert the current fcp and return
        if length(fcp_list_new) == 0
            push!(fcp_list_new, fcp_current)
        end
        # return fcp list
        return fcp_list_new
    end
    # define a function to give connections for a certain fcp
    function getFCPConnections(fcp)
        # list
        fcp_connections = []
        # find all connections
        for site1 in fcp
        for site2 in fcp
            # add connection from 1 to 2 if found
            for c in connectionlist[Int(site1[1])]
                if c[1] == Int(site1[1]) && c[2] == Int(site2[1]) && collect(c[4]) == site2[2] .- site1[2]
                    push!(fcp_connections, c)
                end
            end
        end
        end
        # return the list of connections
        return fcp_connections
    end
    # define a function to compare to lists
    function sameConnections(list1, list2)
        # check if every element in list 1 is in list 2
        for element1 in list1
            # not found yet
            found = false
            # check all elements of list2
            for element2 in list2
                # check if equal
                if Int(element1[1]) == Int(element2[1]) && Int(element1[2]) == Int(element2[2]) && element1[4] == element2[4]
                    found = true
                    break
                end
            end
            # if not found, return false
            if !found
                return false
            end
        end
        # check if every element in list 2 is in list 1
        for element2 in list2
            # not found yet
            found = false
            # check all elements of list1
            for element1 in list1
                # check if equal
                if Int(element1[1]) == Int(element2[1]) && Int(element1[2]) == Int(element2[2]) && element1[4] == element2[4]
                    found = true
                    break
                end
            end
            # if not found, return false
            if !found
                return false
            end
        end
        # every test was successfull: 2 lists identical
        return true
    end
    # define a function to get the center of the fcp
    function getFCPCenter(fcp)
        # add all sites
        center = unitcell.basis[Int(fcp[1][1])]
        for i in 2:size(fcp,1)
            center = center .+ unitcell.basis[Int(fcp[i][1])]
        end
        # add all wraps
        for i in 1:size(fcp, 1)
            for l in 1:size(unitcell.lattice_vectors,1)
                center = center .+ (fcp[i][2][l] .* unitcell.lattice_vectors[l])
            end
        end
        # devide by number of sites
        center = center ./ size(fcp,1)
    end
    # iterate for all sites
    for s in 1:size(unitcell.basis,1)
        # looking for fully connected plaquettes of site s
        # check if neighbors are part of bigger plaquettes
        for neighbor_connection in connectionlist[s]
            # ask if the neighbor is part of an FCP
            fcp_new = [
                [s, collect(neutral_wrap)],
                [Int(neighbor_connection[2]), collect(neighbor_connection[4])]
            ]
            fcp_list_of_neighbor_and_site = getFCPs(fcp_new)
            # go through all fcps that have been identified and add new ones if necessary
            for fcp in fcp_list_of_neighbor_and_site
                # just add to global list, filtering later on
                push!(fcp_list, fcp)
            end
        end
    end
    # up to here: found all fcps
    fcp_list_raw = fcp_list
    fcp_list = []
    fcp_list_connections = []
    # now: filter the identical fcps out
    for fcp_test in fcp_list_raw
        # boolean if found in the ready to use list
        found = false
        # get the connections
        fcp_test_connections = getFCPConnections(fcp_test)
        # check all elements of the ready to use list
        for (i,fcp) in enumerate(fcp_list)
            # first test: dimensions identical?
            if size(fcp, 1) != size(fcp_test, 1)
                continue
            end
            # second test: get the connections associated with this certain fcp
            fcp_connections = fcp_list_connections[i]
            # test if the same lists, i.e. if both are in the other
            if sameConnections(fcp_connections, fcp_test_connections)
                found = true
                break
            end
        end
        # if not found, add to the list
        if !found
            push!(fcp_list, fcp_test)
            push!(fcp_list_connections, fcp_test_connections)
        end
    end
    # maybe print the FCP list
    if printFCP
        println("In total: $(size(fcp_list,1)) FCPs found")
        for i in 1:size(fcp_list,1)
            println("FCP $(i) with $(size(fcp_list[i], 1)) sites:")
            println(" - sites: $(fcp_list[i])")
            println(" - connections: $(fcp_list_connections[i])")
        end
    end
    # replace all fcps by inserting new sites
    positions_new = []
    for p in unitcell.basis
        push!(positions_new, p)
    end
    # push new sites
    for fcp in fcp_list
        # get the center of the fcp and add it as a new position
        push!(positions_new, getFCPCenter(fcp))
    end
    # define new connections
    connections_new = []
    # iterate over all sites / all fcps
    for (i, fcp) in enumerate(fcp_list)
        # new position index
        pos_index = size(unitcell.basis, 1) + i
        # define new connections for every site in the FCP (new position wrap is neutral wrap)
        for (s,site) in enumerate(fcp)
            # build the fields of the new 2 connections to insert
            con_field_1 = Int(site[1])
            con_field_2 = pos_index
            if typeof(connection_strength) == String && connection_strength == "AUTO"
                con_field_3 = "fcp$(i)con$(s)"
            else
                con_field_3 = connection_strength
            end
            con_field_4_array = collect(neutral_wrap) .- site[2]
            if length(con_field_4_array) == 3
                con_field_4_1 = (con_field_4_array[1], con_field_4_array[2], con_field_4_array[3])
                con_field_4_2 = (-con_field_4_array[1], -con_field_4_array[2], -con_field_4_array[3])
            elseif length(con_field_4_array) == 2
                con_field_4_1 = (con_field_4_array[1], con_field_4_array[2])
                con_field_4_2 = (-con_field_4_array[1], -con_field_4_array[2])
            elseif length(con_field_4_array) == 1
                con_field_4_1 = (con_field_4_array[1])
                con_field_4_2 = (-con_field_4_array[1])
            else
                println("wrong length for field 4: $(con_field_4_array)")
            end
            # build new connections
            con_1 = [con_field_1; con_field_2; con_field_3; con_field_4_1]
            con_2 = [con_field_2; con_field_1; con_field_3; con_field_4_2]
            push!(connections_new, con_1)
            push!(connections_new, con_2)
        end
    end
    # build new unitcell
    unitcell_fcp_replaced = Unitcell(
        unitcell.lattice_vectors,
        positions_new,
        connections_new,
        replace(unitcell.filename, ".jld", "_fcp_to_site.jld"))
    # return the new unitcell
    return unitcell_fcp_replaced
end
export getTransformedUnitcellFCPToSite

#function getTransformedLatticeFCPToSite(lattice::Lattice)
#
#end
#export getTransformedLatticeFCPToSite



#-----------------------------------------------------------------------------------------------------------------------------
#
#   Square the lattice / unitcell
#   This maps all connections to NNN connections by multiplication
#
#-----------------------------------------------------------------------------------------------------------------------------
function getSquaredLattice(lattice::Lattice)
    # get the connectivity matrix of the lattice
    connectivity = getConnectionList(lattice)
    # define a list of new connections
    connections_new = Array[]
    # iterate over all sites and check if they host a NNN connection
    for i in 1:size(lattice.positions,1)
        # check the neighbors of site i
        for (i1,c1) in enumerate(connectivity[i])
        for (i2,c2) in enumerate(connectivity[i])
            # propose a connection c1+c2 if they are different
            if i1 < i2
                # build the offsets
                if length(c1[4]) == 1
                    off_1 = (c2[4] - c1[4])
                    off_2 = (c1[4] - c2[4])
                elseif length(c1[4]) == 2
                    off_1 = (c2[4][1] - c1[4][1], c2[4][2] - c1[4][2])
                    off_2 = (c1[4][1] - c2[4][1], c1[4][2] - c2[4][2])
                elseif length(c1[4]) == 3
                    off_1 = (c2[4][1] - c1[4][1], c2[4][2] - c1[4][2], c2[4][3] - c1[4][3])
                    off_2 = (c1[4][1] - c2[4][1], c1[4][2] - c2[4][2], c1[4][3] - c2[4][3])
                end
                # build two new connections
                if typeof(c1[3]) == String || typeof(c2[3]) == String
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); "$(c1[3])*$(c2[3])"; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); "$(c1[3])*$(c2[3])"; off_2]
                else
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); (c1[3])*(c2[3]); off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); (c1[3])*(c2[3]); off_2]
                end
                # push them to the add list
                push!(connections_new, connection_new_1)
                push!(connections_new, connection_new_2)
            end
        end
        end
    end
    # create new Lattice with only new connections (and leave out old ones)
    return getLatticeWithOptimizedConnections(
        Lattice(
            lattice.unitcell,
            lattice.unitcellRepetitions,
            lattice.lattice_vectors,
            lattice.positions,
            lattice.positions_indices,
            connections_new,
            replace(lattice.filename, ".jld", "_squared.jld"))
    )
end
export getSquaredLattice

function getSquaredUnitcell(unitcell::Unitcell)
    # get the connectivity matrix of the unitcell
    connectivity = getConnectionList(unitcell)
    # define a list of new connections
    connections_new = Array[]
    # iterate over all sites and check if they host a NNN connection
    for i in 1:size(unitcell.basis,1)
        # check the neighbors of site i
        for (i1,c1) in enumerate(connectivity[i])
        for (i2,c2) in enumerate(connectivity[i])
            # propose a connection c1+c2 if the connections c1 and c2 are different
            if i1 < i2
                # build the offsets
                if length(c1[4]) == 1
                    off_1 = (c2[4] - c1[4])
                    off_2 = (c1[4] - c2[4])
                elseif length(c1[4]) == 2
                    off_1 = (c2[4][1] - c1[4][1], c2[4][2] - c1[4][2])
                    off_2 = (c1[4][1] - c2[4][1], c1[4][2] - c2[4][2])
                elseif length(c1[4]) == 3
                    off_1 = (c2[4][1] - c1[4][1], c2[4][2] - c1[4][2], c2[4][3] - c1[4][3])
                    off_2 = (c1[4][1] - c2[4][1], c1[4][2] - c2[4][2], c1[4][3] - c2[4][3])
                end
                # build two new connections
                if typeof(c1[3]) == String || typeof(c2[3]) == String
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); "$(c1[3])*$(c2[3])"; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); "$(c1[3])*$(c2[3])"; off_2]
                else
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); (c1[3])*(c2[3]); off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); (c1[3])*(c2[3]); off_2]
                end
                # push them to the add list
                push!(connections_new, connection_new_1)
                push!(connections_new, connection_new_2)
            end
        end
        end
    end
    # create new Lattice with only new connections (and leave out old ones)
    return getUnitcellWithOptimizedConnections(
        Unitcell(
            unitcell.lattice_vectors,
            unitcell.basis,
            connections_new,
            replace(unitcell.filename, ".jld", "_squared.jld"))
    )
end
export getSquaredUnitcell






#-----------------------------------------------------------------------------------------------------------------------------
#
#   Get a list of all independent sublattices / sublattices of a lattice / unitcell (i.e. not connected parts of a lattice)
#   Not connected parts can occur when "squaring" a bipartite lattice
#
#-----------------------------------------------------------------------------------------------------------------------------
function getIndependentSublattices(lattice::Lattice)
    # save relevant lattice data
    connections = lattice.connections
    positions = lattice.positions
    positions_indices = lattice.positions_indices
    # generate a list for the sublattices
    sublattices = []
    # define an array of labels
    labels = zeros(Int64, size(lattice.positions,1))
    currentlabel = 1
    # get the connectivity list of the lattice
    connectivity = getConnectionList(lattice)
    # iterate over all positions
    for i in 1:size(lattice.positions, 1)
        # get the label of site i
        if labels[i] == 0
            labels[i] = currentlabel
            currentlabel = currentlabel + 1
        end
        # iterate over all connections
        for c in connectivity[i]
            # get the connected site
            j = Int(c[2])
            # check the label of j
            if labels[j] == 0
                # set as connected to i
                labels[j] = labels[i]
            elseif labels[j] == labels[i]
                # just ignore
            else
                # merging of two clusters here!
                newlabel = min(labels[i], labels[j])
                # scan through the whole lattice
                for n in 1:length(labels)
                    if labels[n] == labels[i] || labels[n] == labels[j]
                        labels[n] = newlabel
                    end
                end
                # get a new currentlabel
                currentlabel = Int(maximum(labels)) + 1
            end
        end
    end
    # get a list of labels
    labellist = []
    for l in labels
        if !(l in labellist)
            push!(labellist, l)
        end
    end
    # print how many sublattices found
    println("$(length(labellist)) independent sublattice(s) found")
    # check how many sublattices
    if length(labellist) == 1
        # just push the original lattice into the sublattices list
        sublattice = Lattice(
            lattice.unitcell,
            lattice.unitcellRepetitions,
            lattice.lattice_vectors,
            lattice.positions,
            lattice.positions_indices,
            lattice.connections,
            lattice.filename_output
        )
        push!(sublattices, sublattice)
        # return
        return sublattices
    end
    # iterate over labels
    for l in labellist
        # find the new positions
        positions_new = Array[]
        positions_indices_new = []
        connections_new = Array[]
        mapping_indices_old = []
        # iterate over all connections to identify relevant positions
        for c in connections
            # check if belonging to the sublattice
            if labels[Int(c[1])] != l
                # not of this sublattice
                continue
            end
            # get the starting and finishing position
            index_from = Int(c[1])
            index_to = Int(c[2])
            # check where they map to
            if index_from in mapping_indices_old
                index_from_new = findfirst(mapping_indices_old, index_from)
            else
                # add the position to the list of positions
                push!(positions_new, positions[index_from])
                push!(positions_indices_new, positions_indices[index_from])
                # add the index mapping to the list
                push!(mapping_indices_old, index_from)
                # get the new index
                index_from_new = size(positions_new,1)
            end
            if index_to in mapping_indices_old
                index_to_new = findfirst(mapping_indices_old, index_to)
            else
                # add the position to the list of positions
                push!(positions_new, positions[index_to])
                push!(positions_indices_new, positions_indices[index_to])
                # add the index mapping to the list
                push!(mapping_indices_old, index_to)
                # get the new index
                index_to_new = size(positions_new,1)
            end
            # push the new connection
            # format [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            connection_new = [index_from_new; index_to_new; c[3]; c[4]]
			# register as connection
			push!(connections_new, connection_new)
        end
        # generate a new sublattice
        sublattice = Lattice(
            lattice.unitcell,
            lattice.unitcellRepetitions,
            lattice.lattice_vectors,
            positions_new,
            positions_indices_new,
            connections_new,
            "$(lattice.filename[1:end-4])_sublattice_$(l).jld"
        )
        # push it the list
        push!(sublattices, sublattice)
    end
    # return the list
    return sublattices
end
export getIndependentSublattices

function getIndependentSubunitcells(unitcell::Unitcell)
    # save relevant lattice data
    connections = unitcell.connections
    positions = unitcell.basis
    # generate a list for the sublattices
    subunitcells = []
    # define an array of labels
    labels = zeros(Int64, size(positions,1))
    currentlabel = 1
    # get the connectivity list of the unitcell
    connectivity = getConnectionList(unitcell)
    # iterate over all positions
    for i in 1:size(positions, 1)
        # get the label of site i
        if labels[i] == 0
            labels[i] = currentlabel
            currentlabel = currentlabel + 1
        end
        # iterate over all connections
        for c in connectivity[i]
            # get the connected site
            j = Int(c[2])
            # check the label of j
            if labels[j] == 0
                # set as connected to i
                labels[j] = labels[i]
            elseif labels[j] == labels[i]
                # just ignore
            else
                # merging of two clusters here!
                newlabel = min(labels[i], labels[j])
                # scan through the whole lattice
                for n in 1:length(labels)
                    if labels[n] == labels[i] || labels[n] == labels[j]
                        labels[n] = newlabel
                    end
                end
                # get a new currentlabel
                currentlabel = Int(maximum(labels)) + 1
            end
        end
    end
    # get a list of labels
    labellist = []
    for l in labels
        if !(l in labellist)
            push!(labellist, l)
        end
    end
    # print how many sublattices found
    println("$(length(labellist)) independent sublattice(s) found")
    # check how many sublattices
    if length(labellist) == 1
        # just push the original lattice into the sublattices list
        subunitcell = Unitcell(
            unitcell.lattice_vectors,
            unitcell.basis,
            unitcell.connections,
            unitcell.filename
        )
        push!(subunitcells, subunitcell)
        # return
        return subunitcells
    end
    # iterate over labels
    for l in labellist
        # find the new positions
        positions_new = Array[]
        connections_new = Array[]
        mapping_indices_old = []
        # iterate over all connections to identify relevant positions
        for c in connections
            # check if belonging to the sublattice
            if labels[Int(c[1])] != l
                # not of this sublattice
                continue
            end
            # get the starting and finishing position
            index_from = Int(c[1])
            index_to = Int(c[2])
            # check where they map to
            if index_from in mapping_indices_old
                index_from_new = findfirst(mapping_indices_old, index_from)
            else
                # add the position to the list of positions
                push!(positions_new, positions[index_from])
                # add the index mapping to the list
                push!(mapping_indices_old, index_from)
                # get the new index
                index_from_new = size(positions_new,1)
            end
            if index_to in mapping_indices_old
                index_to_new = findfirst(mapping_indices_old, index_to)
            else
                # add the position to the list of positions
                push!(positions_new, positions[index_to])
                # add the index mapping to the list
                push!(mapping_indices_old, index_to)
                # get the new index
                index_to_new = size(positions_new,1)
            end
            # push the new connection
            # format [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            connection_new = [index_from_new; index_to_new; c[3]; c[4]]
			# register as connection
			push!(connections_new, connection_new)
        end
        # generate a new sublattice
        subunitcell = Unitcell(
            unitcell.lattice_vectors,
            positions_new,
            connections_new,
            replace(unitcell.filename, ".jld", "_subunitcell_$(l).jld")
        )
        # push it the list
        push!(subunitcells, subunitcell)
    end
    # return the list
    return subunitcells
end
export getIndependentSubunitcells






#-----------------------------------------------------------------------------------------------------------------------------
#
#   Optimize the lattice / unitcell by collapsing redundant connections into fewer connections and return a new object
#
#-----------------------------------------------------------------------------------------------------------------------------
function getLatticeWithOptimizedConnections(lattice::Lattice)
    # build up new connections
    connections_new = Array[]
    # go through all old connections
    for c in lattice.connections
        # check if the connection already is present in the list
        found = false
        for c_new in connections_new
            if c[1] == c_new[1] && c[2] == c_new[2] && c[4] == c_new[4]
                # found a redundent connection
                found = true
                # add connection strength
                if typeof(c[3]) == String || typeof(c_new[3]) == String
                    c_new[3] = "$(c_new[3])+$(c[3])"
                else
                    c_new[3] = c_new[3] + c[3]
                end
                # break the inner loop
                break
            end
        end
        # if not, add the current connection
        if !found
            push!(connections_new, c)
        end
    end
    # check if the connection strength is 0
    connections_new_2 = Array[]
    for c in connections_new
        # check if close to 0 connection strength
        if typeof(c[3]) == Float64 || typeof(c[3]) == Int64
            if abs(c[3]) > 1e-18
                push!(connections_new_2, c)
            else
                # do not push! c[3] too small!
            end
        elseif isnull(tryparse(Float64,c[3]))
            # not a float and cannot be parsed to float, just push in list as it is
            push!(connections_new_2, c)
        else
            # check if the strength is non-zero
            if abs(parse(Float64,c[3])) > 1e-18
                push!(connections_new_2, c)
            else
                # do not push! c[3] too small!
            end
        end
    end
    # return a new lattice
    return Lattice(
        lattice.unitcell,
        lattice.unitcellRepetitions,
        lattice.lattice_vectors,
        lattice.positions,
        lattice.positions_indices,
        connections_new_2,
        lattice.filename)
end
export getLatticeWithOptimizedConnections

function getUnitcellWithOptimizedConnections(unitcell::Unitcell)
    # build up new connections
    connections_new = Array[]
    # go through all old connections
    for c in unitcell.connections
        # check if the connection already is present in the list
        found = false
        for c_new in connections_new
            if c[1] == c_new[1] && c[2] == c_new[2] && c[4] == c_new[4]
                # found a redundent connection
                found = true
                # add connection strength
                if typeof(c[3]) == String || typeof(c_new[3]) == String
                    c_new[3] = "$(c_new[3])+$(c[3])"
                else
                    c_new[3] = c_new[3] + c[3]
                end
                # break the inner loop
                break
            end
        end
        # if not, add the current connection
        if !found
            push!(connections_new, c)
        end
    end
    # check if the connection strength is 0
    connections_new_2 = Array[]
    for c in connections_new
        # check if close to 0 connection strength
        if typeof(c[3]) == Float64 || typeof(c[3]) == Int64
            if abs(c[3]) > 1e-18
                push!(connections_new_2, c)
            else
                # do not push! c[3] too small!
            end
        elseif isnull(tryparse(Float64,c[3]))
            # not a float and cannot be parsed to float, just push in list as it is
            push!(connections_new_2, c)
        else
            # check if the strength is non-zero
            if abs(parse(Float64,c[3])) > 1e-18
                push!(connections_new_2, c)
            else
                # do not push! c[3] too small!
            end
        end
    end
    # return a new unitcell
    return Unitcell(
        unitcell.lattice_vectors,
        unitcell.basis,
        connections_new_2,
        unitcell.filename)
end
export getUnitcellWithOptimizedConnections




#-----------------------------------------------------------------------------------------------------------------------------
#
#   Optimize the lattice / unitcell by collapsing redundant connections into fewer connections and modify given object
#
#-----------------------------------------------------------------------------------------------------------------------------
function optimizeConnections!(lattice::Lattice)
    # build up new connections
    connections_new = Array[]
    # go through all old connections
    for c in lattice.connections
        # check if the connection already is present in the list
        found = false
        for c_new in connections_new
            if c[1] == c_new[1] && c[2] == c_new[2] && c[4] == c_new[4]
                # found a redundent connection
                found = true
                # add connection strength
                if typeof(c[3]) == String || typeof(c_new[3]) == String
                    c_new[3] = "$(c_new[3])+$(c[3])"
                else
                    c_new[3] = c_new[3] + c[3]
                end
                # break the inner loop
                break
            end
        end
        # if not, add the current connection
        if !found
            push!(connections_new, c)
        end
    end
    # check if the connection strength is 0
    connections_new_2 = Array[]
    for c in connections_new
        # check if close to 0 connection strength
        if typeof(c[3]) == Float64 || typeof(c[3]) == Int64
            if abs(c[3]) > 1e-18
                push!(connections_new_2, c)
            else
                # do not push! c[3] too small!
            end
        elseif isnull(tryparse(Float64,c[3]))
            # not a float and cannot be parsed to float, just push in list as it is
            push!(connections_new_2, c)
        else
            # check if the strength is non-zero
            if abs(parse(Float64,c[3])) > 1e-18
                push!(connections_new_2, c)
            else
                # do not push! c[3] too small!
            end
        end
    end
    # overwrite the connections in the given lattice
    lattice.connections = connections_new_2
end
function optimizeConnections!(unitcell::Unitcell)
    # build up new connections
    connections_new = Array[]
    # go through all old connections
    for c in unitcell.connections
        # check if the connection already is present in the list
        found = false
        for c_new in connections_new
            if c[1] == c_new[1] && c[2] == c_new[2] && c[4] == c_new[4]
                # found a redundent connection
                found = true
                # add connection strength
                if typeof(c[3]) == String || typeof(c_new[3]) == String
                    c_new[3] = "$(c_new[3])+$(c[3])"
                else
                    c_new[3] = c_new[3] + c[3]
                end
                # break the inner loop
                break
            end
        end
        # if not, add the current connection
        if !found
            push!(connections_new, c)
        end
    end
    # check if the connection strength is 0
    connections_new_2 = Array[]
    for c in connections_new
        # check if close to 0 connection strength
        if typeof(c[3]) == Float64 || typeof(c[3]) == Int64
            if abs(c[3]) > 1e-18
                push!(connections_new_2, c)
            else
                # do not push! c[3] too small!
            end
        elseif isnull(tryparse(Float64,c[3]))
            # not a float and cannot be parsed to float, just push in list as it is
            push!(connections_new_2, c)
        else
            # check if the strength is non-zero
            if abs(parse(Float64,c[3])) > 1e-18
                push!(connections_new_2, c)
            else
                # do not push! c[3] too small!
            end
        end
    end
    # overwrite the connections in the given lattice
    unitcell.connections = connections_new_2
end
export optimizeConnections!








#-----------------------------------------------------------------------------------------------------------------------------
#
#   Set / Map interaction strengths of the unitcell / lattice
#
#-----------------------------------------------------------------------------------------------------------------------------
function setAllInteractionStrengths!(unitcell::Unitcell, strengthNew)
    # go through all connections and modify their strength
    for c in unitcell.connections
        c[3] = strengthNew
    end
end
function setAllInteractionStrengths!(lattice::Lattice, strengthNew)
    # go through all connections and modify their strength
    for c in lattice.connections
        c[3] = strengthNew
    end
end
export setAllInteractionStrengths!



function mapInteractionStrengths!(unitcell::Unitcell, mapping; replace_in_strings=true, evaluate=false)
    # get the old strengths
    old_strengths = keys(mapping)
    # iterate for replacement
    for old_strength in old_strengths
        # new strength is given by mapping
        new_strength = mapping[old_strength]
        # check all connections that are identical
        for c in unitcell.connections
            if c[3] == old_strength
                c[3] = new_strength
            end
        end
        # check for string replacement
        if replace_in_strings
        for c in unitcell.connections
            if typeof(c[3]) == String && contains(c[3], old_strength)
                c[3] = replace(c[3], old_strength, new_strength)
            end
        end
        end
    end
    # maybe even evaluate
    if evaluate
        for c in unitcell.connections
            if typeof(c[3]) == String
                c[3] = eval(parse(c[3]))
            end
        end
    end
end
function mapInteractionStrengths!(lattice::Lattice, mapping; replace_in_strings=true, evaluate=false)
    # get the old strengths
    old_strengths = keys(mapping)
    # iterate for replacement
    for old_strength in old_strengths
        # new strength is given by mapping
        new_strength = mapping[old_strength]
        # check all connections that are identical
        for c in lattice.connections
            if c[3] == old_strength
                c[3] = new_strength
            end
        end
        # check for string replacement
        if replace_in_strings
        for c in lattice.connections
            if typeof(c[3]) == String && contains(c[3], old_strength)
                c[3] = replace(c[3], old_strength, new_strength)
            end
        end
        end
    end
    # maybe even evaluate
    if evaluate
        for c in lattice.connections
            if typeof(c[3]) == String
                c[3] = eval(parse(c[3]))
            end
        end
    end
end
export mapInteractionStrengths!




#-----------------------------------------------------------------------------------------------------------------------------
#
#   Remove interactions with given strengths of the unitcell / lattice
#
#-----------------------------------------------------------------------------------------------------------------------------
function removeConnections!(lattice::Lattice, strengthToRemove=0.0)
    connections_new = Array[]
    for c in lattice.connections
        if c[3] != strengthToRemove
            push!(connections_new, c)
        end
    end
    lattice.connections = connections_new;
end
function removeConnections!(unitcell::Unitcell, strengthToRemove=0.0)
    connections_new = Array[]
    for c in unitcell.connections
        if c[3] != strengthToRemove
            push!(connections_new, c)
        end
    end
    unitcell.connections = connections_new;
end

export removeConnections!



# TODO remove sites



#-----------------------------------------------------------------------------------------------------------------------------
#
#   Add interactions along next-nearest neighbors
#
#-----------------------------------------------------------------------------------------------------------------------------
function addNextNearestNeighborsToConnections!(lattice::Lattice; strengthNNN="AUTO", connectionsNN=["ALL"])
    # get the connectivity matrix of the lattice
    connectivity = getConnectionList(lattice)
    # define a list of new connections
    connections_new = Array[]
    # iterate over all sites and check if they host a NNN connection
    for i in 1:size(lattice.positions,1)
        # check the neighbors of site i
        for (i1,c1) in enumerate(connectivity[i])
        for (i2,c2) in enumerate(connectivity[i])
            # check if the connections are valid
            if !(in("ALL",connectionsNN) || (in(c1[3],connectionsNN) && in(c2[3],connectionsNN)))
                continue
            end
            # propose a connection c1+c2 if they are different
            if i1 < i2
                # build the offsets
                if length(c1[4]) == 1
                    off_1 = (c2[4] - c1[4])
                    off_2 = (c1[4] - c2[4])
                elseif length(c1[4]) == 2
                    off_1 = (c2[4][1] - c1[4][1], c2[4][2] - c1[4][2])
                    off_2 = (c1[4][1] - c2[4][1], c1[4][2] - c2[4][2])
                elseif length(c1[4]) == 3
                    off_1 = (c2[4][1] - c1[4][1], c2[4][2] - c1[4][2], c2[4][3] - c1[4][3])
                    off_2 = (c1[4][1] - c2[4][1], c1[4][2] - c2[4][2], c1[4][3] - c2[4][3])
                end
                # build two new connections
                if strengthNNN == "AUTO"
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); "NNN$(c1[2])to$(c2[2])"; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); "NNN$(c1[2])to$(c2[2])"; off_2]
                elseif strengthNNN == "MULTIPLY"
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); "$(c1[3])*$(c2[3])"; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); "$(c1[3])*$(c2[3])"; off_2]
                else
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); strengthNNN; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); strengthNNN; off_2]
                end
                # push them to the add list
                push!(connections_new, connection_new_1)
                push!(connections_new, connection_new_2)
            end
        end
        end
    end
    # push all connections to the lattice
    for c in connections_new
        push!(lattice.connections, c)
    end
end
function addNextNearestNeighborsToConnections!(unitcell::Unitcell; strengthNNN="AUTO", connectionsNN=["ALL"])
    # get the connectivity matrix of the unitcell
    connectivity = getConnectionList(unitcell)
    # define a list of new connections
    connections_new = Array[]
    # iterate over all sites and check if they host a NNN connection
    for i in 1:size(unitcell.basis,1)
        # check the neighbors of site i
        for (i1,c1) in enumerate(connectivity[i])
        for (i2,c2) in enumerate(connectivity[i])
            # check if the connections are valid
            if !(in("ALL",connectionsNN) || (in(c1[3],connectionsNN) && in(c2[3],connectionsNN)))
                continue
            end
            # propose a connection c1+c2 if they are different
            if i1 < i2
                # build the offsets
                if length(c1[4]) == 1
                    off_1 = (c2[4] - c1[4])
                    off_2 = (c1[4] - c2[4])
                elseif length(c1[4]) == 2
                    off_1 = (c2[4][1] - c1[4][1], c2[4][2] - c1[4][2])
                    off_2 = (c1[4][1] - c2[4][1], c1[4][2] - c2[4][2])
                elseif length(c1[4]) == 3
                    off_1 = (c2[4][1] - c1[4][1], c2[4][2] - c1[4][2], c2[4][3] - c1[4][3])
                    off_2 = (c1[4][1] - c2[4][1], c1[4][2] - c2[4][2], c1[4][3] - c2[4][3])
                end
                # build two new connections
                if strengthNNN == "AUTO"
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); "NNN$(c1[2])to$(c2[2])"; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); "NNN$(c1[2])to$(c2[2])"; off_2]
                elseif strengthNNN == "MULTIPLY"
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); "$(c1[3])*$(c2[3])"; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); "$(c1[3])*$(c2[3])"; off_2]
                else
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); strengthNNN; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); strengthNNN; off_2]
                end
                # push them to the add list
                push!(connections_new, connection_new_1)
                push!(connections_new, connection_new_2)
            end
        end
        end
    end
    # push all connections to the unitcell
    for c in connections_new
        push!(unitcell.connections, c)
    end
end

export addNextNearestNeighborsToConnections!
