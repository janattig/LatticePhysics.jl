################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT LATTICE CONSTRUCTION RELATED STUFF
#
#   STRUCTURE OF THE FILE
#
#   1) CONSTRUCTION BASED ON UNITCELLS (2D & 3D)
#   - Periodic Boundaries
#   - Open Boundaries
#   - Semi-Periodic Boundaries
#   - General Boundaries
#
#   2) CONSTRUCTION BASED ON BOND DISTANCE (2D & 3D)
#
#   3) CONSTRUCTION BASED ON SHAPE (2D & 3D)
#   - General Shape
#   - Sphere
#   - Box
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
        filename = "$(FOLDER_LATTICES)$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
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
	positions = Array{Float64, 1}[]
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
	connections = Array{Any,1}[]

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
    lattice_vectors = Array{Float64, 1}[]
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
        filename = "$(FOLDER_LATTICES)$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
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
	positions = Array{Float64, 1}[]
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
	connections = Array{Any,1}[]

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
    lattice_vectors = Array{Float64, 1}[]
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
julia> lattice = getLatticePeriodic(unitcell, [10, 20])
LatticePhysics.Lattice(...)

julia> lattice = getLatticePeriodic(unitcell, [10, 20], load=true)
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
function getLatticeOpen2D(unitcell::Unitcell, repetition_array::Array{Int64}; save::Bool=false, load::Bool=false)

    # extract the cardinal directions of the lattice from the array
    N_a1 = abs(repetition_array[1])
    N_a2 = abs(repetition_array[2])

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_open_$(N_a1)_$(N_a2).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "$(FOLDER_LATTICES)$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
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
	positions = Array{Float64, 1}[]
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
	connections = Array{Any,1}[]

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
    lattice_vectors = Array{Float64, 1}[]

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
function getLatticeOpen3D(unitcell::Unitcell, repetition_array::Array{Int64}; save::Bool=false, load::Bool=false)

    # extract the cardinal directions of the lattice from the array
    N_a1 = abs(repetition_array[1])
    N_a2 = abs(repetition_array[2])
    N_a3 = abs(repetition_array[3])

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_open_$(N_a1)_$(N_a2)_$(N_a3).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "$(FOLDER_LATTICES)$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
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
	positions = Array{Float64, 1}[]
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
		positions[index(i,j,k,a)] = uc_basis[a] + i*uc_lattice_vectors[1] + j*uc_lattice_vectors[2] + k*uc_lattice_vectors[3]
        positions_indices[index(i,j,k,a)] = a
	end
	end
    end
    end

    # GENERATE NEW CONNECTIONS
	connections = Array{Any,1}[]

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
    lattice_vectors = Array{Float64, 1}[]

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
    getLatticeOpen(unitcell::Unitcell, repetition_array::Array{Int64} [; save::Bool, load::Bool])

Function to construct a finite lattice with open boundary conditions
(bonds are cut when connecting to sites outside the wanted area) on all sides out of a `Unitcell` object.
The number of unitcells that are put together in each elementery direction is passed in the
`repetition_array`.

Additionally, the newly created `Lattice` object can directly be saved. If this has been done before, passing a `load=true`
will allow to load the object instead of creating it again.

Note that this function works for both 2D and 3D unitcells.



# Examples

```julia-repl
julia> lattice = getLatticeOpen(unitcell, [10, 20])
LatticePhysics.Lattice(...)

julia> lattice = getLatticeOpen(unitcell, [10, 20], load=true)
LatticePhysics.Lattice(...)
```
"""
function getLatticeOpen(unitcell::Unitcell, repetition_array::Array{Int64}; save::Bool=false, load::Bool=false)

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
function getLatticeSemiperiodic2D(unitcell::Unitcell, repetition_array::Array{Int64}; save::Bool=false, load::Bool=false)

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
        filename = "$(FOLDER_LATTICES)$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
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
	positions = Array{Float64, 1}[]
    positions_indices = Int64[]

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
	connections = Array{Any,1}[]

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
    lattice_vectors = Array{Float64, 1}[]
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
function getLatticeSemiperiodic3D(unitcell::Unitcell, repetition_array::Array{Int64}; save::Bool=false, load::Bool=false)

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
        filename = "$(FOLDER_LATTICES)$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
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
	positions = Array{Float64, 1}[]
	positions_indices = Int64[]

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
	connections = Array{Any,1}[]

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
    lattice_vectors = Array{Float64, 1}[]
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




"""
    getLatticeSemiperiodic(unitcell::Unitcell, repetition_array::Array{Int64} [; save::Bool, load::Bool])

Function to construct a finite lattice with a mix of open and periodic boundary conditions out of a `Unitcell` object.
The number of unitcells that are put together in each elementery direction is passed in the
`repetition_array`, where negative entries denote periodic boundary conditions whereas positive entries denote open boundary conditions.
Note that this function can only handle a mix of periodic and open but not purely open or periodic boundary conditions.

Also, the newly created `Lattice` object can directly be saved. If this has been done before, passing a `load=true`
will allow to load the object instead of creating it again.

Note that this function works for both 2D and 3D unitcells.



# Examples

```julia-repl
julia> lattice = getLatticeSemiperiodic(unitcell, [10, -20])
LatticePhysics.Lattice(...)

julia> lattice = getLatticeSemiperiodic(unitcell, [10, -20], load=true)
LatticePhysics.Lattice(...)
```
"""
function getLatticeSemiperiodic(unitcell::Unitcell, repetition_array::Array{Int64}; save::Bool=false, load::Bool=false)

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





"""
    getLattice(unitcell::Unitcell, repetition_array::Array{Int64} [; save::Bool, load::Bool])
    getLattice(unitcell::Unitcell, repetitions::Int64 [; save::Bool, load::Bool])

Function to construct a finite lattice with a mix of open and periodic boundary conditions out of a `Unitcell` object.
The number of unitcells that are put together in each elementery direction is passed in the
`repetition_array`, where negative entries denote periodic boundary conditions whereas positive entries denote open boundary conditions.
Note that this function CAN handle any mix of periodic and open including purely open or periodic boundary conditions.

Also, the newly created `Lattice` object can directly be saved. If this has been done before, passing a `load=true`
will allow to load the object instead of creating it again.

Note that this function works for both 2D and 3D unitcells.



# Examples

```julia-repl
julia> lattice = getLattice(unitcell, [10, -20])
LatticePhysics.Lattice(...)

julia> lattice = getLattice(unitcell, [10, -20], load=true)
LatticePhysics.Lattice(...)

julia> lattice = getLattice(unitcell, -20, load=true)
LatticePhysics.Lattice(...)
```
"""
function getLattice(unitcell::Unitcell, repetition_array::Array{Int64}; save=false, load=false)

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
function getLattice(unitcell::Unitcell, repetitions::Int64; save=false, load=false)
    # just parse through to more general method
    return getLattice(unitcell, repetitions.*ones(Int64, size(unitcell.lattice_vectors,1)), save=save, load=load)
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
function getLatticeByBondDistance2D(unitcell::Unitcell, bonddistance::Int64; origin::Int64=1, load::Bool=false, save::Bool=false)

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_by_bonddistance_$(bonddistance)_from_$(origin).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "$(FOLDER_LATTICES)$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
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
    connections = Array{Any,1}[]

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
    positions = Array{Float64, 1}[]
    positions_indices = Int64[]

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
            Int64[],
            Array{Float64, 1}[],
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
function getLatticeByBondDistance3D(unitcell::Unitcell, bonddistance::Int64; origin::Int64=1, load::Bool=false, save::Bool=false)

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_by_bonddistance_$(bonddistance)_from_$(origin).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "$(FOLDER_LATTICES)$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
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
    connections = Array{Any, 1}[]

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
    positions = Array{Float64, 1}[]
    positions_indices = Int64[]

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
            Int64[],
            Array{Float64, 1}[],
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
    getLatticeByBondDistance(unitcell::Unitcell, bonddistance::Int64 [; origin::Int64=1, load::Bool=false, save::Bool=false])

Function to construct a finite lattice (with open boundary conditions) as spreading around an origin site up to a certain bond distance. The lattice is
created from a passed `Unitcell` object.

Also, the newly created `Lattice` object can directly be saved. If this has been done before, passing a `load=true`
will allow to load the object instead of creating it again.

Note that this function works for both 2D and 3D unitcells.



# Examples

```julia-repl
julia> lattice = getLatticeByBondDistance(unitcell, 5)
LatticePhysics.Lattice(...)

julia> lattice = getLatticeByBondDistance(unitcell, 8, origin=2)
LatticePhysics.Lattice(...)
```
"""
function getLatticeByBondDistance(unitcell::Unitcell, bonddistance::Int64; origin::Int64=1, load=false, save=false)

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
function getLatticeInShape2D(unitcell::Unitcell, shape::Function, shapename::String="unknown"; origin::Int64=1, load::Bool=false, save::Bool=false)

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_in_shape_$(shapename).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "$(FOLDER_LATTICES)$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
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
    connections = Array{Any,1}[]

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
    positions = Array{Float64, 1}[]
    positions_indices = Int64[]

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
            Int64[],
            Array{Float64, 1}[],
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
function getLatticeInShape3D(unitcell::Unitcell, shape::Function, shapename::String="unknown"; origin::Int64=1, load::Bool=false, save::Bool=false)

    # generate the filename of the output
    if contains(unitcell.filename, FOLDER_UNITCELLS)
        filename = replace(unitcell.filename, FOLDER_UNITCELLS, FOLDER_LATTICES)
        filename = replace(filename, ".jld", "_in_shape_$(shapename).jld")
        filename = replace(filename, "_unitcell_", "_lattice_")
    else
        filename = "$(FOLDER_LATTICES)$(split(unitcell.filename, FOLDER_UNITCELLS[end])[end])"
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
    connections = Array{Any,1}[]

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
    positions = Array{Float64, 1}[]
    positions_indices = Int64[]

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
            Int64[],
            Array{Float64, 1}[],
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
    getLatticeInShape(unitcell::Unitcell, shape::Function, [ shapename::String="unknown" ; origin::Int64=1, load::Bool=false, save::Bool=false])

Function to construct a finite lattice (with open boundary conditions)
as spreading around an origin site up to a certain shape which is given by the passed `Function` `shape`.
The lattice is created from a passed `Unitcell` object.

Also, the newly created `Lattice` object can directly be saved. If this has been done before, passing a `load=true`
will allow to load the object instead of creating it again.

Note that this function works for both 2D and 3D unitcells.



# Examples

```julia-repl
julia> lattice = getLatticeInShape(unitcell, shape_function)
LatticePhysics.Lattice(...)

julia> lattice = getLatticeInShape(unitcell, shape_function, "myshape")
LatticePhysics.Lattice(...)

julia> lattice = getLatticeInShape(unitcell, shape_function, "myshape", origin=2)
LatticePhysics.Lattice(...)
```
"""
function getLatticeInShape(unitcell::Unitcell, shape::Function, shapename::String="unknown"; origin::Int64=1, load::Bool=false, save::Bool=false)

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
export getLatticeInShape






#-----------------------------------------------------------------------------------------------------------------------------
#
#   Special cases for building in shape for 2D and 3D and general lattices
#
#-----------------------------------------------------------------------------------------------------------------------------


"""
    getLatticeInSphere(unitcell::Unitcell, radius::Float64 [; origin::Int64=1, load::Bool=false, save::Bool=false])

Function to construct a finite lattice (with open boundary conditions)
as spreading around an origin site up to a radius (forming a sphere / circle).
The lattice is created from a passed `Unitcell` object.

Also, the newly created `Lattice` object can directly be saved. If this has been done before, passing a `load=true`
will allow to load the object instead of creating it again.

Note that this function works for both 2D and 3D unitcells.



# Examples

```julia-repl
julia> lattice = getLatticeInSphere(unitcell, 3.0)
LatticePhysics.Lattice(...)

julia> lattice = getLatticeInSphere(unitcell, 3.0, origin=2)
LatticePhysics.Lattice(...)
```
"""
function getLatticeInSphere(unitcell::Unitcell, radius::Float64; origin::Int64=1, load::Bool=false, save::Bool=false)

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


"""
    getLatticeInBox(unitcell::Unitcell, extent_array::Array{Float64} [; origin::Int64=1, load::Bool=false, save::Bool=false])

Function to construct a finite lattice (with open boundary conditions)
as spreading around an origin site up to the boundaries of a given box / rectangle.
The lattice is created from a passed `Unitcell` object.

The box dimensions are given by the passed `extent_array` which contains the length in all elementary
directions as `Float64` numbers. The origin site is place exactly at the center of the box.

The newly created `Lattice` object can directly be saved. If this has been done before, passing a `load=true`
will allow to load the object instead of creating it again.

Note that this function works for both 2D and 3D unitcells.



# Examples

```julia-repl
julia> lattice = getLatticeInSphere(unitcell, 3.0)
LatticePhysics.Lattice(...)

julia> lattice = getLatticeInSphere(unitcell, 3.0, origin=2)
LatticePhysics.Lattice(...)
```
"""
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
