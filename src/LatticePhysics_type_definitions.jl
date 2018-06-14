################################################################################
#
#   TYPE DEFINITIONS OF OBJECT CLASSES IN JULIA
#
#   STRUCTURE OF THE FILE:
#
#   1) Definition of UNITCELL OBJECT
#       - Saving
#       - Loading
#
#   2) Definition of LATTICE OBJECT
#       - Saving
#       - Loading
#
#   3) CONVERSION LATTICE --> UNITCELL
#
#   4) OBTAINING GENERAL INFORMATION
#       - Print Information (printInfo(..))
#       - obtain the connection strengths of an object
#
#   5) OBTAINING CONNECTIVITY INFORMATION
#       - getConnectivityList(..)
#       - getConnectionList(..)
#
################################################################################






#-----------------------------------------------------------------------------------------------------------------------------
#
#   THE UNIT CELL OBJECT CLASS IN JULIA
#
#   FIELDS:
#   - lattice_vectors   lattice vectors of the underlying bravais lattice
#   - basis             basis sites of the unit cell
#   - connections       the connetions between sites of the unit cells
#   - filename          The filename that the unitcell should be stored into
#
#-----------------------------------------------------------------------------------------------------------------------------
"""
    mutable struct Unitcell

The type that contains information on a unitcell. Fields are

    lattice_vectors :: Array{Array{Float64, 1}, 1}
    basis           :: Array{Array{Float64, 1}, 1}
    connections     :: Array{Array{Any, 1}, 1}
    filename        :: String

Note that connections have the format

    [index_from::Int, index_to::Int, strength::Any, unitcell_offset::Tuple]


New objects can be created with

    Unitcell(filename::String)

to load from a given file with name `filename` or with

    Unitcell()

if the constructed object serves as a placeholder. In this case, the filename
is set to `UNITCELL_DUMMY_FILENAME`.



# Examples

```julia-repl
julia> Unitcell("myunitcellfile.jld")
```
"""
mutable struct Unitcell

    # FIELDS

    # basis vectors of the Bravais lattice
    lattice_vectors::Array{Array{Float64, 1}, 1}

    # basis sites within the unitcell
    basis::Array{Array{Float64, 1}, 1}

    # list of connections of the format
    # [index_from, index_to, strength, unitcell_offset]
    connections::Array{Array{Any, 1}, 1}

    # filename if saved
    filename::String



    # the general constructor
    function Unitcell(lattice_vectors::Array{Array{Float64, 1}, 1}, basis::Array{Array{Float64, 1}, 1}, connections::Array{Array,1}, filename::String)
        return new(lattice_vectors, basis, connections, filename)
    end
    function Unitcell(lattice_vectors::Array{Array, 1}, basis::Array{Array, 1}, connections::Array{Array,1}, filename::String)
        return new(lattice_vectors, basis, connections, filename)
    end
    # the custom constructor when loading a unitcell from a file
    function Unitcell(filename::String)
        # define an empty unitcell
        uc = new(nothing,nothing,nothing,filename)
        # load the uc
        loadUnitcell(uc)
        # return the uc
        return uc
    end
    # the dummy constructor (for lattices that do not need a UC)
    function Unitcell()
        # just initialize everything
        return new([], [], [], UNITCELL_DUMMY_FILENAME)
    end

end







#-----------------------------------------
#
#   METHODS FOR SAVING UNITCELLS
#
#-----------------------------------------


"""
    saveUnitcell(uc::Unitcell [, filename::String])

Save the `Unitcell` object `unitcell` to a file by using the module `JLD`.
A new filename can be given by passing an additional argument.

Note that the filename has to end with `.jld` or the ending will be added automatically.

Folders that lead to the specified file will be created on the way.



# Examples

```julia-repl
julia> saveUnitcell(unitcell)
"filename.jld"

julia> saveUnitcell(unitcell, "testfile.jld")
"testfile.jld"
```
"""
function saveUnitcell(uc::Unitcell, filename::String="NONE")
    # maybe set the new filename
    if filename != "NONE"
        uc.filename = filename
    end
    # make sure that the file ends with .jld
    if uc.filename[end-4:end] != ".jld"
        uc.filename = "$(filename).jld"
    end
    # ensure default path is build
    if contains(uc.filename, "/")
		# get the containing folder
		folder = uc.filename[1:findlast(uc.filename, '/')]
		# build the path to that folder
		mkpath(folder)
	end
    # open the file and write the fields into it
    save(uc.filename,
        "basis",            uc.basis,
        "connections",      uc.connections,
        "lattice vectors",  uc.lattice_vectors,
        "filename",         uc.filename
    )
    # return the filename
    return uc.filename
end



#-----------------------------------------
#
#   METHODS FOR LOADING UNITCELLS
#
#-----------------------------------------

"""
    loadUnitcell(uc::Unitcell [, filename::String])
    loadUnitcell(filename::String)

Loads a `Unitcell` object `unitcell` from its file by using the module `JLD`.
A new filename can be given by passing an additional argument.
Moreover, if there is no `Unitcell` object which should be loaded but one wants to obtain a new
object, one can only pass the filename of the unitcell file.



# Examples

```julia-repl
julia> loadUnitcell(unitcell)
"filename.jld"

julia> loadUnitcell(unitcell, "testfile.jld")
"testfile.jld"

julia> loadUnitcell("testfile.jld")
LatticePhysics.Unitcell[...]
```
"""
function loadUnitcell(uc::Unitcell, filename::String="NONE")
    # maybe set the new filename
    if filename != "NONE"
        uc.filename = filename
    end
    # open the file and load all fields from it
    uc.basis            = load(uc.filename, "basis")
    uc.connections      = load(uc.filename, "connections")
    uc.lattice_vectors  = load(uc.filename, "lattice vectors")
    # return the filename
    return uc.filename
end
function loadUnitcell(filename::String)
    # return the filename
    return Unitcell(filename)
end





# make the type public accessible
export Unitcell
export loadUnitcell, saveUnitcell



















#-----------------------------------------------------------------------------------------------------------------------------
#
#   THE LATTICE OBJECT CLASS IN JULIA
#
#   FIELDS:
#   - unitcell              The unit cell of the lattice (if present)
#   - unitcellRepetitions   Number of repetitions for each lattice direction (integer array with length of lattice vectors of UC)
#   - lattice_vectors       all periodic lattice vectors that span the complete lattice (empty list if finite lattice)
#   - positions             all sites of the lattice
#   - connections           the connetions between sites
#   - filename              The filename that the lattice should be stored into
#
#-----------------------------------------------------------------------------------------------------------------------------
"""
    mutable struct Lattice

The type that contains information on a lattice. Fields are

	unitcell            :: Unitcell
	unitcellRepetitions :: Array{Int64, 1}
    lattice_vectors     :: Array{Array{Float64, 1}, 1}
    positions           :: Array{Array{Float64, 1}, 1}
	positions_indices   :: Array{Int64, 1}
    connections         :: Array{Array{Any, 1}, 1}
    filename            :: String

Note that connections have the format

    [index_from::Int, index_to::Int, strength::Any, lattice_offset::Tuple]

where `lattice_offset` only describes the offset along the periodic directions

New objects can be created with

    Lattice(filename::String)

to load from a given file with name `filename` or with

    Lattice()

if the constructed object serves as a placeholder. In this case, the filename
is set to `LATTICE_DUMMY_FILENAME`.

For a complete list of all constructors please consult the documentation.



# Examples

```julia-repl
julia> Lattice("mylatticefile.jld")
```
"""
type Lattice

    # FIELDS

	# The unitcell of the lattice if applicable
    unitcell::Unitcell
	# how often the unitcell is repeated in the elementary directions
    unitcellRepetitions::Array{Int64, 1}

	# The lattice vectors that span the LATTICE in periodic dimensions
    lattice_vectors::Array{Array{Float64, 1}, 1}

	# Array with the positions of all sites within the lattice
    positions::Array{Array{Float64, 1}, 1}
	# array of ints that give the index of the site in the original unit cell
    positions_indices::Array{Int64, 1}

	# Array with connections between sites (same format as in Unitcell)
    connections::Array{Array{Any, 1}, 1}

	# Filename of the jld file
    filename::String


	# general Constructor if everything is given
    function Lattice(
				unitcell::Unitcell,
				unitcellRepetitions::Array{Int64, 1},
				lattice_vectors::Array{Array{Float64, 1}, 1},
				positions::Array{Array{Float64, 1}, 1},
				positions_indices::Array{Int64, 1},
				connections::Array{Array{Any, 1}, 1},
				filename::String
			)
        # just initialize everything and set the position indices to 1
        return new(unitcell, unitcellRepetitions, lattice_vectors, positions, positions_indices, connections, filename)
    end
    function Lattice(
				unitcell::Unitcell,
				unitcellRepetitions::Array{Int64, 1},
				lattice_vectors::Array{Array, 1},
				positions::Array{Array, 1},
				positions_indices::Array{Int64, 1},
				connections::Array{Array, 1},
				filename::String
			)
        # just initialize everything and set the position indices to 1
        return new(unitcell, unitcellRepetitions, lattice_vectors, positions, positions_indices, connections, filename)
    end



	# Constructor if no position indices should be given
    function Lattice(
				unitcell::Unitcell,
				unitcellRepetitions::Array{Int64, 1},
				lattice_vectors::Array{Array{Float64, 1}, 1},
				positions::Array{Array{Float64, 1}, 1},
				connections::Array{Array{Any, 1}, 1},
				filename::String
			)
        # just initialize everything and set the position indices to 1
        return new(unitcell, unitcellRepetitions, lattice_vectors, positions,ones(size(positions,1)), connections, filename)
    end

    # Constructor if the lattice is not constructed from a unitcell
    function Lattice(
				lattice_vectors::Array{Array{Float64, 1}, 1},
				positions::Array{Array{Float64, 1}, 1},
				connections::Array{Array{Any, 1}, 1},
				filename::String
			)
        # just initialize everything
        return new(Unitcell(), [], lattice_vectors, positions,ones(size(positions,1)), connections, filename)
    end

    # Constructor when loading from a file
    function Lattice(filename::String)
        # define an empty unitcell
        lattice = new(Unitcell(),[], nothing,nothing,[],nothing,filename)
        # load the uc
        loadLattice(lattice)
        # return the uc
        return lattice
    end

    # Dummy constructor
    function Lattice()
        # define an empty unitcell
        lattice = new(Unitcell(),[], nothing,nothing,[],nothing,LATTICE_DUMMY_FILENAME)
        # return the uc
        return lattice
    end

end



#-----------------------------------------
#
#   METHODS FOR SAVING LATTICES
#
#-----------------------------------------

"""
    saveLattice(lattice::Lattice [, filename::String])

Save the `Lattice` object `lattice` to a file by using the module `JLD`.
A new filename can be given by passing an additional argument.

Note that the filename has to end with `.jld` or the ending will be added automatically.

Folders that lead to the specified file will be created on the way.



# Examples

```julia-repl
julia> saveLattice(lattice)
"filename.jld"

julia> saveLattice(lattice, "testfile.jld")
"testfile.jld"
```
"""
function saveLattice(lattice::Lattice, filename::String="NONE")
    # maybe set the new filename
    if filename != "NONE"
        lattice.filename = filename
    end
    # make sure that the file ends with .jld
    if lattice.filename[end-4:end] != ".jld"
        lattice.filename = "$(filename).jld"
    end
    # ensure default path is build
    if contains(lattice.filename, "/")
		# get the containing folder
		folder = lattice.filename[1:findlast(lattice.filename, '/')]
		# build the path to that folder
		mkpath(folder)
	end
    # open the file and write the fields into it
    save(lattice.filename,
        "unitcell",             lattice.unitcell,
        "unitcell repetitions", lattice.unitcellRepetitions,
        "positions",            lattice.positions,
        "positions indices",    lattice.positions_indices,
        "connections",          lattice.connections,
        "lattice vectors",      lattice.lattice_vectors,
        "filename",             lattice.filename
    )
    # return the filename
    return lattice.filename
end






#-----------------------------------------
#
#   METHODS FOR LOADING LATTICES
#
#-----------------------------------------


"""
    loadLattice(lattice::Lattice [, filename::String])
    loadLattice(filename::String)

Loads a `Lattice` object `lattice` from its file by using the module `JLD`.
A new filename can be given by passing an additional argument.
Moreover, if there is no `Lattice` object which should be loaded but one wants to obtain a new
object, one can only pass the filename of the unitcell file. (Note, that this can also be done
by using the constructor `Lattice(filename)`)



# Examples

```julia-repl
julia> loadLattice(lattice)
"filename.jld"

julia> loadLattice(lattice, "testfile.jld")
"testfile.jld"

julia> loadLattice("testfile.jld")
LatticePhysics.Lattice[...]
```
"""
function loadLattice(lattice::Lattice, filename::String="NONE")
    # maybe set the new filename
    if filename != "NONE"
        lattice.filename = filename
    end
    # open the file and load all fields from it
    lattice.unitcell            = load(lattice.filename, "unitcell")
    lattice.unitcellRepetitions = load(lattice.filename, "unitcell repetitions")
    lattice.positions           = load(lattice.filename, "positions")
    lattice.positions_indices   = load(lattice.filename, "positions indices")
    lattice.connections         = load(lattice.filename, "connections")
    lattice.lattice_vectors     = load(lattice.filename, "lattice vectors")
    # return the filename
    return lattice.filename
end
function loadLattice(filename::String)
    # just return the object that the constructor would return
	return Lattice(filename)
end




# EXPORT the relevant types and methods
export Lattice
export saveLattice, loadLattice
























#-----------------------------------------------------------------------------
#
#   METHODS FOR CONVERTING LATTICE --> UNITCELL
#
#-----------------------------------------------------------------------------

"""
	toUnitcell(lattice::Lattice [, filename::String])

Converts a `Lattice` object `lattice` to a `Unitcell` object by keeping
connections and lattice vectors and setting all sites as basis sites.
Optionally, a further argument can specify a new filename for the newly
created unitcell.


# Examples

```julia-repl
julia> toUnitcell(lattice)
LatticePhysics.Unitcell[...]

julia> toUnitcell(lattice, "unitcell_test.jld")
LatticePhysics.Unitcell[...]
```
"""
function toUnitcell(lattice::Lattice, filename::String="AUTOMATIC")
    # maybe set the new filename
    if filename == "AUTOMATIC"
        filename = replace(lattice.filename, ".jld", "_to_unitcell.jld")
    end
	# return a newly built unitcell object
    return Unitcell(
        lattice.lattice_vectors,
        lattice.positions,
        lattice.connections,
		filename
            )
end

# Export the conversion method
export toUnitcell




















#-----------------------------------------------------------------------------
#
#   METHODS FOR PRINTING OR OBTAINING GENERAL INFORMATION
#
#-----------------------------------------------------------------------------


"""
    printInfo(unitcell::Unitcell [; detailed::Bool=false])
    printInfo(lattice::Lattice [; detailed::Bool=false])

Prints (detailed) information about a `Unitcell` object or a `Lattice` object.
Information Includes Bravais lattice vectors, number of sites and number of bonds.
When detailed output is enabled, furthermore all sites and bonds are printed.
Last but now least, statistics is printed about how many bonds per site are found.


# Examples

```julia-repl
julia> printInfo(unitcell)
...

julia> printInfo(unitcell, detailed=true)
...

```
"""
function printInfo(unitcell::Unitcell; detailed::Bool=false)

    # Header of the information
    println("Information on the unitcell stored in file \"$(unitcell.filename)\":")

    # General information on Bravais lattice vectors
    println(" - periodicity given by $(size(unitcell.lattice_vectors,1)) lattice vectors:")
    for l in unitcell.lattice_vectors
        println("     - $(l)")
    end

    # Information on basis sites
    if detailed
        println(" - $(size(unitcell.basis,1)) sites in unitcell of dimension $(length(unitcell.basis[1])):")
        for site in unitcell.basis
            println("     - $(site)")
        end
    else
        println(" - $(size(unitcell.basis,1)) sites in unitcell of dimension $(length(unitcell.basis[1]))")
    end

    # Information on bonds / connections
    if detailed
        println(" - $(size(unitcell.connections,1)) connections in the unitcell:")
        for c in unitcell.connections
            if typeof(c[3]) == String
                println("     - from $(c[1]) to $(c[2]) (with warp $(c[4])): \"$(c[3])\"")
            else
                println("     - from $(c[1]) to $(c[2]) (with warp $(c[4])): $(c[3])")
            end
        end
    else
        println(" - $(size(unitcell.connections,1)) connections in the unitcell")
    end

    # Information on bond / connection statistics
    println(" - $(size(unitcell.connections,1)/size(unitcell.basis,1)) connections per site")
    # check statistics of connections
    nc = zeros(Int64, size(unitcell.basis,1))
    for i in 1:length(nc)
        for c in unitcell.connections
            if Int(c[1]) == i
                nc[i] += 1
            end
        end
    end
    tc = 0
    c = 0
    print(" - statistics of connections per site:")
    while tc < size(unitcell.basis, 1)
        cc = 0
        for ncc in nc
            if ncc == c
                cc += 1
            end
        end
        if cc != 0
            print(" $(c)($(cc))")
        end
        c = c+1
        tc += cc
    end
    println("")
    broken = false
    for c1 in unitcell.connections
        counterpart = false
        for c2 in unitcell.connections
            if c1==c2
                continue
            end
            if c1[1] != c2[2] || c1[2] != c2[1]
                continue # site indices not correct
            end
            if c1[3] != c2[3]
                continue # connection strength not equal
            end
            if sum([abs(element) for element in collect(c1[4]) .+ collect(c2[4])]) > 1e-9
                continue # wrap not equal
            end
            counterpart = true
            break
        end
        if counterpart == false
            broken = true
        end
    end
    if broken
        print(" - connectivity of unitcell is broken (connections not vice versa)")
    else
        print(" - connectivity of unitcell is okay (including periodic BC)")
    end
    println("")
end
function printInfo(lattice::Lattice; detailed=false)
	# Print the header informations
    println("Information on the lattice stored in file \"$(lattice.filename)\":")

	# Information on lattice vectors, print all of them
    println(" - periodicity given by $(size(lattice.lattice_vectors,1)) lattice vectors:")
    for l in lattice.lattice_vectors
        println("     - $(l)")
    end

	# Information on sites, depending on detailed print information
    if detailed
        println(" - $(size(lattice.positions,1)) sites in lattice of dimension $(length(lattice.positions[1])):")
        for site in lattice.positions
            println("     - $(site)")
        end
    else
        println(" - $(size(lattice.positions,1)) sites in lattice of dimension $(length(lattice.positions[1]))")
    end

	# Information on connectioins, depending on detailed print information
    if detailed
        println(" - $(size(lattice.connections,1)) connections in the lattice:")
        for c in lattice.connections
            if typeof(c[3]) == String
                println("     - from $(c[1]) to $(c[2]) (with warp $(c[4])): \"$(c[3])\"")
            else
                println("     - from $(c[1]) to $(c[2]) (with warp $(c[4])): $(c[3])")
            end
        end
    else
        println(" - $(size(lattice.connections,1)) connections in the lattice")
    end
	# Number of connections per site
    println(" - $(size(lattice.connections,1)/size(lattice.positions,1)) connections per site")
    # check statistics of connections
    nc = zeros(Int64, size(lattice.positions,1))
    for i in 1:length(nc)
        for c in lattice.connections
            if Int(c[1]) == i
                nc[i] += 1
            end
        end
    end
    tc = 0
    c = 0
	# print the statistics of the connections
    print(" - statistics of connections per site:")
    while tc < size(lattice.positions, 1)
        cc = 0
        for ncc in nc
            if ncc == c
                cc += 1
            end
        end
        if cc != 0
            print(" $(c)($(cc))")
        end
        c = c+1
        tc += cc
    end
    println("")
	# Test the connectivity of the lattice
    broken = false
    for c1 in lattice.connections
        counterpart = false
        for c2 in lattice.connections
            if c1==c2
                continue
            end
            if c1[1] != c2[2] || c1[2] != c2[1]
                continue # site indices not correct
            end
            if c1[3] != c2[3]
                continue # connection strength not equal
            end
            if sum([abs(element) for element in collect(c1[4]) .+ collect(c2[4])]) > 1e-9
                continue # wrap not equal
            end
            counterpart = true
            break
        end
        if counterpart == false
            broken = true
        end
    end
	# Print the connectivity of the lattice
    if broken
        println(" - connectivity of lattice is broken (connections not vice versa)")
    else
        println(" - connectivity of lattice is okay (including periodic BC)")
    end
end
export printInfo





"""
    getConnectionStrengthList(unitcell::Unitcell)
	getConnectionStrengthList(lattice::Lattice)

Obtains which connection strength values are taken in the definition of a given
`Lattice` object `lattice` or a `Unitcell` object `unitcell`.
The functions return a list which contains all connection strengths in the order
of appearence.


# Examples

```julia-repl
julia> connection_strength_list = getConnectionStrengthList(lattice)
["tx", "ty", "tz"]
```
"""
function getConnectionStrengthList(unitcell::Unitcell)
    # list of connetion strengths
    cs_list = []
    # iterate over all connections
    for c in unitcell.connections
        if !(c[3] in cs_list)
            push!(cs_list, c[3])
        end
    end
    # return the list
    return cs_list
end
function getConnectionStrengthList(lattice::Lattice)
    # list of connetion strengths
    cs_list = []
    # iterate over all connections
    for c in lattice.connections
        if !(c[3] in cs_list)
            push!(cs_list, c[3])
        end
    end
    # return the list
    return cs_list
end

# export the FUNCTIONS
export getConnectionStrengthList
























#-----------------------------------------------------------------------------
#
#   METHODS FOR CONNECTIVITY INFORMATION
#
#-----------------------------------------------------------------------------


"""
	getConnectivityList(lattice::Lattice)
    getConnectivityList(unitcell::Unitcell)

Obtains connectivity information of a `Lattice` object `lattice` or a `Unitcell` object `unitcell`.
The functions return a list which contains a list of Tuples for every site in the
passed objects. Therefore, `connectivity_list[i][j]` can be used to access the connection `j`
from site `i` which will give a Tuple `(to_index, strength)`.


# Examples

```julia-repl
julia> connectivity_list = getConnectivityList(lattice)
Array{...}[...]

julia> connectivity_list[2][1]   # connection 1 outgoing from site 2
(12, "tx")
```
"""
function getConnectivityList(lattice::Lattice)
    # list of lists
    lists = Array[]
    for p in lattice.positions
        push!(lists, [])
    end
    for c in lattice.connections
        push!(lists[Int(c[1])], (c[2], c[3]) )
    end
    # return the lists
    return lists
end
function getConnectivityList(unitcell::Unitcell)
    # list of lists
    lists = Array[]
    for p in unitcell.basis
        push!(lists, [])
    end
    for c in unitcell.connections
        push!(lists[Int(c[1])], (c[2], c[3]) )
    end
    # return the lists
    return lists
end

# export the FUNCTIONS
export getConnectivityList





"""
	getConnectionList(lattice::Lattice)
    getConnectionList(unitcell::Unitcell)

Obtains connectivity information of a `Lattice` object `lattice` or a `Unitcell` object `unitcell`.
The functions return a list which contains a list of connection Arrays for every site in the
passed objects. Therefore, `connectivity_list[i][j]` can be used to access the connection `j`
from site `i` which will give an Array `[from_index, to_index, strength, wrap_tuple]`
like it is saved in the list of all connections of the respective object.


# Examples

```julia-repl
julia> connection_list = getConnectionList(lattice)
Array{...}[...]

julia> connection_list[2][1]   # connection 1 outgoing from site 2
[1, 12, "tx", (0,0)]
```
"""
function getConnectionList(lattice::Lattice)
    # list of lists
    lists = Array[]
    for p in lattice.positions
        push!(lists, Array[])
    end
    for c in lattice.connections
        push!(lists[Int(c[1])], c )
    end
    # return the lists
    return lists
end
function getConnectionList(unitcell::Unitcell)
    # list of lists
    lists = Array[]
    for p in unitcell.basis
        push!(lists, Array[])
    end
    for c in unitcell.connections
        push!(lists[Int(c[1])], c )
    end
    # return the lists
    return lists
end

# export the FUNCTIONS
export getConnectionList
