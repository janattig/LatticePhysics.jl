#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#
#   TYPE DEFINITIONS OF OBJECT CLASSES IN JULIA
#
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------






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

    lattice_vectors::Array{Array{Float64, 1}, 1}
    basis::Array{Array{Float64, 1}, 1}
    connections::Array{Array{Any, 1}, 1}
    filename::String

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
    mkpath(dirname(uc.filename))
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
Lattice.Unitcell[...]
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




#-----------------------------------------
#
#   METHODS FOR PRINTING INFORMATION
#
#-----------------------------------------


"""
    printInfo(uc::Unitcell [; detailed::Bool=false])

Prints (detailed) information about a `Unitcell` object. Information Includes
Bravais lattice vectors, number of sites and number of bonds. When detailed output is
enabled, furthermore all sites and bonds are printed. Last but now least, statistics
is printed about how many bonds per site are found.


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
type Lattice

    # FIELDS
    unitcell::Unitcell
    unitcellRepetitions::Array{Int64, 1}
    lattice_vectors     # that span the lattice in periodic dimensions
    positions           # positions of all sites
    positions_indices   # array of ints that give the index of the site in the original unit cell
    connections         # connections between all sites
    filename

    # the overall usual constructor
    function Lattice(unitcell::Unitcell, unitcellRepetitions, lattice_vectors, positions, positions_indices, connections, filename)
        # just initialize everything
        return new(unitcell, unitcellRepetitions, lattice_vectors, positions, positions_indices, connections, filename)
    end
    function Lattice(unitcell::Unitcell, unitcellRepetitions, lattice_vectors, positions, connections, filename)
        # just initialize everything
        return new(unitcell, unitcellRepetitions, lattice_vectors, positions,ones(size(positions,1)), connections, filename)
    end
    # the custom constructor
    function Lattice(lattice_vectors, positions, connections, filename)
        # just initialize everything
        return new(Unitcell(), [], lattice_vectors, positions,ones(size(positions,1)), connections, filename)
    end
    # the constructor when loading from a file
    function Lattice(filename)
        # define an empty unitcell
        lattice = new(Unitcell(),[], nothing,nothing,[],nothing,filename)
        # load the uc
        loadLattice(lattice)
        # return the uc
        return lattice
    end
    # dummy constructor
    function Lattice()
        # define an empty unitcell
        lattice = new(Unitcell(),[], nothing,nothing,[],nothing,LATTICE_DUMMY_FILENAME)
        # return the uc
        return lattice
    end

end



# METHODS FOR SAVING

# save a lattice to the lattice file specified
function saveLattice(lattice::Lattice)
    # ensure default path is build
    buildFolderLattices()
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
# save a lattice to a new lattice file specified
function saveLattice(lattice::Lattice, filename)
    # set the new filename
    lattice.filename = filename
    # ensure default path is build
    buildFolderLattices()
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



# METHODS FOR LOADING

# method for loading the lattice
function loadLattice(lattice::Lattice)
    # ensure default path is build
    buildFolderLattices()
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
# method for loading the unit cell from a different filename
function loadLattice(lattice::Lattice, filename_new)
    # overwrite the existing filename
    lattice.filename            = filename_new
    # ensure default path is build
    buildFolderLattices()
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




# METHOD FOR PRINTING INFORMATION ABOUT THE LATTICE
function printInfo(lattice::Lattice; detailed=false)
    println("Information on the lattice stored in file \"$(lattice.filename)\":")
    println(" - periodicity given by $(size(lattice.lattice_vectors,1)) lattice vectors:")
    for l in lattice.lattice_vectors
        println("     - $(l)")
    end
    if detailed
        println(" - $(size(lattice.positions,1)) sites in lattice of dimension $(length(lattice.positions[1])):")
        for site in lattice.positions
            println("     - $(site)")
        end
    else
        println(" - $(size(lattice.positions,1)) sites in lattice of dimension $(length(lattice.positions[1]))")
    end
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
    if broken
        print(" - connectivity of lattice is broken (connections not vice versa)")
    else
        print(" - connectivity of lattice is okay (including periodic BC)")
    end
    println("")
end




# METHODS FOR OBTAINING STRUCTURED CONNECTION INFORMATION

# get for every site a list with connected sites and connection strengths
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
# get for every site a list with all connections
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

# get a list of all connection strengths
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


# METHOD FOR DEFINING A LATTICE AS A UNITCELL
function toUnitcell(lattice::Lattice)
    return Unitcell(
        lattice.lattice_vectors,
        lattice.positions,
        lattice.connections,
        replace(lattice.filename, ".jld", "_to_unitcell.jld")
            )
end



# EXPORT the relevant types and methods
export Lattice
export saveLattice, loadLattice
export printInfo
export getConnectivityList, getConnectionList, getConnectionStrengthList
export toUnitcell
