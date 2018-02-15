#-----------------------------------------------------------------------------------------------------------------------------
#
#   The unified Module for all different sorts of lattices
#   Includes the following structure:
#
#   -   Definition of UNITCELL TYPE
#       -   saving / loading Unitcell
#   -   Definition of LATTICE TYPE
#       -   saving / loading Lattice
#       -   information on Lattice
#       -   connectivities of Lattice
#
#   -   Functions to generate Unitcell objects for 2D lattices
#   -   Functions to generate Unitcell objects for 3D lattices
#   -   TODO Function to generate Unitcell object from collection of points
#   -   Function to print generating code for a given Unitcell object
#
#   -   Building Lattices as (periodic/open/semiperiodic) repeated patterns of Unitcells for 2D unitcells
#   -   Building Lattices as (periodic/open/semiperiodic) repeated patterns of Unitcells for 3D unitcells
#   -   Building Lattices as (periodic/open/semiperiodic) repeated patterns of Unitcells for any unitcell
#   -   Building Lattices by bond distance to an origin site
#   -   Building Lattices in a shape around an origin site
#
#   -   Modifying Lattices and Interaction Strengths
#
#   -   TODO SVG help methods
#   -   Plotting Lattices as SVG files
#   -   TODO Calculating dispersions for 1D and 2D fields of k values
#   -   TODO Calculation of Fermi surfaces or surfaces of some energy value
#
#
#   AUTHOR:             Jan Attig
#   DATE started:       2017-08-16
#   DATE last version:  2017-09-12
#
#-----------------------------------------------------------------------------------------------------------------------------


# Start of module
module LatticePhysics








#-----------------------------------------------------------------------------------------------------------------------------
#
#   DEPENDENCIES AND USED LIBRARIES
#
#-----------------------------------------------------------------------------------------------------------------------------

# JLD is used for all save / load applications
using JLD

# PyPlot is used for all plotting of band structures (plotting of lattices is done via SVG file creation)
using PyPlot
using PyCall

# Optim is used for minimizing the band structure to find the ground state energy of the system
using Optim







#-----------------------------------------------------------------------------------------------------------------------------
#
#   SOME FILENAMES AND DEFAULT PATHS
#
#-----------------------------------------------------------------------------------------------------------------------------

# THE FOLDER FOR UNITCELLS
FOLDER_UNITCELLS = "unitcells/"
# THE FOLDER FOR LATICES
FOLDER_LATTICES = "lattices/"

# export the folders
export FOLDER_UNITCELLS
export FOLDER_LATTICES

# Dummy filenames
UNITCELL_DUMMY_FILENAME = "THIS IS NOT A UNIT CELL"
LATTICE_DUMMY_FILENAME = "THIS IS NOT A LATTICE"
export UNITCELL_DUMMY_FILENAME
export LATTICE_DUMMY_FILENAME




# THE FOLDER FOR BAND STRUCTURES
FOLDER_SPECTRA = "band_structures/"

# export the folder
export FOLDER_SPECTRA






#-----------------------------------------------------------------------------------------------------------------------------
#
#   FUNCTIONS TO ENSURE DEFAULT PATHS EXIST
#
#-----------------------------------------------------------------------------------------------------------------------------
function buildFolderUnitcells()
    if !isdir(FOLDER_UNITCELLS)
        mkdir(FOLDER_UNITCELLS)
    end
end
export buildFolderUnitcells
function buildFolderLattices()
    if !isdir(FOLDER_LATTICES)
        mkdir(FOLDER_LATTICES)
    end
end
export buildFolderLattices
function buildFolderSpectra()
    if !isdir(FOLDER_SPECTRA)
        mkdir(FOLDER_SPECTRA)
    end
end
export buildFolderSpectra


function buildFolders()
    buildFolderUnitcells()
    buildFolderLattices()
    buildFolderSpectra()
end
export buildFolders



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
type Unitcell

    # FIELDS
    lattice_vectors
    basis
    connections
    filename

    # the custom constructor
    function Unitcell(lattice_vectors, basis, connections, filename)
        # just initialize everything
        return new(lattice_vectors, basis, connections, filename)
    end
    # the constructor when loading from a file
    function Unitcell(filename)
        # define an empty unitcell
        uc = new(nothing,nothing,nothing,filename)
        # load the uc
        loadUnitcell(uc)
        # return the uc
        return uc
    end
    # the dummy constructor for lattices that do not need a UC
    function Unitcell()
        # just initialize everything
        return new([], [], [], UNITCELL_DUMMY_FILENAME)
    end
end


# methods for saving the unit cell
function saveUnitcell(uc::Unitcell)
    # ensure default path is build
    buildFolderUnitcells()
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
function saveUnitcell(uc::Unitcell, filename_new)
    # set the new filename
    uc.filename = filename_new
    # ensure default path is build
    buildFolderUnitcells()
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

# methods for loading the unit cell
function loadUnitcell(uc::Unitcell)
    # ensure default path is build
    buildFolderUnitcells()
    # open the file and load all fields from it
    uc.basis            = load(uc.filename, "basis")
    uc.connections      = load(uc.filename, "connections")
    uc.lattice_vectors  = load(uc.filename, "lattice vectors")
    # return the filename
    return uc.filename
end
function loadUnitcell(uc::Unitcell, filename_new)
    # overwrite the existing filename
    uc.filename         = filename_new
    # ensure default path is build
    buildFolderUnitcells()
    # open the file and load all fields from it
    uc.basis            = load(uc.filename, "basis")
    uc.connections      = load(uc.filename, "connections")
    uc.lattice_vectors  = load(uc.filename, "lattice vectors")
    # return the filename
    return uc.filename
end

# method for loading the unit cell and returning it
function loadUnitcell(filename)
    # return the filename
    return Unitcell(filename)
end


# print information about unitcell (export with lattice printInfo)
function printInfo(unitcell::Unitcell; detailed=false)
    println("Information on the unitcell stored in file \"$(unitcell.filename)\":")
    println(" - periodicity given by $(size(unitcell.lattice_vectors,1)) lattice vectors:")
    for l in unitcell.lattice_vectors
        println("     - $(l)")
    end
    if detailed
        println(" - $(size(unitcell.basis,1)) sites in unitcell of dimension $(length(unitcell.basis[1])):")
        for site in unitcell.basis
            println("     - $(site)")
        end
    else
        println(" - $(size(unitcell.basis,1)) sites in unitcell of dimension $(length(unitcell.basis[1]))")
    end    
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




















#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#
#   UNITCELL DEFINITIONS
#
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------------------------------
#
#   Individual FUNCTIONS FOR 2D UNITCELLS
#
#   They all generate the special filenames for the unitcells and return the unitcells as objects
#   Functions can be given a version integer to distinguish several implementations of the same lattice
#   Functions can be specified to already save the unitcell to file
#
#-----------------------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------------------
# SQUARE LATTICE
# 1 - simple, 1 site per UC
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellSquare(version=1; save=true)
    # SIMPLE SQUARE LATTICE
    if version == 1
        # the lattice vectors
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 1; 1.0; (0, 1)],
            [1; 1; 1.0; (0, -1)],
            [1; 1; 1.0; (1, 0)],
            [1; 1; 1.0; (-1, 0)],
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_square_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellSquare

#-----------------------------------------------------------------------------------------------------------------------------
# EXTENDED SQUARE LATTICE
# 1 - simple, 3 sites per UC
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellExtendedSquare(version=1; save=true)
    # SIMPLE EXTENDED SQUARE LATTICE
    if version == 1
        # the lattice vectors
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [0.5, 0.0],
            [0.0, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; 1.0; (0, 0)],
            [1; 2; 1.0; (-1, 0)],
            [1; 3; 1.0; (0, 0)],
            [1; 3; 1.0; (0, -1)],

            [2; 1; 1.0; (0, 0)],
            [2; 1; 1.0; (1, 0)],

            [3; 1; 1.0; (0, 0)],
            [3; 1; 1.0; (0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_extended_square_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellExtendedSquare

#-----------------------------------------------------------------------------------------------------------------------------
# CHECKERBOARD LATTICE
# 1 - simple, 2 sites per UC
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellCheckerboard(version=1; save=true)
    # SIMPLE CHECKERBOARD LATTICE
    if version == 1
        # the lattice vectors
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [0.5, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 1; 1.0; (1, 0)],
            [1; 1; 1.0; (-1, 0)],
            [1; 2; 1.0; (0, 0)],
            [1; 2; 1.0; (0, -1)],
            [1; 2; 1.0; (-1, 0)],
            [1; 2; 1.0; (-1, -1)],

            [2; 2; 1.0; (0, 1)],
            [2; 2; 1.0; (0, -1)],
            [2; 1; 1.0; (0, 0)],
            [2; 1; 1.0; (0, 1)],
            [2; 1; 1.0; (1, 0)],
            [2; 1; 1.0; (1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_checkerboard_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellCheckerboard


#-----------------------------------------------------------------------------------------------------------------------------
# SHASTRY SUTHERLAND LATTICE
# 1 - simple, 4 sites per UC
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellShastrySutherland(version=1; save=true)
    # SIMPLE SHASTRY SUTHERLAND LATTICE
    if version == 1
        # the lattice vectors
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [0.5, 0.0],
            [0.5, 0.5],
            [0.0, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 3; "t1"; (0, 0)],
            [1; 2; "t2"; (-1,0)],
            [1; 4; "t3"; (0,-1)],
            [1; 2; "t4"; (0, 0)],
            [1; 4; "t5"; (0, 0)],

            [2; 4; "t1"; (1,-1)],
            [2; 1; "t2"; (1, 0)],
            [2; 3; "t3"; (0, 0)],
            [2; 1; "t4"; (0, 0)],
            [2; 3; "t5"; (0,-1)],

            [3; 1; "t1"; (0, 0)],
            [3; 4; "t2"; (0, 0)],
            [3; 2; "t3"; (0, 0)],
            [3; 4; "t4"; (1, 0)],
            [3; 2; "t5"; (0, 1)],

            [4; 2; "t1"; (-1,1)],
            [4; 3; "t2"; (0, 0)],
            [4; 1; "t3"; (0, 1)],
            [4; 3; "t4"; (-1,0)],
            [4; 1; "t5"; (0, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_shastry_sutherland_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellShastrySutherland


#-----------------------------------------------------------------------------------------------------------------------------
# ADVANCED SQUARE LATTICE
# 1 - simple, 4 sites per UC, 2 different couplings
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellAdvancedSquare(version=1; save=true)
    # SIMPLE ADVANCED SQUARE LATTICE
    if version == 1
        # the lattice vectors
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [0.5, 0.0],
            [0.5, 0.5],
            [0.0, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; "t2"; (-1,0)],
            [1; 4; "t2"; (0,-1)],
            [1; 2; "t1"; (0, 0)],
            [1; 4; "t1"; (0, 0)],

            [2; 1; "t2"; (1, 0)],
            [2; 3; "t1"; (0, 0)],
            [2; 1; "t1"; (0, 0)],
            [2; 3; "t2"; (0,-1)],

            [3; 4; "t1"; (0, 0)],
            [3; 2; "t1"; (0, 0)],
            [3; 4; "t2"; (1, 0)],
            [3; 2; "t2"; (0, 1)],

            [4; 3; "t1"; (0, 0)],
            [4; 1; "t2"; (0, 1)],
            [4; 3; "t2"; (-1,0)],
            [4; 1; "t1"; (0, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_advanced_square_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellAdvancedSquare
#-----------------------------------------------------------------------------------------------------------------------------
# SQUARE OCTAGON LATTICE
# 1 - simple, 4 sites per UC
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellSquareOctagon(version=1; save=true)
    # SQUARE OCTAGON LATTICE
    if version == 1
        # the lattice vectors
        a1 = [3*sqrt(3.0)/4., -3*sqrt(3.0)/4.]
        a2 = [3*sqrt(3.0)/4., 3*sqrt(3.0)/4.]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [0.0, 1/sqrt(3.0)],
			[-1/sqrt(3.0), 0.0],
			[-1/sqrt(3.0), 1/sqrt(3.0)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; 1.0; (0, 0)],
			[2; 4; 1.0; (0, 0)],
			[3; 4; 1.0; (0, 0)],
			[3; 1; 1.0; (0, 0)],
			[2; 1; 1.0; (0, 0)],
			[4; 2; 1.0; (0, 0)],
			[4; 3; 1.0; (0, 0)],
			[1; 3; 1.0; (0, 0)],

            [3; 2; 1.0; (0, -1)],
            [2; 3; 1.0; (0, 1)],
            [4; 1; 1.0; (-1, 0)],
            [1; 4; 1.0; (1, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_squareoctagon_unitcell.jld"
    elseif version == 2
        # the lattice vectors
        a1 = [1.0 + 1./sqrt(2.0), -(1.0 + 1./sqrt(2.0))]
        a2 = [1.0 + 1./sqrt(2.0),  (1.0 + 1./sqrt(2.0))]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [ 0.5, -0.5],
            [ 0.5,  0.5],
            [-0.5, -0.5],
            [-0.5,  0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; 1.0; (0, 0)],
            [2; 4; 1.0; (0, 0)],
            [3; 4; 1.0; (0, 0)],
            [3; 1; 1.0; (0, 0)],
            [2; 1; 1.0; (0, 0)],
            [4; 2; 1.0; (0, 0)],
            [4; 3; 1.0; (0, 0)],
            [1; 3; 1.0; (0, 0)],
            [3; 2; 1.0; (0, -1)],
            [2; 3; 1.0; (0, 1)],
            [4; 1; 1.0; (-1, 0)],
            [1; 4; 1.0; (1, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_squareoctagon_alternative_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [3*sqrt(3.0)/4., -3*sqrt(3.0)/4.]
        a2 = [3*sqrt(3.0)/4., 3*sqrt(3.0)/4.]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [0.0, 1/sqrt(3.0)],
			[-1/sqrt(3.0), 0.0],
			[-1/sqrt(3.0), 1/sqrt(3.0)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[          
            [1; 2; "tx"; (0,0)],
            [2; 4; "ty"; (0,0)],
            [3; 4; "tx"; (0,0)],
            [3; 1; "ty"; (0,0)],
            [2; 1; "tx"; (0,0)],
            [4; 2; "ty"; (0,0)],
            [4; 3; "tx"; (0,0)],
            [1; 3; "ty"; (0,0)],
            [3; 2; "tz"; (0,-1)],
            [2; 3; "tz"; (0,1)],
            [4; 1; "tz"; (-1,0)],
            [1; 4; "tz"; (1,0)],
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_squareoctagon_kitaev_unitcell.jld"
    elseif version == 5
        # the lattice vectors
        a1 = [1.0 + 1./sqrt(2.0), -(1.0 + 1./sqrt(2.0))]
        a2 = [1.0 + 1./sqrt(2.0),  (1.0 + 1./sqrt(2.0))]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [ 0.5, -0.5],
            [ 0.5,  0.5],
            [-0.5, -0.5],
            [-0.5,  0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[          
            [1; 2; "tx"; (0,0)],
            [2; 4; "ty"; (0,0)],
            [3; 4; "tx"; (0,0)],
            [3; 1; "ty"; (0,0)],
            [2; 1; "tx"; (0,0)],
            [4; 2; "ty"; (0,0)],
            [4; 3; "tx"; (0,0)],
            [1; 3; "ty"; (0,0)],
            [3; 2; "tz"; (0,-1)],
            [2; 3; "tz"; (0,1)],
            [4; 1; "tz"; (-1,0)],
            [1; 4; "tz"; (1,0)],
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_squareoctagon_alternative_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellSquareOctagon

#-----------------------------------------------------------------------------------------------------------------------------
# BCC LATTICE in 2D
# just another representation of the square lattice
# 1 - simple, 2 sites per UC
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellBCC2D(version=1; save=true)
    if version == 1
        # the lattice vectors
        a1 = [1, 0]
        a2 = [0, 1]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [0.5, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; 1.0; (0, 0)],
            [1; 2; 1.0; (-1, 0)],
            [1; 2; 1.0; (0, -1)],
            [1; 2; 1.0; (-1, -1)],

            [2; 1; 1.0; (0, 0)],
            [2; 1; 1.0; (1, 0)],
            [2; 1; 1.0; (0, 1)],
            [2; 1; 1.0; (1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_unitcell_bcc_2d.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellBCC2D

#-----------------------------------------------------------------------------------------------------------------------------
# FULLY CONNECTED SQUARE LATTICE
# same as square lattice, but with additional X like couplings on the square plaquettes
# 1 - simple, 1 site per UC (same as square), coupling ratio fixed 2:1 or adjustable
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellFullyConnectedSquare(version=1; save=true, J1=1.0, J1X=0.5)
    if version == 1
        if J1==1.0 && J1X==0.5
            # the lattice vectors
            a1 = [1.0, 0.0]
            a2 = [0.0, 1.0]
            lattice_vectors = Array[]
            push!(lattice_vectors, a1)
            push!(lattice_vectors, a2)
            # Basis Definition
            basis = Array[
                [0.0, 0.0]
            ]
            # Connection Definition
            # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            connections = Array[
                [1; 1; 1.0; (0, 1)],
                [1; 1; 1.0; (0, -1)],
                [1; 1; 1.0; (1, 0)],
                [1; 1; 1.0; (-1, 0)],
                [1; 1; 0.5; (1, 1)],
                [1; 1; 0.5; (-1, -1)],
                [1; 1; 0.5; (1, -1)],
                [1; 1; 0.5; (-1, 1)]
            ]
            # filename
            filename = "$(FOLDER_UNITCELLS)2d_fully_connected_square_unitcell.jld"
        else
            # the lattice vectors
            a1 = [1.0, 0.0]
            a2 = [0.0, 1.0]
            lattice_vectors = Array[]
            push!(lattice_vectors, a1)
            push!(lattice_vectors, a2)
            # Basis Definition
            basis = Array[
                [0.0, 0.0]
            ]
            # Connection Definition
            # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            connections = Array[
                [1; 1; J1; (0, 1)],
                [1; 1; J1; (0, -1)],
                [1; 1; J1; (1, 0)],
                [1; 1; J1; (-1, 0)],
                [1; 1; J1X; (1, 1)],
                [1; 1; J1X; (-1, -1)],
                [1; 1; J1X; (1, -1)],
                [1; 1; J1X; (-1, 1)]
            ]
            # filename
            filename = "$(FOLDER_UNITCELLS)2d_fully_connected_square_$(J1)_$(J1X)_unitcell.jld"
        end
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellFullyConnectedSquare




#-----------------------------------------------------------------------------------------------------------------------------
# TRIANGULAR LATTICE
# 1 - simple, 1 site per UC (symmetric around x axis)
# 3 - anisotropic coupling, 1 site per UC (symmetric around x axis)
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellTriangular(version=1; save=true)
    # SIMPLE TRIANGULAR LATTICE
    if version == 1
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 1; 1.0; (0, 1)],
            [1; 1; 1.0; (0, -1)],
            [1; 1; 1.0; (1, 0)],
            [1; 1; 1.0; (-1, 0)],
            [1; 1; 1.0; (1, -1)],
            [1; 1; 1.0; (-1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_triangular_unitcell.jld"
    elseif version == 3
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 1; "1"; (0, 1)],
            [1; 1; "1"; (0, -1)],
            [1; 1; "1"; (1, 0)],
            [1; 1; "1"; (-1, 0)],
            [1; 1; "2"; (1, -1)],
            [1; 1; "2"; (-1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_triangular_aniso_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellTriangular

#-----------------------------------------------------------------------------------------------------------------------------
# HONEYCOMB LATTICE
# 1 - simple, 2 sites per UC (symmetric around x axis, gives ZZ edge in strip)
# 2 - simple, 2 sites per UC (symmetric around y axis, gives AC edge in strip)
# 3 - anisotropic hopping, 2 sites per UC (symmetric around x axis, gives ZZ edge in strip)
# 4 - Kitaev couplings
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellHoneycomb(version=1; save=true)
    # SIMPLE HONEYCOMB LATTICE
    if version == 1
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; 1.0; (0, 0)],
            [1; 2; 1.0; (-1, 0)],
            [1; 2; 1.0; (0, -1)],
            [2; 1; 1.0; (0, 0)],
            [2; 1; 1.0; (1, 0)],
            [2; 1; 1.0; (0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_honeycomb_ZZ_unitcell.jld"
    elseif version == 2
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0), 0.0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; 1.0; (0, 0)],
            [1; 2; 1.0; (-1, 0)],
            [1; 2; 1.0; (1, -1)],
            [2; 1; 1.0; (0, 0)],
            [2; 1; 1.0; (1, 0)],
            [2; 1; 1.0; (-1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_honeycomb_AC_unitcell.jld"
    elseif version == 3
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; "1"; (0, 0)],
            [1; 2; "2"; (-1, 0)],
            [1; 2; "2"; (0, -1)],
            [2; 1; "1"; (0, 0)],
            [2; 1; "2"; (1, 0)],
            [2; 1; "2"; (0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_honeycomb_aniso_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; "tx"; (0, 0)],
            [1; 2; "ty"; (-1, 0)],
            [1; 2; "tz"; (0, -1)],
            [2; 1; "tx"; (0, 0)],
            [2; 1; "ty"; (1, 0)],
            [2; 1; "tz"; (0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_honeycomb_ZZ_kitaev_unitcell.jld"
    elseif version == 5
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0), 0.0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; "tx"; (0, 0)],
            [1; 2; "ty"; (-1, 0)],
            [1; 2; "tz"; (1, -1)],
            [2; 1; "tx"; (0, 0)],
            [2; 1; "ty"; (1, 0)],
            [2; 1; "tz"; (-1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_honeycomb_AC_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellHoneycomb


#-----------------------------------------------------------------------------------------------------------------------------
# KAGOME LATTICE
# 1 - simple, 3 sites per UC (symmetric around x axis)
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellKagome(version=1; save=true)
    # SIMPLE KAGOME LATTICE
    if version == 1
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [sqrt(3.0)/4, -0.25],
            [sqrt(3.0)/4, +0.25]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; 1.0; (0, 0)],
            [1; 3; 1.0; (0, 0)],
            [1; 2; 1.0; (-1, 0)],
            [1; 3; 1.0; (0, -1)],

            [2; 1; 1.0; (0, 0)],
            [2; 3; 1.0; (0, 0)],
            [2; 1; 1.0; (1, 0)],
            [2; 3; 1.0; (1, -1)],

            [3; 1; 1.0; (0, 0)],
            [3; 2; 1.0; (0, 0)],
            [3; 1; 1.0; (0, 1)],
            [3; 2; 1.0; (-1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_kagome_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellKagome


#-----------------------------------------------------------------------------------------------------------------------------
# KAGOME MINUS LATTICE
# 1 - simple, 6 sites per UC (symmetric around x axis)
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellKagomeMinus(version=1; save=true)
    # SIMPLE KAGOME MINUS LATTICE
    if version == 1
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        b1 = [0.0, 0.0]
        b2 = [1/sqrt(3.0), 0.0]
        basis = Array[
            b1 .+ 1/4 * (b2 .- b1),
            b1 .+ 1/4 * (b2 .- a2),
            b1 .+ 1/4 * (b2 .- a1),
            b2 .+ 1/4 * (b1 .- b2),
            b2 .+ 1/4 * (a2 .- b2),
            b2 .+ 1/4 * (a1 .- b2)
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; 1.0; (0,0)],
            [2; 1; 1.0; (0,0)],
            [1; 3; 1.0; (0,0)],
            [3; 1; 1.0; (0,0)],
            [2; 3; 1.0; (0,0)],
            [3; 2; 1.0; (0,0)],
            [4; 5; 1.0; (0,0)],
            [5; 4; 1.0; (0,0)],
            [4; 6; 1.0; (0,0)],
            [6; 4; 1.0; (0,0)],
            [5; 6; 1.0; (0,0)],
            [6; 5; 1.0; (0,0)],
            [1; 4; 1.0; (0,0)],
            [4; 1; 1.0; (0,0)],
            [2; 5; 1.0; (0,-1)],
            [5; 2; 1.0; (0,1)],
            [3; 6; 1.0; (-1,0)],
            [6; 3; 1.0; (1,0)],
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_kagome_minus_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        b1 = [0.0, 0.0]
        b2 = [1/sqrt(3.0), 0.0]
        basis = Array[
            b1 .+ 1/4 * (b2 .- b1),
            b1 .+ 1/4 * (b2 .- a2),
            b1 .+ 1/4 * (b2 .- a1),
            b2 .+ 1/4 * (b1 .- b2),
            b2 .+ 1/4 * (a2 .- b2),
            b2 .+ 1/4 * (a1 .- b2)
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; "tzp"; (0,0)],
            [2; 1; "tzp"; (0,0)],
            [1; 3; "txp"; (0,0)],
            [3; 1; "txp"; (0,0)],
            [2; 3; "typ"; (0,0)],
            [3; 2; "typ"; (0,0)],
            [4; 5; "tzp"; (0,0)],
            [5; 4; "tzp"; (0,0)],
            [4; 6; "txp"; (0,0)],
            [6; 4; "txp"; (0,0)],
            [5; 6; "typ"; (0,0)],
            [6; 5; "typ"; (0,0)],
            [1; 4; "ty"; (0,0)],
            [4; 1; "ty"; (0,0)],
            [2; 5; "tx"; (0,-1)],
            [5; 2; "tx"; (0,1)],
            [3; 6; "tz"; (-1,0)],
            [6; 3; "tz"; (1,0)],
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_kagome_minus_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellKagomeMinus

#-----------------------------------------------------------------------------------------------------------------------------
# HONEYCOMB-XXX LATTICE
# 1 - simple, 11 sites per UC (symmetric around x axis, gives ZZ edge in strip, all couplings identical)
# 2 - simple, 11 sites per UC (symmetric around x axis, gives ZZ edge in strip, couplings fine tuned to be square root)
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellHoneycombXXX(version=1; save=true)
    # distinguish by version
    if version == 1
       # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0,0.0],
            [0.57735,0.0],
            [0.288675,0.0],
            [-0.144338,0.25],
            [-0.144338,-0.25],
            [0.144338,0.0],
            [0.433013,0.0],
            [-0.0721688,0.125],
            [-0.216506,0.375],
            [-0.0721688,-0.125],
            [-0.216506,-0.375]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 6; 1.0; (0,0)],
            [6; 1; 1.0; (0,0)],
            [6; 3; 1.0; (0,0)],
            [3; 6; 1.0; (0,0)],
            [3; 7; 1.0; (0,0)],
            [7; 3; 1.0; (0,0)],
            [7; 2; 1.0; (0,0)],
            [2; 7; 1.0; (0,0)],
            [1; 8; 1.0; (0,0)],
            [8; 1; 1.0; (0,0)],
            [8; 4; 1.0; (0,0)],
            [4; 8; 1.0; (0,0)],
            [4; 9; 1.0; (0,0)],
            [9; 4; 1.0; (0,0)],
            [9; 2; 1.0; (-1,0)],
            [2; 9; 1.0; (1,0)],
            [1; 10; 1.0; (0,0)],
            [10; 1; 1.0; (0,0)],
            [10; 5; 1.0; (0,0)],
            [5; 10; 1.0; (0,0)],
            [5; 11; 1.0; (0,0)],
            [11; 5; 1.0; (0,0)],
            [11; 2; 1.0; (0,-1)],
            [2; 11; 1.0; (0,1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_honeycomb-XXX_unitcell.jld"
    elseif version == 2
       # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0,0.0],
            [0.57735,0.0],
            [0.288675,0.0],
            [-0.144338,0.25],
            [-0.144338,-0.25],
            [0.144338,0.0],
            [0.433013,0.0],
            [-0.0721688,0.125],
            [-0.216506,0.375],
            [-0.0721688,-0.125],
            [-0.216506,-0.375]
        ]
        # Connection Definition
        a = sqrt(sqrt(2.0/3.0))
        b = sqrt(sqrt(3.0/2.0))
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 6; a; (0,0)],
            [6; 1; a; (0,0)],
            [6; 3; b; (0,0)],
            [3; 6; b; (0,0)],
            [3; 7; b; (0,0)],
            [7; 3; b; (0,0)],
            [7; 2; a; (0,0)],
            [2; 7; a; (0,0)],
            [1; 8; a; (0,0)],
            [8; 1; a; (0,0)],
            [8; 4; b; (0,0)],
            [4; 8; b; (0,0)],
            [4; 9; b; (0,0)],
            [9; 4; b; (0,0)],
            [9; 2; a; (-1,0)],
            [2; 9; a; (1,0)],
            [1; 10; a; (0,0)],
            [10; 1; a; (0,0)],
            [10; 5; b; (0,0)],
            [5; 10; b; (0,0)],
            [5; 11; b; (0,0)],
            [11; 5; b; (0,0)],
            [11; 2; a; (0,-1)],
            [2; 11; a; (0,1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_honeycomb-XXX_a_b_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellHoneycombXXX





#-----------------------------------------------------------------------------------------------------------------------------
#
#   Individual FUNCTIONS FOR 3D UNITCELLS
#
#   They all generate the special filenames for the unitcells and return the unitcells as objects
#   Functions can be given a version integer to distinguish several implementations of the same lattice
#   Functions can be specified to already save the unitcell to file
#
#-----------------------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------------------
# DIAMOND LATTICE
# 1 - simple, 2 sites per UC, all connections have interaction strength J1
# 2 - 2 sites per UC, additional next-nearest neighbor connections, interaction strengths are J1 and J2
# 3 - 2 sites per UC, additional anisotropic next-nearest neighbor connections, interaction strengths are J1 and J21, J22
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellDiamond(version=1; save=true, J1=1.0, J2=0.5, J21="J21", J22="J22")
    if version == 1
        # the lattice vectors
        a1 = 0.5 .* [0, 1, 1]
        a2 = 0.5 .* [1, 0, 1]
        a3 = 0.5 .* [1, 1, 0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [0.25, 0.25, 0.25]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)],
            [1; 2; J1; (-1, 0, 0)],
            [1; 2; J1; (0, -1, 0)],
            [1; 2; J1; (0, 0, -1)],

            [2; 1; J1; (0, 0, 0)],
            [2; 1; J1; (1, 0, 0)],
            [2; 1; J1; (0, 1, 0)],
            [2; 1; J1; (0, 0, 1)],
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_diamond_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_diamond_$(J1)_unitcell.jld"
        end
    elseif version == 2
        # the lattice vectors
        a1 = 0.5 .* [0, 1, 1]
        a2 = 0.5 .* [1, 0, 1]
        a3 = 0.5 .* [1, 1, 0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [0.25, 0.25, 0.25]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)],
            [1; 2; J1; (-1, 0, 0)],
            [1; 2; J1; (0, -1, 0)],
            [1; 2; J1; (0, 0, -1)],

            [2; 1; J1; (0, 0, 0)],
            [2; 1; J1; (1, 0, 0)],
            [2; 1; J1; (0, 1, 0)],
            [2; 1; J1; (0, 0, 1)],


            [1; 1; J2; (1, 0, 0)],
            [1; 1; J2; (-1, 0, 0)],
            [1; 1; J2; (0, 1, 0)],
            [1; 1; J2; (0, -1, 0)],
            [1; 1; J2; (0, 0, 1)],
            [1; 1; J2; (0, 0, -1)],

            [1; 1; J2; (1, -1, 0)],
            [1; 1; J2; (-1, 1, 0)],
            [1; 1; J2; (0, 1, -1)],
            [1; 1; J2; (0, -1, 1)],
            [1; 1; J2; (-1, 0, 1)],
            [1; 1; J2; (1, 0, -1)],


            [2; 2; J2; (1, 0, 0)],
            [2; 2; J2; (-1, 0, 0)],
            [2; 2; J2; (0, 1, 0)],
            [2; 2; J2; (0, -1, 0)],
            [2; 2; J2; (0, 0, 1)],
            [2; 2; J2; (0, 0, -1)],

            [2; 2; J2; (1, -1, 0)],
            [2; 2; J2; (-1, 1, 0)],
            [2; 2; J2; (0, 1, -1)],
            [2; 2; J2; (0, -1, 1)],
            [2; 2; J2; (-1, 0, 1)],
            [2; 2; J2; (1, 0, -1)]
        ]
        # filename
        if J1==1.0 && J2==0.5
            filename = "$(FOLDER_UNITCELLS)3d_diamond_NN_NNN_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_diamond_NN_NNN_$(J1)_$(J2)_unitcell.jld"
        end
    elseif version == 3
        # the lattice vectors
        a1 = 0.5 .* [0, 1, 1]
        a2 = 0.5 .* [1, 0, 1]
        a3 = 0.5 .* [1, 1, 0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [0.25, 0.25, 0.25]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)],
            [1; 2; J1; (-1, 0, 0)],
            [1; 2; J1; (0, -1, 0)],
            [1; 2; J1; (0, 0, -1)],

            [2; 1; J1; (0, 0, 0)],
            [2; 1; J1; (1, 0, 0)],
            [2; 1; J1; (0, 1, 0)],
            [2; 1; J1; (0, 0, 1)],


            [1; 1; J21; (1, 0, 0)],
            [1; 1; J21; (-1, 0, 0)],
            [1; 1; J21; (0, 1, 0)],
            [1; 1; J21; (0, -1, 0)],
            [1; 1; J22; (0, 0, 1)],
            [1; 1; J22; (0, 0, -1)],

            [1; 1; J22; (1, -1, 0)],
            [1; 1; J22; (-1, 1, 0)],
            [1; 1; J21; (0, 1, -1)],
            [1; 1; J21; (0, -1, 1)],
            [1; 1; J21; (-1, 0, 1)],
            [1; 1; J21; (1, 0, -1)],


            [2; 2; J21; (1, 0, 0)],
            [2; 2; J21; (-1, 0, 0)],
            [2; 2; J21; (0, 1, 0)],
            [2; 2; J21; (0, -1, 0)],
            [2; 2; J22; (0, 0, 1)],
            [2; 2; J22; (0, 0, -1)],

            [2; 2; J22; (1, -1, 0)],
            [2; 2; J22; (-1, 1, 0)],
            [2; 2; J21; (0, 1, -1)],
            [2; 2; J21; (0, -1, 1)],
            [2; 2; J21; (-1, 0, 1)],
            [2; 2; J21; (1, 0, -1)]
        ]
        # filename
        if J1==1.0 && J21=="J21" && J22=="J22"
            filename = "$(FOLDER_UNITCELLS)3d_diamond_aniso_NNN_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_diamond_aniso_NNN_$(J1)_$(J21)_$(J22)_unitcell.jld"
        end
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellDiamond

#-----------------------------------------------------------------------------------------------------------------------------
# BCC LATTICE
# 1 - simple, 2 sites per UC, all connections have interaction strength J1
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellBCC(version=1; save=true, J1=1.0)
    if version == 1
        # the lattice vectors
        a1 = [1, 0, 0]
        a2 = [0, 1, 0]
        a3 = [0, 0, 1]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)],
            [1; 2; J1; (-1, 0, 0)],
            [1; 2; J1; (0, -1, 0)],
            [1; 2; J1; (-1, -1, 0)],
            [1; 2; J1; (0, 0, -1)],
            [1; 2; J1; (-1, 0, -1)],
            [1; 2; J1; (0, -1, -1)],
            [1; 2; J1; (-1, -1, -1)],

            [2; 1; J1; (0, 0, 0)],
            [2; 1; J1; (1, 0, 0)],
            [2; 1; J1; (0, 1, 0)],
            [2; 1; J1; (1, 1, 0)],
            [2; 1; J1; (0, 0, 1)],
            [2; 1; J1; (1, 0, 1)],
            [2; 1; J1; (0, 1, 1)],
            [2; 1; J1; (1, 1, 1)]
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_bcc_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_bcc_$(J1)_unitcell.jld"
        end
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellBCC

#-----------------------------------------------------------------------------------------------------------------------------
# PYROCHLORE LATTICE
# 1 - simple, 4 sites per UC, all connections have interaction strength J1
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellPyrochlore(version=1; save=true, J1=1.0)
    if version == 1
        # the lattice vectors
        a1 = [0, 0.5, 0.5]
        a2 = [0.5, 0, 0.5]
        a3 = [0.5, 0.5, 0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0., 0., 0.],
            [0., 0.25, 0.25],
            [0.25, 0., 0.25],
            [0.25, 0.25, 0.]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)],
            [1; 3; J1; (0, 0, 0)],
            [1; 4; J1; (0, 0, 0)],
            [2; 1; J1; (0, 0, 0)],
            [2; 3; J1; (0, 0, 0)],
            [2; 4; J1; (0, 0, 0)],
            [3; 1; J1; (0, 0, 0)],
            [3; 2; J1; (0, 0, 0)],
            [3; 4; J1; (0, 0, 0)],
            [4; 1; J1; (0, 0, 0)],
            [4; 2; J1; (0, 0, 0)],
            [4; 3; J1; (0, 0, 0)],

            [1; 4; J1; (0, 0, -1)],
            [4; 1; J1; (0, 0, 1)],
            [1; 2; J1; (-1, 0, 0)],
            [2; 1; J1; (1, 0, 0)],
            [1; 3; J1; (0, -1, 0)],
            [3; 1; J1; (0, 1, 0)],

            [2; 3; J1; (1, -1, 0)],
            [3; 2; J1; (-1, 1, 0)],
            [2; 4; J1; (1, 0, -1)],
            [4; 2; J1; (-1, 0, 1)],

            [3; 4; J1; (0, 1, -1)],
            [4; 3; J1; (0, -1, 1)]
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_pyrochlore_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_pyrochlore_$(J1)_unitcell.jld"
        end
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellPyrochlore


#-----------------------------------------------------------------------------------------------------------------------------
# LATTICE (8,3)a
# 1 - simple, 6 sites per UC, all connections have interaction strength J1
# 4 - Kitaev
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcell_8_3_a(version=1; save=true, J1=1.0)
    if version == 1
        # the lattice vectors
        a1 = [1.0, 0.0, 0.0]
        a2 = [-0.5, sqrt(3)/2., 0.0]
        a3 = [0.0, 0.0, (3*sqrt(2))/5.]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.5, sqrt(3)/10., 0.0],
            [3/5., sqrt(3)/5., (2*sqrt(2))/5.],
            [0.1, (3*sqrt(3))/10., sqrt(2)/5.],
            [0.4, sqrt(3)/5., sqrt(2)/5.],
            [0.0, (2*sqrt(3))/5., 0.0],
            [-0.1, (3*sqrt(3))/10., (2*sqrt(2))/5.]   
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 4; J1; (0, 0, 0)], 
            [4; 2; J1; (0, 0, 0)], # zz
            [4; 3; J1; (0, 0, 0)], 
            [5; 3; J1; (0, 0, 0)], 
            [3; 6; J1; (0, 0, 0)], # zz

            [4; 1; J1; (0, 0, 0)], 
            [2; 4; J1; (0, 0, 0)], # zz
            [3; 4; J1; (0, 0, 0)], 
            [3; 5; J1; (0, 0, 0)], 
            [6; 3; J1; (0, 0, 0)], # zz

            [5; 1; J1; (0, 1, 0)], # zz
            [1; 5; J1; (0, -1, 0)], # zz

            [2; 6; J1; (1, 0, 0)], 
            [6; 2; J1; (-1, 0, 0)],

            [1; 2; J1; (0, 0, -1)], 
            [2; 1; J1; (0, 0, 1)], 
            
            [5; 6; J1; (0, 0, -1)], 
            [6; 5; J1; (0, 0, 1)]           
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_8_3_a_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_8_3_a_$(J1)_unitcell.jld"
        end
    elseif version == 4
        # the lattice vectors
        a1 = [1.0, 0.0, 0.0]
        a2 = [-0.5, sqrt(3)/2., 0.0]
        a3 = [0.0, 0.0, (3*sqrt(2))/5.]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.5, sqrt(3)/10., 0.0],
            [3/5., sqrt(3)/5., (2*sqrt(2))/5.],
            [0.1, (3*sqrt(3))/10., sqrt(2)/5.],
            [0.4, sqrt(3)/5., sqrt(2)/5.],
            [0.0, (2*sqrt(3))/5., 0.0],
            [-0.1, (3*sqrt(3))/10., (2*sqrt(2))/5.]   
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 4; "ty"; (0, 0, 0)], 
            [4; 2; "tz"; (0, 0, 0)], # zz
            [4; 3; "tx"; (0, 0, 0)], 
            [5; 3; "ty"; (0, 0, 0)], 
            [3; 6; "tz"; (0, 0, 0)], # zz

            [4; 1; "ty"; (0, 0, 0)], 
            [2; 4; "tz"; (0, 0, 0)], # zz
            [3; 4; "tx"; (0, 0, 0)], 
            [3; 5; "ty"; (0, 0, 0)], 
            [6; 3; "tz"; (0, 0, 0)], # zz

            [5; 1; "tz"; (0, 1, 0)], # zz
            [1; 5; "tz"; (0, -1, 0)], # zz

            [2; 6; "ty"; (1, 0, 0)], 
            [6; 2; "ty"; (-1, 0, 0)],

            [1; 2; "tx"; (0, 0, -1)], 
            [2; 1; "tx"; (0, 0, 1)], 
            
            [5; 6; "tx"; (0, 0, -1)], 
            [6; 5; "tx"; (0, 0, 1)]             
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_8_3_a_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcell_8_3_a

#-----------------------------------------------------------------------------------------------------------------------------
# LATTICE (8,3)b
# 1 - simple, 6 sites per UC, all connections have interaction strength J1
# 4 - Kitaev
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcell_8_3_b(version=1; save=true, J1=1.0)
    if version == 1
        # the lattice vectors
        a1 = [1/2., 1/(2*sqrt(3)), sqrt(2)/(5*sqrt(3))]
        a2 = [0, 1/sqrt(3), (2*sqrt(2))/(5*sqrt(3))]
        a3 = [0, 0, sqrt(6)/5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [1/10., 1/(2*sqrt(3)), sqrt(2)/(5*sqrt(3))],
            [1/5., sqrt(3)/5, sqrt(6)/5],
            [3/10., 11/(10*sqrt(3)), (4*sqrt(2))/(5*sqrt(3))],
            [1/5., 2/(5*sqrt(3)), (2*sqrt(2))/(5*sqrt(3))],
            [3/10., (3*sqrt(3))/10., sqrt(6)/5],
            [2/5., 1/sqrt(3), sqrt(2)/sqrt(3)]   
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 4; J1; (0, 0, 0)], # zz
            [4; 2; J1; (0, 0, 0)], 
            [2; 5; J1; (0, 0, 0)], # zz
            [5; 3; J1; (0, 0, 0)], 
            [3; 6; J1; (0, 0, 0)], # zz

            [4; 1; J1; (0, 0, 0)], # zz
            [2; 4; J1; (0, 0, 0)], 
            [5; 2; J1; (0, 0, 0)], # zz
            [3; 5; J1; (0, 0, 0)], 
            [6; 3; J1; (0, 0, 0)], # zz

            [6; 1; J1; (1, 0, 1)], 
            [1; 6; J1; (-1, 0, -1)], 

            [4; 3; J1; (0, -1, 0)], 
            [3; 4; J1; (0, 1, 0)], 

            [1; 2; J1; (0, 0, -1)], 
            [2; 1; J1; (0, 0, 1)], 

            [5; 6; J1; (0, 0, -1)], 
            [6; 5; J1; (0, 0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_8_3_b_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [1/2., 1/(2*sqrt(3)), sqrt(2)/(5*sqrt(3))]
        a2 = [0, 1/sqrt(3), (2*sqrt(2))/(5*sqrt(3))]
        a3 = [0, 0, sqrt(6)/5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [1/10., 1/(2*sqrt(3)), sqrt(2)/(5*sqrt(3))],
            [1/5., sqrt(3)/5, sqrt(6)/5],
            [3/10., 11/(10*sqrt(3)), (4*sqrt(2))/(5*sqrt(3))],
            [1/5., 2/(5*sqrt(3)), (2*sqrt(2))/(5*sqrt(3))],
            [3/10., (3*sqrt(3))/10., sqrt(6)/5],
            [2/5., 1/sqrt(3), sqrt(2)/sqrt(3)]   
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 4; "tz"; (0, 0, 0)], # zz
            [4; 2; "ty"; (0, 0, 0)], 
            [2; 5; "tz"; (0, 0, 0)], # zz
            [5; 3; "ty"; (0, 0, 0)], 
            [3; 6; "tz"; (0, 0, 0)], # zz

            [4; 1; "tz"; (0, 0, 0)], # zz
            [2; 4; "ty"; (0, 0, 0)], 
            [5; 2; "tz"; (0, 0, 0)], # zz
            [3; 5; "ty"; (0, 0, 0)], 
            [6; 3; "tz"; (0, 0, 0)], # zz

            [6; 1; "ty"; (1, 0, 1)], 
            [1; 6; "ty"; (-1, 0, -1)], 

            [4; 3; "tx"; (0, -1, 0)], 
            [3; 4; "tx"; (0, 1, 0)], 

            [1; 2; "tx"; (0, 0, -1)], 
            [2; 1; "tx"; (0, 0, 1)], 

            [5; 6; "tx"; (0, 0, -1)], 
            [6; 5; "tx"; (0, 0, 1)]    
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_8_3_b_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    saveUnitcell(uc)
    # return the unitcell
    return uc
end
export getUnitcell_8_3_b

#-----------------------------------------------------------------------------------------------------------------------------
# LATTICE (8,3)c
# 1 - simple, 8 sites per UC, all connections have interaction strength J1
# 4 - Kitaev
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcell_8_3_c(version=1; save=true, J1=1.0)
    if version == 1
        # the lattice vectors
        a1 = [1., 0., 0.]
        a2 = [-1/2., sqrt(3)/2., 0.]
        a3 = [0., 0., 2/5.]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [-1/5., 4/(5*sqrt(3)), 1/10.],  #1
            [ 0., 7/(5*sqrt(3)), 1/10.],    #2		
            [ 1/5., 4/(5*sqrt(3)), 1/10.],  #3
            [ 1/2., 1/(2*sqrt(3)), 3/10.],  #4
            [ 0., 1/sqrt(3), 1/10.],        #5		
            [ 3/10., 7/(10*sqrt(3)), 3/10.],#6
            [ 1/2., 1/(10*sqrt(3)), 3/10.], #7
            [ 7/10., 7/(10*sqrt(3)), 3/10.],#8    
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 5; J1; (0, 0, 0)],
            [2; 5; J1; (0, 0, 0)], # zz
            [5; 3; J1; (0, 0, 0)],
            [3; 6; J1; (0, 0, 0)], # zz
            [6; 4; J1; (0, 0, 0)],
            [4; 7; J1; (0, 0, 0)], # zz
            [4; 8; J1; (0, 0, 0)],

            [5; 1; J1; (0, 0, 0)],
            [5; 2; J1; (0, 0, 0)], # zz
            [3; 5; J1; (0, 0, 0)],
            [6; 3; J1; (0, 0, 0)], # zz
            [4; 6; J1; (0, 0, 0)],
            [7; 4; J1; (0, 0, 0)], # zz
            [8; 4; J1; (0, 0, 0)],

            [8; 1; J1; (1, 0, 0)],
            [1; 8; J1; (-1, 0, 0)],

            [8; 1; J1; (1, 0, 1)],   # zz
            [1; 8; J1; (-1, 0, -1)], # zz 

            [2; 7; J1; (0, 1, 0)],
            [7; 2; J1; (0, -1, 0)],

            [2; 7; J1; (0, 1, -1)],
            [7; 2; J1; (0, -1, 1)],

            [3; 6; J1; (0, 0, -1)],
            [6; 3; J1; (0, 0, 1)]        
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_8_3_c_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_8_3_c_$(J1)_unitcell.jld"
        end
    elseif version == 4
        # the lattice vectors
        a1 = [1., 0., 0.]
        a2 = [-1/2., sqrt(3)/2., 0.]
        a3 = [0., 0., 2/5.]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [-1/5., 4/(5*sqrt(3)), 1/10.], #1
            [ 0., 7/(5*sqrt(3)), 1/10.], #2		
            [ 1/5., 4/(5*sqrt(3)), 1/10.], #3
            [ 1/2., 1/(2*sqrt(3)), 3/10.], #4
            [ 0., 1/sqrt(3), 1/10.], #5		
            [ 3/10., 7/(10*sqrt(3)), 3/10.], #6
            [ 1/2., 1/(10*sqrt(3)), 3/10.], #7
            [ 7/10., 7/(10*sqrt(3)), 3/10.], #8    
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 5; "tx"; (0, 0, 0)],
            [2; 5; "tz"; (0, 0, 0)], # zz
            [5; 3; "ty"; (0, 0, 0)],
            [3; 6; "tz"; (0, 0, 0)], # zz
            [6; 4; "ty"; (0, 0, 0)],
            [4; 7; "tz"; (0, 0, 0)], # zz
            [4; 8; "tx"; (0, 0, 0)],

            [5; 1; "tx"; (0, 0, 0)],
            [5; 2; "tz"; (0, 0, 0)], # zz
            [3; 5; "ty"; (0, 0, 0)],
            [6; 3; "tz"; (0, 0, 0)], # zz
            [4; 6; "ty"; (0, 0, 0)],
            [7; 4; "tz"; (0, 0, 0)], # zz
            [8; 4; "tx"; (0, 0, 0)],

            [8; 1; "ty"; (1, 0, 0)],
            [1; 8; "ty"; (-1, 0, 0)],

            [8; 1; "tz"; (1, 0, 1)],   # zz
            [1; 8; "tz"; (-1, 0, -1)], # zz 

            [2; 7; "tx"; (0, 1, 0)],
            [7; 2; "tx"; (0, -1, 0)],

            [2; 7; "ty"; (0, 1, -1)],
            [7; 2; "ty"; (0, -1, 1)],

            [3; 6; "tx"; (0, 0, -1)],
            [6; 3; "tx"; (0, 0, 1)]        
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_8_3_c_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcell_8_3_c

#-----------------------------------------------------------------------------------------------------------------------------
# LATTICE (8,3)n
# 1 - simple, 16 sites per UC, all connections have interaction strength J1
# 4 - Kitaev
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcell_8_3_n(version=1; save=true, J1=1.0)
    if version == 1
        # the lattice vectors
        a = [1.0, 0.0, 0.0]
        b = [0.0, 1.0, 0.0]
        c = [0.0, 0.0, 4/(2*sqrt(3) + sqrt(2))]
        a1 = a
        a2 = b
        a3 = 0.5*(a+b+c)
        x = (sqrt(3) + sqrt(2))/(2*(2*sqrt(3) + sqrt(2)))
        z = 0.125
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            x*a + (0.5 - x)*b + c/4.,
            (1-x)*a + (0.5 - x)*b + c/4.,
            (0.5 + x)*a + b/2. + (0.5 - z)*c,
            (1-x)*a + (0.5 + x)*b + c/4.,
            x*a + (0.5 + x)*b + c/4.,
            (0.5 - x)*a + b/2. + (0.5 - z)*c,
            (1-x)*b + z*c,
            x*b + z*c,
            (0.5 - x)*a + x*b + c/4.,
            a/2. + (0.5 - x)*b + (0.5 - z)*c,
            (0.5 + x)*a + x*b + c/4.,
            (0.5 + x)*a + (1 - x)*b + c/4.,
            a/2. + (0.5 + x)*b + (0.5 - z)*c,
            (0.5 - x)*a + (1 - x)*b + c/4.,
            x*a + z*c,
            (1-x)*a + z*c		
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 10; J1; (0, 0, 0)], 
            [10; 2; J1; (0, 0, 0)], 
            [2; 11; J1; (0, 0, 0)], # zz
            [11; 3; J1; (0, 0, 0)], 
            [3; 12; J1; (0, 0, 0)], 
            [12; 4; J1; (0, 0, 0)], # zz
            [4; 13; J1; (0, 0, 0)], 
            [13; 5; J1; (0, 0, 0)], 
            [5; 14; J1; (0, 0, 0)], # zz
            [14; 6; J1; (0, 0, 0)], 
            [6; 9; J1; (0, 0, 0)], 
            [9; 8; J1; (0, 0, 0)], 
            [9; 1; J1; (0, 0, 0)], # zz
            [1; 15; J1; (0, 0, 0)], 
            [2; 16; J1; (0, 0, 0)], 
            [14; 7; J1; (0, 0, 0)], 

            [10; 1; J1; (0, 0, 0)], 
            [2; 10; J1; (0, 0, 0)], 
            [11; 2; J1; (0, 0, 0)], # zz
            [3; 11; J1; (0, 0, 0)], 
            [12; 3; J1; (0, 0, 0)], 
            [4; 12; J1; (0, 0, 0)], # zz
            [13; 4; J1; (0, 0, 0)], 
            [5; 13; J1; (0, 0, 0)], 
            [14; 5; J1; (0, 0, 0)], # zz
            [6; 14; J1; (0, 0, 0)], 
            [9; 6; J1; (0, 0, 0)], 
            [8; 9; J1; (0, 0, 0)], 
            [1; 9; J1; (0, 0, 0)], # zz
            [15; 1; J1; (0, 0, 0)], 
            [16; 2; J1; (0, 0, 0)], 
            [7; 14; J1; (0, 0, 0)], 

            [11; 8; J1; (1, 0, 0)],
            [8; 11; J1; (-1, 0, 0)],

            [12; 7; J1; (1, 0, 0)],
            [7; 12; J1; (-1, 0, 0)],

            [7; 10; J1; (0, 1, -1)], # zz
            [10; 7; J1; (0, -1, 1)], # zz

            [8; 13; J1; (0, 0, -1)], # zz
            [13; 8; J1; (0, 0, 1)],  # zz

            [16; 4; J1; (0, -1, 0)],
            [4; 16; J1; (0, 1, 0)],

            [15; 5; J1; (0, -1, 0)],
            [5; 15; J1; (0, 1, 0)],

            [15; 3; J1; (0, 0, -1)], # zz
            [3; 15; J1; (0, 0, 1)],  # zz

            [16; 6; J1; (1, 0, -1)], # zz
            [6; 16; J1; (-1, 0, 1)]  # zz            
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_8_3_n_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_8_3_n_$(J1)_unitcell.jld"
        end
    elseif version == 4
        # the lattice vectors
        a = [1.0, 0.0, 0.0]
        b = [0.0, 1.0, 0.0]
        c = [0.0, 0.0, 4/(2*sqrt(3) + sqrt(2))]
        a1 = a
        a2 = b
        a3 = 0.5*(a+b+c)
        x = (sqrt(3) + sqrt(2))/(2*(2*sqrt(3) + sqrt(2)))
        z = 0.125
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            x*a + (0.5 - x)*b + c/4.,
            (1-x)*a + (0.5 - x)*b + c/4.,
            (0.5 + x)*a + b/2. + (0.5 - z)*c,
            (1-x)*a + (0.5 + x)*b + c/4.,
            x*a + (0.5 + x)*b + c/4.,
            (0.5 - x)*a + b/2. + (0.5 - z)*c,
            (1-x)*b + z*c,
            x*b + z*c,
            (0.5 - x)*a + x*b + c/4.,
            a/2. + (0.5 - x)*b + (0.5 - z)*c,
            (0.5 + x)*a + x*b + c/4.,
            (0.5 + x)*a + (1 - x)*b + c/4.,
            a/2. + (0.5 + x)*b + (0.5 - z)*c,
            (0.5 - x)*a + (1 - x)*b + c/4.,
            x*a + z*c,
            (1-x)*a + z*c
            
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 10; "tx"; (0, 0, 0)], 
            [10; 2; "ty"; (0, 0, 0)], 
            [2; 11; "tz"; (0, 0, 0)], # zz
            [11; 3; "tx"; (0, 0, 0)], 
            [3; 12; "ty"; (0, 0, 0)], 
            [12; 4; "tz"; (0, 0, 0)], # zz
            [4; 13; "tx"; (0, 0, 0)], 
            [13; 5; "ty"; (0, 0, 0)], 
            [5; 14; "tz"; (0, 0, 0)], # zz
            [14; 6; "tx"; (0, 0, 0)], 
            [6; 9; "ty"; (0, 0, 0)], 
            [9; 8; "tx"; (0, 0, 0)], 
            [9; 1; "tz"; (0, 0, 0)], # zz
            [1; 15; "ty"; (0, 0, 0)], 
            [2; 16; "tx"; (0, 0, 0)], 
            [14; 7; "ty"; (0, 0, 0)], 

            [10; 1; "tx"; (0, 0, 0)], 
            [2; 10; "ty"; (0, 0, 0)], 
            [11; 2; "tz"; (0, 0, 0)], # zz
            [3; 11; "tx"; (0, 0, 0)], 
            [12; 3; "ty"; (0, 0, 0)], 
            [4; 12; "tz"; (0, 0, 0)], # zz
            [13; 4; "tx"; (0, 0, 0)], 
            [5; 13; "ty"; (0, 0, 0)], 
            [14; 5; "tz"; (0, 0, 0)], # zz
            [6; 14; "tx"; (0, 0, 0)], 
            [9; 6; "ty"; (0, 0, 0)], 
            [8; 9; "tx"; (0, 0, 0)], 
            [1; 9; "tz"; (0, 0, 0)], # zz
            [15; 1; "ty"; (0, 0, 0)], 
            [16; 2; "tx"; (0, 0, 0)], 
            [7; 14; "ty"; (0, 0, 0)], 

            [11; 8; "ty"; (1, 0, 0)],
            [8; 11; "ty"; (-1, 0, 0)],

            [12; 7; "tx"; (1, 0, 0)],
            [7; 12; "tx"; (-1, 0, 0)],

            [7; 10; "tz"; (0, 1, -1)], # zz
            [10; 7; "tz"; (0, -1, 1)], # zz

            [8; 13; "tz"; (0, 0, -1)], # zz
            [13; 8; "tz"; (0, 0, 1)],  # zz

            [16; 4; "ty"; (0, -1, 0)],
            [4; 16; "ty"; (0, 1, 0)],

            [15; 5; "tx"; (0, -1, 0)],
            [5; 15; "tx"; (0, 1, 0)],

            [15; 3; "tz"; (0, 0, -1)], # zz
            [3; 15; "tz"; (0, 0, 1)],  # zz

            [16; 6; "tz"; (1, 0, -1)], # zz
            [6; 16; "tz"; (-1, 0, 1)] # zz            
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_8_3_n_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcell_8_3_n


#-----------------------------------------------------------------------------------------------------------------------------
# LATTICE (9,3)a
# 1 - simple, 12 sites per UC, all connections have interaction strength J1
# 2 - simple shifted
# 4 - Kitaev
# 5 - Kitaev shifted
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcell_9_3_a(version=1; save=true, J1=1.0)
    if version==1
        # the lattice vectors
        a = [1.0, 0.0, 0.0]
        b = [-0.5, sqrt(3)/2., 0.0]
        c = [0.0, 0.0, sqrt(6*(4 + sqrt(3)))/(1 + 2*sqrt(3))]
        a1 = -a/3. + b/3. + c/3.
        a2 = -a/3. -2*b/3. + c/3.
        a3 = 2*a/3. + b/3. + c/3.
        d_f = sqrt(3)/(1+2*sqrt(3))
        d_h = (29 - 3*sqrt(3))/132.
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            d_f * a,
            2*d_h * a + d_h * b + c/12.,
            d_f * (a + b),
            d_h * a + 2*d_h * b - c/12.,
            d_f * b,
            -d_h *a + d_h * b + c/12.,
            -d_f * a,
            -2*d_h * a - d_h * b - c/12.,
            -d_f * (a + b),
            -d_h * a - 2*d_h*b + c/12.,
            -d_f * b,
            d_h * a - d_h * b - c/12.
            
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)], 
            [2; 3; J1; (0, 0, 0)], 
            [3; 4; J1; (0, 0, 0)], 
            [4; 5; J1; (0, 0, 0)], 
            [5; 6; J1; (0, 0, 0)], 
            [6; 7; J1; (0, 0, 0)], 
            [7; 8; J1; (0, 0, 0)],
            [8; 9; J1; (0, 0, 0)], 
            [9; 10; J1; (0, 0, 0)], 
            [10; 11; J1; (0, 0, 0)], 
            [11; 12; J1; (0, 0, 0)], 
            [12; 1; J1; (0, 0, 0)], 

            [2; 1; J1; (0, 0, 0)], 
            [3; 2; J1; (0, 0, 0)], 
            [4; 3; J1; (0, 0, 0)], 
            [5; 4; J1; (0, 0, 0)], 
            [6; 5; J1; (0, 0, 0)], 
            [7; 6; J1; (0, 0, 0)], 
            [8; 7; J1; (0, 0, 0)],
            [9; 8; J1; (0, 0, 0)], 
            [10; 9; J1; (0, 0, 0)], 
            [11; 10; J1; (0, 0, 0)], 
            [12; 11; J1; (0, 0, 0)], 
            [1; 12; J1; (0, 0, 0)], 


            [3; 9; J1; (0, -1, 1)], # zz
            [9; 3; J1; (0, 1, -1)], # zz

            [1; 7; J1; (-1, 0, 1)], # zz
            [7; 1; J1; (1, 0, -1)], # zz

            [5; 11; J1; (1, -1, 0)], # zz
            [11; 5; J1; (-1, 1, 0)], # zz
            
            [12; 6; J1; (-1, 0, 0)], # zz
            [6; 12; J1; (1, 0, 0)], # zz

            [8; 2; J1; (0, 0, -1)], # zz
            [2; 8; J1; (0, 0, 1)], # zz
    
            [4; 10; J1; (0, -1, 0)], # zz
            [10; 4; J1; (0, 1, 0)]   # zz 
            
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_9_3_a_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_9_3_a_$(J1)_unitcell.jld"
        end
    elseif version == 2
        # the lattice vectors
        a1 = [-sqrt(3)/2., 1/2., 1/sqrt(3)]
        a2 = [0, -1, 1/sqrt(3)]
        a3 = [sqrt(3)/2, 1/2, 1/sqrt(3)]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0, 0, 0],
            a1/6. - a2/6.,
            a1/3. - a2/3.,
            a1/2. - a2/3. - a3/6.,
            2/3.*a1 - a2/3. - a3/3.,
            2/3.*a1 - a2/6. - a3/2.,
            2/3.*a1 - 2/3.*a3,
            a1/2. + a2/6. - 2/3.*a3,
            a1/3. + a2/3. - 2/3.*a3,
            a1/6. + a2/3. - a3/2.,
            a2/3. - a3/3.,
            a2/6. - a3/6.		
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)], 
            [2; 3; J1; (0, 0, 0)], 
            [3; 4; J1; (0, 0, 0)], 
            [4; 5; J1; (0, 0, 0)], 
            [5; 6; J1; (0, 0, 0)], 
            [6; 7; J1; (0, 0, 0)], 
            [7; 8; J1; (0, 0, 0)],
            [8; 9; J1; (0, 0, 0)], 
            [9; 10; J1; (0, 0, 0)], 
            [10; 11; J1; (0, 0, 0)], 
            [11; 12; J1; (0, 0, 0)], 
            [12; 1; J1; (0, 0, 0)], 

            [2; 1; J1; (0, 0, 0)], 
            [3; 2; J1; (0, 0, 0)], 
            [4; 3; J1; (0, 0, 0)], 
            [5; 4; J1; (0, 0, 0)], 
            [6; 5; J1; (0, 0, 0)], 
            [7; 6; J1; (0, 0, 0)], 
            [8; 7; J1; (0, 0, 0)],
            [9; 8; J1; (0, 0, 0)], 
            [10; 9; J1; (0, 0, 0)], 
            [11; 10; J1; (0, 0, 0)], 
            [12; 11; J1; (0, 0, 0)], 
            [1; 12; J1; (0, 0, 0)], 

            [1; 7; J1; (-1, 0, 1)],
            [7; 1; J1; (1, 0, -1)],
            [2; 8; J1; (0, 0, 1)],
            [8; 2; J1; (0, 0, -1)],
            [3; 9; J1; (0, -1, 1)],
            [9; 3; J1; (0, 1, -1)],
            [4; 10; J1; (0, -1, 0)],
            [10; 4; J1; (0, 1, 0)],
            [5; 11; J1; (1, -1, 0)],
            [11; 5; J1; (-1, 1, 0)],
            [6; 12; J1; (1, 0, 0)],
            [12; 6; J1; (-1, 0, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_9_3_a_V2_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a = [1.0, 0.0, 0.0]
        b = [-0.5, sqrt(3)/2., 0.0]
        c = [0.0, 0.0, sqrt(6*(4 + sqrt(3)))/(1 + 2*sqrt(3))]
        a1 = -a/3. + b/3. + c/3.
        a2 = -a/3. -2*b/3. + c/3.
        a3 = 2*a/3. + b/3. + c/3.
        d_f = sqrt(3)/(1+2*sqrt(3))
        d_h = (29 - 3*sqrt(3))/132.
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            d_f * a,
            2*d_h * a + d_h * b + c/12.,
            d_f * (a + b),
            d_h * a + 2*d_h * b - c/12.,
            d_f * b,
            -d_h *a + d_h * b + c/12.,
            -d_f * a,
            -2*d_h * a - d_h * b - c/12.,
            -d_f * (a + b),
            -d_h * a - 2*d_h*b + c/12.,
            -d_f * b,
            d_h * a - d_h * b - c/12.
            
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; "tx"; (0, 0, 0)], 
            [2; 3; "ty"; (0, 0, 0)], 
            [3; 4; "tx"; (0, 0, 0)], 
            [4; 5; "ty"; (0, 0, 0)], 
            [5; 6; "tx"; (0, 0, 0)], 
            [6; 7; "ty"; (0, 0, 0)], 
            [7; 8; "tx"; (0, 0, 0)],
            [8; 9; "ty"; (0, 0, 0)], 
            [9; 10; "tx"; (0, 0, 0)], 
            [10; 11; "ty"; (0, 0, 0)], 
            [11; 12; "tx"; (0, 0, 0)], 
            [12; 1; "ty"; (0, 0, 0)], 

            [2; 1; "tx"; (0, 0, 0)], 
            [3; 2; "ty"; (0, 0, 0)], 
            [4; 3; "tx"; (0, 0, 0)], 
            [5; 4; "ty"; (0, 0, 0)], 
            [6; 5; "tx"; (0, 0, 0)], 
            [7; 6; "ty"; (0, 0, 0)], 
            [8; 7; "tx"; (0, 0, 0)],
            [9; 8; "ty"; (0, 0, 0)], 
            [10; 9; "tx"; (0, 0, 0)], 
            [11; 10; "ty"; (0, 0, 0)], 
            [12; 11; "tx"; (0, 0, 0)], 
            [1; 12; "ty"; (0, 0, 0)], 


            [3; 9; "tz"; (0, -1, 1)], # zz
            [9; 3; "tz"; (0, 1, -1)], # zz

            [1; 7; "tz"; (-1, 0, 1)], # zz
            [7; 1; "tz"; (1, 0, -1)], # zz

            [5; 11; "tz"; (1, -1, 0)], # zz
            [11; 5; "tz"; (-1, 1, 0)], # zz
            
            [12; 6; "tz"; (-1, 0, 0)], # zz
            [6; 12; "tz"; (1, 0, 0)], # zz

            [8; 2; "tz"; (0, 0, -1)], # zz
            [2; 8; "tz"; (0, 0, 1)], # zz
    
            [4; 10; "tz"; (0, -1, 0)], # zz
            [10; 4; "tz"; (0, 1, 0)]   # zz 
            
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_9_3_a_kitaev_unitcell.jld"
    elseif version == 5
        # the lattice vectors
        a1 = [-sqrt(3)/2., 1/2., 1/sqrt(3)]
        a2 = [0, -1, 1/sqrt(3)]
        a3 = [sqrt(3)/2, 1/2, 1/sqrt(3)]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0, 0, 0],
            a1/6. - a2/6.,
            a1/3. - a2/3.,
            a1/2. - a2/3. - a3/6.,
            2/3.*a1 - a2/3. - a3/3.,
            2/3.*a1 - a2/6. - a3/2.,
            2/3.*a1 - 2/3.*a3,
            a1/2. + a2/6. - 2/3.*a3,
            a1/3. + a2/3. - 2/3.*a3,
            a1/6. + a2/3. - a3/2.,
            a2/3. - a3/3.,
            a2/6. - a3/6.		
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; "tx"; (0, 0, 0)], 
            [2; 3; "ty"; (0, 0, 0)], 
            [3; 4; "tx"; (0, 0, 0)], 
            [4; 5; "ty"; (0, 0, 0)], 
            [5; 6; "tx"; (0, 0, 0)], 
            [6; 7; "ty"; (0, 0, 0)], 
            [7; 8; "tx"; (0, 0, 0)],
            [8; 9; "ty"; (0, 0, 0)], 
            [9; 10; "tx"; (0, 0, 0)], 
            [10; 11; "ty"; (0, 0, 0)], 
            [11; 12; "tx"; (0, 0, 0)], 
            [12; 1; "ty"; (0, 0, 0)], 

            [2; 1; "tx"; (0, 0, 0)], 
            [3; 2; "ty"; (0, 0, 0)], 
            [4; 3; "tx"; (0, 0, 0)], 
            [5; 4; "ty"; (0, 0, 0)], 
            [6; 5; "tx"; (0, 0, 0)], 
            [7; 6; "ty"; (0, 0, 0)], 
            [8; 7; "tx"; (0, 0, 0)],
            [9; 8; "ty"; (0, 0, 0)], 
            [10; 9; "tx"; (0, 0, 0)], 
            [11; 10; "ty"; (0, 0, 0)], 
            [12; 11; "tx"; (0, 0, 0)], 
            [1; 12; "ty"; (0, 0, 0)], 

            [1; 7; "tz"; (-1, 0, 1)],
            [7; 1; "tz"; (1, 0, -1)],
            [2; 8; "tz"; (0, 0, 1)],
            [8; 2; "tz"; (0, 0, -1)],
            [3; 9; "tz"; (0, -1, 1)],
            [9; 3; "tz"; (0, 1, -1)],
            [4; 10; "tz"; (0, -1, 0)],
            [10; 4; "tz"; (0, 1, 0)],
            [5; 11; "tz"; (1, -1, 0)],
            [11; 5; "tz"; (-1, 1, 0)],
            [6; 12; "tz"; (1, 0, 0)],
            [12; 6; "tz"; (-1, 0, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_9_3_a_V2_kitaev_unitcell.jld"    
    end    
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcell_9_3_a


#-----------------------------------------------------------------------------------------------------------------------------
# HYPEROCTAGON LATTICE (10,3)a
# 1 - simple, 4 sites per UC, all connections have interaction strength J1
# 4 - Kitaev
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellHyperoctagon(version=1; save=true, J1=1.0)
    if version == 1
        # the lattice vectors
        a1 = [1, 0, 0]
        a2 = [0.5, 0.5, -0.5]
        a3 = [0.5, 0.5, 0.5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.125, 0.125, 0.125],
            [5/8., 3/8., -1/8.],
            [3/8., 1/8., -1/8.],
            [7/8., 3/8., 1/8.]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 3; J1; (0, 0, 0)],
            [3; 2; J1; (0, 0, 0)],  # zz
            [2; 4; J1; (0, 0, 0)],

            [3; 1; J1; (0, 0, 0)],
            [2; 3; J1; (0, 0, 0)],  # zz
            [4; 2; J1; (0, 0, 0)],

            [4; 1; J1; (1, 0, 0)],  # zz
            [1; 4; J1; (-1, 0, 0)], # zz

            [2; 1; J1; (0, 1, 0)],
            [1; 2; J1; (0, -1, 0)],

            [3; 4; J1; (0, 0, -1)],
            [4; 3; J1; (0, 0, 1)],
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_hyperoctagon_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_hyperoctagon_$(J1)_unitcell.jld"
        end
    elseif version == 4
        # the lattice vectors
        a1 = [1, 0, 0]
        a2 = [0.5, 0.5, -0.5]
        a3 = [0.5, 0.5, 0.5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.125, 0.125, 0.125],
            [5/8., 3/8., -1/8.],
            [3/8., 1/8., -1/8.],
            [7/8., 3/8., 1/8.]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 3; "tx"; (0, 0, 0)],
            [3; 2; "tz"; (0, 0, 0)],  # zz
            [2; 4; "tx"; (0, 0, 0)],

            [3; 1; "tx"; (0, 0, 0)],
            [2; 3; "tz"; (0, 0, 0)],  # zz
            [4; 2; "tx"; (0, 0, 0)],

            [4; 1; "tz"; (1, 0, 0)],  # zz
            [1; 4; "tz"; (-1, 0, 0)], # zz

            [2; 1; "ty"; (0, 1, 0)],
            [1; 2; "ty"; (0, -1, 0)],

            [3; 4; "ty"; (0, 0, -1)],
            [4; 3; "ty"; (0, 0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_hyperoctagon_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
function getUnitcell_10_3_a(version=1; save=true, J1=1.0)
    return getUnitcellHyperoctagon(version, save=save, J1=J1)
end
export getUnitcellHyperoctagon
export getUnitcell_10_3_a

#-----------------------------------------------------------------------------------------------------------------------------
# HYPERHONEYCOMB LATTICE (10,3)b
# 1 - simple, 4 sites per UC, all connections have interaction strength J1
# 2 - simple, 4 sites per UC, all connections have interaction strength J1, site 4 shifted through the unitcell
# 4 - Kitaev
# 5 - Kitaev shifted unitcell
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellHyperhoneycomb(version=1; save=true, J1=1.0)
    if version == 1
        # the lattice vectors
        a1 = [-1, 1, -2]
        a2 = [-1, 1, 2]
        a3 = [2, 4, 0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 2.0, 1.0],
            [0.0, -1.0, 1.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)],
            [1; 4; J1; (0, 0, 0)],
            [1; 4; J1; (1, 0, 0)],

            [2; 1; J1; (0, 0, 0)],
            [2; 3; J1; (0, 0, 0)],
            [2; 3; J1; (0, -1, 0)],

            [3; 2; J1; (0, 0, 0)],
            [3; 2; J1; (0, 1, 0)],
            [3; 4; J1; (0, 0, 1)],

            [4; 1; J1; (0, 0, 0)],
            [4; 1; J1; (-1, 0, 0)],
            [4; 3; J1; (0, 0, -1)]        
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_hyperhoneycomb_v1_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_hyperhoneycomb_v1_$(J1)_unitcell.jld"
        end
    elseif version == 2
        # the lattice vectors
        a1 = [-1, 1, -2]
        a2 = [-1, 1, 2]
        a3 = [2, 4, 0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 2.0, 1.0],
            [2.0, 3.0, 1.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)],
            [1; 4; J1; (0, 0, -1)],
            [1; 4; J1; (1, 0, -1)],

            [2; 1; J1; (0, 0, 0)],
            [2; 3; J1; (0, 0, 0)],
            [2; 3; J1; (0, -1, 0)],

            [3; 2; J1; (0, 0, 0)],
            [3; 2; J1; (0, 1, 0)],
            [3; 4; J1; (0, 0, 0)],

            [4; 1; J1; (0, 0, 1)],
            [4; 1; J1; (-1, 0, 1)],
            [4; 3; J1; (0, 0, 0)]        
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_hyperhoneycomb_v2_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_hyperhoneycomb_v2_$(J1)_unitcell.jld"
        end
    elseif version == 4
        # the lattice vectors
        a1 = [-1, 1, -2]
        a2 = [-1, 1, 2]
        a3 = [2, 4, 0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 2.0, 1.0],
            [0.0, -1.0, 1.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; "tz"; (0, 0, 0)],
            [1; 4; "tx"; (0, 0, 0)],
            [1; 4; "ty"; (1, 0, 0)],

            [2; 1; "tz"; (0, 0, 0)],
            [2; 3; "tx"; (0, 0, 0)],
            [2; 3; "ty"; (0, -1, 0)],

            [3; 2; "tx"; (0, 0, 0)],
            [3; 2; "ty"; (0, 1, 0)],
            [3; 4; "tz"; (0, 0, 1)],

            [4; 1; "tx"; (0, 0, 0)],
            [4; 1; "ty"; (-1, 0, 0)],
            [4; 3; "tz"; (0, 0, -1)]        
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_hyperhoneycomb_kitaev_unitcell.jld"
    elseif version == 5
        # the lattice vectors
        a1 = [-1, 1, -2]
        a2 = [-1, 1, 2]
        a3 = [2, 4, 0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 2.0, 1.0],
            [2.0, 3.0, 1.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; "tz"; (0, 0, 0)],
            [1; 4; "tx"; (0, 0, -1)],
            [1; 4; "ty"; (1, 0, -1)],

            [2; 1; "tz"; (0, 0, 0)],
            [2; 3; "tx"; (0, 0, 0)],
            [2; 3; "ty"; (0, -1, 0)],

            [3; 2; "tx"; (0, 0, 0)],
            [3; 2; "ty"; (0, 1, 0)],
            [3; 4; "tz"; (0, 0, 0)],

            [4; 1; "tx"; (0, 0, 1)],
            [4; 1; "ty"; (-1, 0, 1)],
            [4; 3; "tz"; (0, 0, 0)]        
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_hyperhoneycomb_v2_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
function getUnitcell_10_3_b(version=1; save=true, J1=1.0)
    return getUnitcellHyperhoneycomb(version, save=save, J1=J1)
end
export getUnitcellHyperhoneycomb
export getUnitcell_10_3_b

#-----------------------------------------------------------------------------------------------------------------------------
# LATTICE (10,3)c
# 1 - simple, 6 sites per UC, all connections have interaction strength J1
# 4 - Kitaev
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcell_10_3_c(version=1; save=true, J1=1.0)
    if version == 1
        # the lattice vectors
        a1 = [1, 0, 0]
        a2 = [-0.5, sqrt(3)/2., 0.0]
        a3 = [0.0, 0.0, (3*sqrt(3))/2.]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.25, 1/(4*sqrt(3)), 1/(2*sqrt(3))],
            [0.75, 1/(4*sqrt(3)), 2/sqrt(3)],
            [0.5, 1/sqrt(3), 7/(2*sqrt(3))],
            [0.75, 1/(4*sqrt(3)), 1/sqrt(3)],
            [0.5, 1/sqrt(3), 5/(2*sqrt(3))],
            [0.25, 1/(4*sqrt(3)), 4/sqrt(3)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 4; J1; (0, 0, 0)],
            [4; 2; J1; (0, 0, 0)],  # zz
            [2; 5; J1; (0, 0, 0)],
            [5; 3; J1; (0, 0, 0)],  # zz
            [3; 6; J1; (0, 0, 0)],

            [4; 1; J1; (0, 0, 0)],
            [2; 4; J1; (0, 0, 0)],  # zz
            [5; 2; J1; (0, 0, 0)],
            [3; 5; J1; (0, 0, 0)],  # zz
            [6; 3; J1; (0, 0, 0)],

            [4; 1; J1; (1, 0, 0)],  
            [1; 4; J1; (-1, 0, 0)], 

            [5; 2; J1; (0, 1, 0)], 
            [2; 5; J1; (0, -1, 0)],

            [3; 6; J1; (1, 1, 0)], 
            [6; 3; J1; (-1, -1, 0)], 

            [1; 6; J1; (0, 0, -1)],  # zz
            [6; 1; J1; (0, 0, 1)], # zz
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_10_3_c_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_10_3_c_$(J1)_unitcell.jld"
        end
    elseif version == 4
        # the lattice vectors
        a1 = [1, 0, 0]
        a2 = [-0.5, sqrt(3)/2., 0.0]
        a3 = [0.0, 0.0, (3*sqrt(3))/2.]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.25, 1/(4*sqrt(3)), 1/(2*sqrt(3))],
            [0.75, 1/(4*sqrt(3)), 2/sqrt(3)],
            [0.5, 1/sqrt(3), 7/(2*sqrt(3))],
            [0.75, 1/(4*sqrt(3)), 1/sqrt(3)],
            [0.5, 1/sqrt(3), 5/(2*sqrt(3))],
            [0.25, 1/(4*sqrt(3)), 4/sqrt(3)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 4; "tx"; (0, 0, 0)],
            [4; 2; "tz"; (0, 0, 0)],  # zz
            [2; 5; "tx"; (0, 0, 0)],
            [5; 3; "tz"; (0, 0, 0)],  # zz
            [3; 6; "tx"; (0, 0, 0)],

            [4; 1; "tx"; (0, 0, 0)],
            [2; 4; "tz"; (0, 0, 0)],  # zz
            [5; 2; "tx"; (0, 0, 0)],
            [3; 5; "tz"; (0, 0, 0)],  # zz
            [6; 3; "tx"; (0, 0, 0)],

            [4; 1; "ty"; (1, 0, 0)],  
            [1; 4; "ty"; (-1, 0, 0)], 

            [5; 2; "ty"; (0, 1, 0)], 
            [2; 5; "ty"; (0, -1, 0)],

            [3; 6; "ty"; (1, 1, 0)], 
            [6; 3; "ty"; (-1, -1, 0)], 

            [1; 6; "tz"; (0, 0, -1)],  # zz
            [6; 1; "tz"; (0, 0, 1)] # zz
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_10_3_c_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcell_10_3_c

#-----------------------------------------------------------------------------------------------------------------------------
# LATTICE (10,3)d
# 1 - simple, 8 sites per UC, all connections have interaction strength J1
# 4 - Kitaev
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcell_10_3_d(version=1; save=true, J1=1.0)
    if version == 1
        # the lattice vectors
        a = 0.25*(2 - sqrt(2))
        c = 0.5
        a1 = [1, 0, 0]
        a2 = [0.5, 0.5, 0.0]
        a3 = [0.0, 0.0, c]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1 - a2)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, -a, 0.75*c],
            [-a, 0.0, 0.5*c],
            [0, a, 0.25*c],
            [a, 0.0, 0.0],
            [-a, -0.5, 0.25*c],
            [0, a - 0.5, 0.5*c],
            [a, -0.5, 0.75*c],
            [0, -a - 0.5, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)],
            [3; 4; J1; (0, 0, 0)],  
            [5; 6; J1; (0, 0, 0)],  

            [2; 1; J1; (0, 0, 0)],
            [4; 3; J1; (0, 0, 0)],  
            [6; 5; J1; (0, 0, 0)],

            [2; 3; J1; (0, 0, 0)], 
            [3; 2; J1; (0, 0, 0)], 

            [6; 7; J1; (0, 0, 0)], 
            [7; 6; J1; (0, 0, 0)],

            [1; 6; J1; (0, 0, 0)], #zz 
            [6; 1; J1; (0, 0, 0)], #zz

            [4; 5; J1; (0, 1, 0)], #zz
            [5; 4; J1; (0, -1, 0)], #zz

            [8; 5; J1; (0, 0, 0)], 
            [5; 8; J1; (0, 0, 0)], 

            [8; 7; J1; (0, 0, -1)], 
            [7; 8; J1; (0, 0, 1)], 

            [2; 7; J1; (-1, 0, 0)], # zz
            [7; 2; J1; (1, 0, 0)],  # zz

            [3; 8; J1; (-1, 1, 0)], # zz
            [8; 3; J1; (1, -1, 0)],  # zz

            [1; 4; J1; (0, 0, 1)], 
            [4; 1; J1; (0, 0, -1)], 

        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_10_3_d_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_10_3_d_$(J1)_unitcell.jld"
        end
    elseif version == 4
        # the lattice vectors
        a = 0.25*(2 - sqrt(2))
        c = 0.5
        a1 = [1, 0, 0]
        a2 = [0.5, 0.5, 0.0]
        a3 = [0.0, 0.0, c]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1 - a2)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, -a, 0.75*c],
            [-a, 0.0, 0.5*c],
            [0, a, 0.25*c],
            [a, 0.0, 0.0],
            [-a, -0.5, 0.25*c],
            [0, a - 0.5, 0.5*c],
            [a, -0.5, 0.75*c],
            [0, -a - 0.5, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; "tx"; (0, 0, 0)],
            [3; 4; "tx"; (0, 0, 0)],  
            [5; 6; "ty"; (0, 0, 0)],  

            [2; 1; "tx"; (0, 0, 0)],
            [4; 3; "tx"; (0, 0, 0)],  
            [6; 5; "ty"; (0, 0, 0)],

            [2; 3; "ty"; (0, 0, 0)], 
            [3; 2; "ty"; (0, 0, 0)], 

            [6; 7; "tx"; (0, 0, 0)], 
            [7; 6; "tx"; (0, 0, 0)],

            [1; 6; "tz"; (0, 0, 0)], #zz 
            [6; 1; "tz"; (0, 0, 0)], #zz

            [4; 5; "tz"; (0, 1, 0)], #zz
            [5; 4; "tz"; (0, -1, 0)], #zz

            [8; 5; "tx"; (0, 0, 0)], 
            [5; 8; "tx"; (0, 0, 0)], 

            [8; 7; "ty"; (0, 0, -1)], 
            [7; 8; "ty"; (0, 0, 1)], 

            [2; 7; "tz"; (-1, 0, 0)], # zz
            [7; 2; "tz"; (1, 0, 0)],  # zz

            [3; 8; "tz"; (-1, 1, 0)], # zz
            [8; 3; "tz"; (1, -1, 0)],  # zz

            [1; 4; "ty"; (0, 0, 1)], 
            [4; 1; "ty"; (0, 0, -1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_10_3_d_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcell_10_3_d






#-----------------------------------------------------------------------------------------------------------------------------
#
#   GET A UNITCELL OBJECT BY ONLY PASSING SITES AND LATTICE VECTORS
#
#   Parameters (required)
#   - sites: an array indicating the position of sites
#   - lattice_vectors: an array indicating the lattice vectors of the lattice
#
#   Parameters (optional)
#   - max_ij(k): the maximum number of unitcells in all direction that is beeing considered to build geometric neighbors
#   - epsilon: The precision up to which distances are considered to be equal
#   - min_NN: minimal number of nearest neighbors that ALL sites have to have (-1: Auto)
#   - max_NN: maximal number of nearest neighbors that ALL sites have to have (-1: Auto)
#   - strength_NN: interaction strength of the new connections
#   - name: name of the lattice (if "AUTO" is selected, a name will be build automatically)
#   - save: if to save the newly constructed unitcell
#
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellFromSites2D(
        sites,
        lattice_vectors;
        max_ij::Int64=3,
        epsilon::Float64=1e-8,
        min_NN=-1,
        max_NN=-1,
        strength_NN=1.0,
        name::String="AUTO",
        save=false
    )
    # make a new list with connections
    connections = Array[]
    # define what "close to" means
    function distance(p1, p2)
        return sqrt(sum((p1.-p2).*(p1.-p2)))
    end
    function closeto(p1, p2)
        return distance(p1, p2) < epsilon
    end
    # function to sort the checklist
    function sortfunction(checklistentry)
        return checklistentry[3]
    end
    # iterate over all sites
    for s in 1:size(sites, 1)
        # build a new checklist for sites
        checklist = Array[]
        # check all other sites in a lattice around this one
        for i in -max_ij:max_ij
        for j in -max_ij:max_ij
        for alpha in 1:size(sites, 1)
            # site at i*a1 + j*a2 + basis[alpha]
            position = i.*lattice_vectors[1] .+ j.*lattice_vectors[2] .+ sites[alpha]
            # check if the position is close to the original position
            if closeto(position, sites[s])
                # just the original site
                continue
            end
            # put the site into the list
            # format: [site_index, (wrap), distance]
            push!(checklist, [alpha, (i,j), distance(position, sites[s])])
        end
        end
        end
        # sort the checklist
        sort!(checklist, by=sortfunction)
        # evaluate the checklist
        # first: find out the number of nearest neighbors
        number_NN = 0
        # if unrestricted: simply go by distance
        if min_NN <= 0 && max_NN <= 0
            for c in checklist
                if closeto(c[3], checklist[1][3])
                    number_NN = number_NN + 1
                else
                    break
                end
            end
        elseif min_NN > 0 && max_NN <=0
            number_NN = min_NN-1
            for c in checklist[min_NN:end]
                if closeto(c[3], checklist[min_NN][3])
                    number_NN = number_NN + 1
                else
                    break
                end
            end
        else   
            if min_NN <= 0
                min_NN = 1
            end
            number_NN = min_NN-1
            for c in checklist[min_NN:max_NN]
                if closeto(c[3], checklist[min_NN][3])
                    number_NN = number_NN + 1
                else
                    break
                end
            end
        end
        println("number of nearest neighbors of site $(s): $(number_NN)")
        if closeto(checklist[number_NN][3],checklist[number_NN+1][3])
            println("cut at the wrong number of nearest neighbors, there are more that are equally close!")
        end
        # restrict the checklist to that length
        checklist = checklist[1:number_NN]
        # parse into connections
        for c in checklist
            # define a new connection of interaction strength type float or int
            if typeof(strength_NN) == Float64 || typeof(strength_NN) == Int64
                push!(connections, [s, c[1], strength_NN, c[2]])
            elseif typeof(strength_NN) == String
                # automatic connection geometry based strength
                if strength_NN == "AUTO"
                    push!(connections, [s, c[1], 1/c[3], c[2]])
                end
            else
                println("unknown connection strength type: $(typeof(strength_NN))")
                return
            end
        end
    end
    # correct all connections - i.e. find out if the returning connection is present, if not: add it
    connections_additional = Array[]
    for c in connections
        # check if "-c" is present
        found = false
        for mc in connections
            if c[1] == mc[2] && c[2] == mc[1] && c[3] == mc[3] && -c[4][1] == mc[4][1] && -c[4][2] == mc[4][2]
                found = true
            end
        end
        if !found
            push!(connections_additional, [c[2], c[1], c[3], (-c[4][1], -c[4][2])])
        end
    end
    # add all additional connections
    for c in connections_additional
        push!(connections, c)
    end
    # build a filename
    if name == "AUTO"
        name = "from_$(size(sites, 1))_sites"
    end
    filename = "$(FOLDER_UNITCELLS)2d_$(name)_unitcell.jld"
    # build the finished unitcell object
    unitcell = Unitcell(
        lattice_vectors,
        sites,
        connections,
        filename
    )
    # maybe save
    if save
        saveUnitcell(unitcell)
    end
    # return the unitcell
    return unitcell
end
export getUnitcellFromSites2D

function getUnitcellFromSites3D(
        sites,
        lattice_vectors;
        max_ijk::Int64=3,
        epsilon::Float64=1e-8,
        min_NN=-1,
        max_NN=-1,
        strength_NN=1.0,
        name::String="AUTO",
        save=false
    )
    # make a new list with connections
    connections = Array[]
    # define what "close to" means
    function distance(p1, p2)
        return sqrt(sum((p1.-p2).*(p1.-p2)))
    end
    function closeto(p1, p2)
        return distance(p1, p2) < epsilon
    end
    # function to sort the checklist
    function sortfunction(checklistentry)
        return checklistentry[3]
    end
    # iterate over all sites
    for s in 1:size(sites, 1)
        # build a new checklist for sites
        checklist = Array[]
        # check all other sites in a lattice around this one
        for i in -max_ijk:max_ijk
        for j in -max_ijk:max_ijk
        for k in -max_ijk:max_ijk
        for alpha in 1:size(sites, 1)
            # site at i*a1 + j*a2 + k*a3 + basis[alpha]
            position = i.*lattice_vectors[1] .+ j.*lattice_vectors[2] .+ k.*lattice_vectors[3] .+ sites[alpha]
            # check if the position is close to the original position
            if closeto(position, sites[s])
                # just the original site
                continue
            end
            # put the site into the list
            # format: [site_index, (wrap), distance]
            push!(checklist, [alpha, (i,j,k), distance(position, sites[s])])
        end
        end
        end
        end
        # sort the checklist
        sort!(checklist, by=sortfunction)
        # evaluate the checklist
        # first: find out the number of nearest neighbors
        number_NN = 0
        # if unrestricted: simply go by distance
        if min_NN <= 0 && max_NN <= 0
            for c in checklist
                if closeto(c[3], checklist[1][3])
                    number_NN = number_NN + 1
                else
                    break
                end
            end
        elseif min_NN > 0 && max_NN <=0
            number_NN = min_NN-1
            for c in checklist[min_NN:end]
                if closeto(c[3], checklist[min_NN][3])
                    number_NN = number_NN + 1
                else
                    break
                end
            end
        else   
            if min_NN <= 0
                min_NN = 1
            end
            number_NN = min_NN-1
            for c in checklist[min_NN:max_NN]
                if closeto(c[3], checklist[min_NN][3])
                    number_NN = number_NN + 1
                else
                    break
                end
            end
        end
        println("number of nearest neighbors of site $(s): $(number_NN)")
        if closeto(checklist[number_NN][3],checklist[number_NN+1][3])
            println("cut at the wrong number of nearest neighbors, there are more that are equally close!")
        end
        # restrict the checklist to that length
        checklist = checklist[1:number_NN]
        # parse into connections
        for c in checklist
            # define a new connection of interaction strength type float or int
            if typeof(strength_NN) == Float64 || typeof(strength_NN) == Int64
                push!(connections, [s, c[1], strength_NN, c[2]])
            elseif typeof(strength_NN) == String
                # automatic connection geometry based strength
                if strength_NN == "AUTO"
                    push!(connections, [s, c[1], 1/c[3], c[2]])
                end
            else
                println("unknown connection strength type: $(typeof(strength_NN))")
                return
            end
        end
    end
    # correct all connections - i.e. find out if the returning connection is present, if not: add it
    connections_additional = Array[]
    for c in connections
        # check if "-c" is present
        found = false
        for mc in connections
            if c[1] == mc[2] && c[2] == mc[1] && c[3] == mc[3] && -c[4][1] == mc[4][1] && -c[4][2] == mc[4][2] && -c[4][3] == mc[4][3]
                found = true
            end
        end
        if !found
            push!(connections_additional, [c[2], c[1], c[3], (-c[4][1], -c[4][2], -c[4][3])])
        end
    end
    # add all additional connections
    for c in connections_additional
        push!(connections, c)
    end
    # build a filename
    if name == "AUTO"
        name = "from_$(size(sites, 1))_sites"
    end
    filename = "$(FOLDER_UNITCELLS)3d_$(name)_unitcell.jld"
    # build the finished unitcell object
    unitcell = Unitcell(
        lattice_vectors,
        sites,
        connections,
        filename
    )
    # maybe save
    if save
        saveUnitcell(unitcell)
    end
    # return the unitcell
    return unitcell
end
export getUnitcellFromSites3D

function getUnitcellFromSites(
        sites,
        lattice_vectors;
        max_ijk::Int64=3,
        epsilon::Float64=1e-8,
        min_NN=-1,
        max_NN=-1,
        strength_NN=1.0,
        name::String="AUTO",
        save=false
    )
    # check the dimension and call the respective function
    if size(lattice_vectors, 1) == 1
        println("building from sites with one lattice vector not supported")
    elseif size(lattice_vectors, 1) == 2
        return getUnitcellFromSites2D(
            sites,
            lattice_vectors;
            max_ij=max_ijk,
            epsilon=epsilon,
            min_NN=min_NN,
            max_NN=max_NN,
            strength_NN=strength_NN,
            name=name,
            save=save
        )
    elseif size(lattice_vectors, 1) == 3
        return getUnitcellFromSites3D(
            sites,
            lattice_vectors;
            max_ijk=max_ijk,
            epsilon=epsilon,
            min_NN=min_NN,
            max_NN=max_NN,
            strength_NN=strength_NN,
            name=name,
            save=save
        )
    elseif size(lattice_vectors, 1) > 3
        println("building from sites with more than 3 lattice vector not supported")
    end
end
export getUnitcellFromSites






#-----------------------------------------------------------------------------------------------------------------------------
#
#   UNITCELL GENERATING CODE
#   Input a Unitcell and optionally a version and the function prints the generating code snippet to the std.out stream
#   version specifies if the code contains a version distinction and which version the printed code has
#
#-----------------------------------------------------------------------------------------------------------------------------
function printUnitcellGeneratingCode(unitcell::Unitcell; version=-1)
    # distinguish if versioning yes or no
    if version != -1
        # header printing
        println("function getUnitcellCustom(version=$(version); save=true)")
        # version distinguish
        println("   # distinguish by version")
        println("   if version == $(version)")
        # lattice vectors
        println("       # the lattice vectors")
        index = 0
        for l in unitcell.lattice_vectors
            index = index+1
            println("       a$(index) = $(l)")
        end
        println("       lattice_vectors = Array[]")
        for index in 1:size(unitcell.lattice_vectors, 1)
            println("       push!(lattice_vectors, a$(index))")
        end
        # Basis
        println("       # Basis Definition")
        println("       basis = Array[")
        for b in unitcell.basis[1:end-1]
            println("           $(b),")
        end
        println("           $(unitcell.basis[end])")
        println("       ]")
        # Connections
        println("       # Connection Definition")
        println("       # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]")
        println("       connections = Array[")
        for c in unitcell.connections[1:end-1]
            if typeof(c[3]) == String
                println("           [$(c[1]); $(c[2]); \"$(c[3])\"; $(c[4])],")
            else
                println("           [$(c[1]); $(c[2]); $(c[3]); $(c[4])],")
            end
        end
        c = unitcell.connections[end]
        if typeof(c[3]) == String
            println("           [$(c[1]); $(c[2]); \"$(c[3])\"; $(c[4])],")
        else
            println("           [$(c[1]); $(c[2]); $(c[3]); $(c[4])],")
        end
        println("       ]")
        # Filename
        println("       # filename")
        if size(unitcell.lattice_vectors,1) == 2
            println("       filename = \"\$(FOLDER_UNITCELLS)2d_custom_unitcell.jld\"")
        elseif size(unitcell.lattice_vectors,1) == 3
            println("       filename = \"\$(FOLDER_UNITCELLS)3d_custom_unitcell.jld\"")
        else
            println("       filename = \"\$(FOLDER_UNITCELLS)Xd_custom_unitcell.jld\"")
        end
        # end
        println("   end")
    else
        # header printing
        println("function getUnitcellCustom(version=1; save=true)")
        # lattice vectors
        println("   # the lattice vectors")
        index = 0
        for l in unitcell.lattice_vectors
            index = index+1
            println("   a$(index) = $(l)")
        end
        println("   lattice_vectors = Array[]")
        for index in 1:size(unitcell.lattice_vectors, 1)
            println("   push!(lattice_vectors, a$(index))")
        end
        # Basis
        println("   # Basis Definition")
        println("   basis = Array[")
        for b in unitcell.basis[1:end-1]
            println("       $(b),")
        end
        println("       $(unitcell.basis[end])")
        println("   ]")
        # Connections
        println("   # Connection Definition")
        println("   # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]")
        println("   connections = Array[")
        for c in unitcell.connections[1:end-1]
            if typeof(c[3]) == String
                println("       [$(c[1]); $(c[2]); \"$(c[3])\"; $(c[4])],")
            else
                println("       [$(c[1]); $(c[2]); $(c[3]); $(c[4])],")
            end
        end
        c = unitcell.connections[end]
        if typeof(c[3]) == String
            println("       [$(c[1]); $(c[2]); \"$(c[3])\"; $(c[4])],")
        else
            println("       [$(c[1]); $(c[2]); $(c[3]); $(c[4])],")
        end
        println("   ]")
        # Filename
        println("   # filename")
        if size(unitcell.lattice_vectors,1) == 2
            println("   filename = \"\$(FOLDER_UNITCELLS)2d_custom_unitcell.jld\"")
        elseif size(unitcell.lattice_vectors,1) == 3
            println("   filename = \"\$(FOLDER_UNITCELLS)3d_custom_unitcell.jld\"")
        else
            println("   filename = \"\$(FOLDER_UNITCELLS)Xd_custom_unitcell.jld\"")
        end
    end
    # part which is not version dependent
    println("   # generate unitcell")
    println("   uc = Unitcell(lattice_vectors, basis, connections, filename)")
    println("   if save")
    println("       saveUnitcell(uc)")
    println("   end")
    println("   # return the unitcell")
    println("   return uc")
    println("end")
end
export printUnitcellGeneratingCode
































#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#
#   LATTICE DEFINITIONS
#
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------




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
function getLatticePeriodic2D(unitcell::Unitcell, repetition_array::Array{Int64}; save=true, load=false)

    # extract the cardinal directions of the lattice from the array
    N_a1 = repetition_array[1]
    N_a2 = repetition_array[2]

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
function getLatticePeriodic3D(unitcell::Unitcell, repetition_array::Array{Int64}; save=true, load=false)

    # extract the cardinal directions of the lattice from the array
    N_a1 = repetition_array[1]
    N_a2 = repetition_array[2]
    N_a3 = repetition_array[3]

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

export getLatticePeriodic2D
export getLatticePeriodic3D


# FOR UNKNOWN DIMENSION
function getLatticePeriodic(unitcell::Unitcell, repetition_array::Array{Int64}; save=true, load=false)

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
    # return a new lattice
    return Lattice(
        lattice.unitcell,
        lattice.unitcellRepetitions,
        lattice.lattice_vectors,
        lattice.positions,
        lattice.positions_indices,
        connections_new,
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
    # return a new unitcell
    return Unitcell(
        unitcell.lattice_vectors,
        unitcell.basis,
        connections_new,
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
    # overwrite the connections in the given lattice
    lattice.connections = connections_new
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
    # overwrite the connections in the given lattice
    unitcell.connections = connections_new
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
            if contains(c[3], old_strength)
                c[3] = replace(c[3], old_strength, new_strength)
            end
        end
        end
    end
    # maybe even evaluate
    if evaluate
        for c in unitcell.connections
            c[3] = eval(parse(c[3]))
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
            if contains(c[3], old_strength)
                c[3] = replace(c[3], old_strength, new_strength)
            end
        end
        end
    end
    # maybe even evaluate
    if evaluate
        for c in lattice.connections
            c[3] = eval(parse(c[3]))
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












#-----------------------------------------------------------------------------------------------------------------------------
#
#   Find Plaquettes of some site in the lattice with desired length len
#
#-----------------------------------------------------------------------------------------------------------------------------
function getPlaquettesOfSite(lattice::Lattice, site, len)
    # get the connectivity of the lattice
    cl = getConnectionList(lattice)
    # history includes all sites that are already visited, including the current site
    function givePlaquettePaths(history, steps, patharray)
        # check if steps = 0, if yes, return history only
        if length(history) == steps+1
            if history[end] == history[1]
                push!(patharray, history)
            end
            return
        elseif history[end] in history[1:end-1]
            # building a closed loop
            return
        end
        # current site
        currentSite = history[end]
        # get all options
        options = []
        for c in cl[currentSite]
            push!(options, Int(c[2]))
        end
        # distinguish between origin and further apart sites
        if length(history) == 1
            # go for all options
            for o in options
                # get the new history
                history_new = copy(history)
                push!(history_new, o)
                # get all subpaths
                givePlaquettePaths(history_new, steps, patharray)
            end
        else
            # go for all options that dont go back
            for o in options
                # check if going back
                if o == history[end-1]
                    continue
                end
                # get the new history
                history_new = copy(history)
                push!(history_new, o)
                # get all subpaths
                givePlaquettePaths(history_new, steps, patharray)
            end
        end
        # return
        return
    end
    # call the function to generate array of plaquettes
    paths = Array[]
    givePlaquettePaths([site], len, paths)
    # check all paths if they are plaquettes
    plaquettes = Array[]
    for p in paths
        if p[1] == p[end]
            # check if the inverse path is already in the plaquettes list
            if !(p[end:-1:1] in plaquettes) 
                push!(plaquettes, p)
            end
        end
    end
    # return the array
    return plaquettes
end
export getPlaquettesOfSite

function getPlaquettesOfLattice(lattice::Lattice, len)
    # get all plaquettes of all sites
    plaquettes = Array[]
    plaquettes_sorted = Array[]
    for i in 1:size(lattice.positions,1)
        # obtain plaquettes
        plaquettes_site = getPlaquettesOfSite(lattice, i, len)
        # push plaquettes into array
        for p in plaquettes_site
            # check if present
            p_sorted = copy(p[1:end-1])
            sort!(p_sorted)
            # search if present
            found = false
            for ps in plaquettes_sorted
                # check if lists equal
                if findfirst([ps[i] == p_sorted[i] for i in 1:len], false) == 0
                    found = true
                    break
                end
            end
            if !(found)
                push!(plaquettes, p)
                push!(plaquettes_sorted, p_sorted)
            end
        end
    end
    # return the list
    return plaquettes
end
export getPlaquettesOfLattice



function printPlaquetteStatisticsOfLattice(lattice::Lattice; detailed=false, max_length=12)
    # print header information
    println("printing plaquette information for lattice:")
    # iterate over all lengths
    for l in 3:max_length
        # get plaquettes
        plaquettes = getPlaquettesOfLattice(lattice, l)
        # print statistics
        println("- length $(l) bonds: $(size(plaquettes,1)) plaquettes")
        # print more information
        if detailed
            for p in plaquettes
                println("    - $(p)")
            end
        end
    end
end
export printPlaquetteStatisticsOfLattice

function printPlaquetteStatisticsOfSite(lattice::Lattice, site; detailed=false, max_length=12)
    # print header information
    println("printing plaquette information for site $(site) of lattice:")
    # iterate over all lengths
    for l in 3:max_length
        # get plaquettes
        plaquettes = getPlaquettesOfSite(lattice, site, l)
        # print statistics
        println("- length $(l) bonds: $(size(plaquettes,1)) plaquettes")
        # print more information
        if detailed
            for p in plaquettes
                println("    - $(p)")
            end
        end
    end
end
export printPlaquetteStatisticsOfSite

















#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#
#   PLOTTING TO SVG IMAGES
#
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------------------------------
#
#   HELPER METHODS FOR SVG CREATION
#   Get the correct SVG code lines for SVG beginning and end of document
#   Get the correct SVG code lines for certain objects
#
#-----------------------------------------------------------------------------------------------------------------------------

# HEADER STRING (must be first in every SVG file)
# width and height denote the dimensions of the image in px
function getSVGHeaderString(width::Int64, height::Int64)
	headerstring = """<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>
<svg
	xmlns:dc=\"http://purl.org/dc/elements/1.1/\"
	xmlns:cc=\"http://creativecommons.org/ns#\"
	xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"
	xmlns:svg=\"http://www.w3.org/2000/svg\"
	xmlns=\"http://www.w3.org/2000/svg\"
	xmlns:xlink=\"http://www.w3.org/1999/xlink\"
	version=\"1.1\"
	id=\"svg2\"
	viewBox=\"0 0 $(width) $(height)\"
	height=\"$(height)\"
	width=\"$(width)\">
	<defs id=\"defs4\">
	</defs>
	<metadata
		id=\"metadata7\">
		<rdf:RDF>
			<cc:Work rdf:about=\"\">
				<dc:format>image/svg+xml</dc:format>
				<dc:type
					rdf:resource=\"http://purl.org/dc/dcmitype/StillImage\" />
				<dc:title></dc:title>
			</cc:Work>
		</rdf:RDF>
	</metadata>

"""
	return headerstring
end

# FOOTER STRING (must be end of every SVG file)
function getSVGFooterString()
    return """
</svg>
"""
end

# STRING FOR AN ELLIPSE
# Parameters are
# - id: The id which the object has later on
# - centerX / centerY: The coordinates of the ellipse center in the canvas
# - radiusX / radiusY: The radii of the ellipse
# - color: hex string of the ellipse color
# formerly: getEllipseString
function getSVGStringEllipse(id, centerX, centerY, radiusX, radiusY, color; label=987654321, labelcolor="#000000")
	es = """
	<ellipse
		style=\"color:$(color);fill:$(color);fill-opacity:1;fill-rule:nonzero\"
		id=\"$(id)\"
		cx=\"$(centerX)\"
		cy=\"$(centerY)\"
		rx=\"$(radiusX)px\"
		ry=\"$(radiusY)px\" />

"""
    # check if to add a label
    if label != 987654321
        es = """
    $(es)

    <text
    x=\"$(centerX)\"
    y=\"$(centerY+radiusY/2)\"
    style=\"text-anchor: middle; font-size: $(radiusX*1.2)px; fill:$(labelcolor)\"   
    >
    $(label)
    </text>
"""
    end
	return es
end

# STRING FOR A STROKED ELLIPSE
# Parameters are
# - id: The id which the object has later on
# - centerX / centerY: The coordinates of the ellipse center in the canvas
# - radiusX / radiusY: The radii of the ellipse
# - colorFill: hex string of the ellipse fill color
# - colorStroke: hex string of the ellipse stroke color
# - strokewidth: Float or Int giving the width of the surrounding stroke
# formerly: getStrokedEllipseString
function getSVGStringEllipseStroked(id, centerX, centerY, radiusX, radiusY, colorFill, colorStroke, strokewidth; label=987654321, labelcolor="#000000")
	es = """
	<ellipse
		style=\"color:$(colorFill);fill:$(colorFill);fill-opacity:1;fill-rule:nonzero;stroke:$(colorStroke);stroke-width:$(strokewidth);stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-dasharray:none;stroke-dashoffset:0;stroke-opacity:1\"
		id=\"$(id)\"
		cx=\"$(centerX)\"
		cy=\"$(centerY)\"
		rx=\"$(radiusX)px\"
		ry=\"$(radiusY)px\" />

"""
    # check if to add a label
    if label != 987654321
        es = """
    $(es)
    
    <text
    x=\"$(centerX)\"
    y=\"$(centerY+radiusY/2)\"
    style=\"text-anchor: middle; font-size: $(radiusX*1.2)px; fill:$(labelcolor)\"   
    >
    $(label)
    </text>
"""
    end
	return es
end

# STRING FOR A LINE
# Parameters are
# - id: The id which the object has later on
# - from: coordinate array [x,y] of the starting point
# - to:   coordinate array [x,y] of the end point
# - colorStroke: hex string of the color that the line has
# - strokewidth: Float or Int giving the width of the line
# - dashed: is the line dashed or not
# formerly: getLineString
function getSVGStringLine(id, from, to, colorStroke, strokewidth; dashed=false, opacity::Float64=1.0)
	if dashed
        ls = """
	<path
		style=\"fill:none;fill-rule:evenodd;stroke:$(colorStroke);stroke-width:$(strokewidth);stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-linecap:butt;stroke-dasharray:$(strokewidth),$(strokewidth*2);stroke-dashoffset:0;stroke-opacity:$(opacity)\"
		id=\"$(id)\"
		d=\"m $(from[1]),$(from[2]) $(to[1]-from[1]),$(to[2]-from[2])\"/>
"""
    else
        ls = """
	<path
		style=\"fill:none;fill-rule:evenodd;stroke:$(colorStroke);stroke-width:$(strokewidth);stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-dasharray:none;stroke-dashoffset:0;stroke-opacity:$(opacity)\"
		id=\"$(id)\"
		d=\"m $(from[1]),$(from[2]) $(to[1]-from[1]),$(to[2]-from[2])\"/>
"""
    end
	return ls
end


# STRING FOR A PLAQUETTE
# Parameters are
# - id: The id which the object has later on
# - plaquette_points: list of point coordinates of points that form the plaquette
# - colorFill: hex string of the color that the plaquette has
# - opacity: Opacity of the plaquette
function getSVGStringPlaquette(id, plaquette_points, colorFill; opacity=0.5)
    plaquette_string = "M"
    for p in plaquette_points
        plaquette_string = "$(plaquette_string) $(p[1]),$(p[2])"
    end
    plaquette_string = "$(plaquette_string) z"
    ps = """
	<path
		style=\"fill:$(colorFill);fill-opacity:$(opacity);fill-rule:evenodd;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"
		d=\"$(plaquette_string)\"
		id=\"$(id)\" />
"""
end

# CONVERSION OF RGB COLORS TO HEX STRINGS
function color_hex(r::Int64,g::Int64,b::Int64)
    return "#$(hex(r,2))$(hex(g,2))$(hex(b,2))"
end
function color_hex(rgb::Array{Int64})
    return color_hex(rgb[1], rgb[2], rgb[3])
end



# GET AUTOMATICAL SEQUENCE OF COLORS
function getColorSequence(len::Int64)
    # define a new list for the sequence
    colors = Array[]
    # fill the sequence depending on the number of requested elements
    if len <= 1
        # only one color --> black
        push!(colors, [0,0,0])
    elseif len == 2
        # only two colors --> black and red
        push!(colors, [0,0,0])
        push!(colors, [255,0,0])
    elseif len == 3
        # only three colors --> red, green, blue (Kitaev)
        push!(colors, [255,0,0])
        push!(colors, [0,255,0])
        push!(colors, [0,0,255])
    else
        # continuous sequence of random colors
        for i in 1:len
            push!(colors, [rand(1:255),rand(1:255),rand(1:255)])
        end
    end
    # return the sequence
    return colors
end
function getGreySequence(len::Int64)
    # define a new list for the sequence
    colors = Array[]
    # fill the sequence depending on the number of requested elements
    if len <= 1
        # only one color --> black
        push!(colors, [0,0,0])
    elseif len == 2
        # only two colors --> black and grey
        push!(colors, [0,0,0])
        push!(colors, [125,125,125])
    elseif len == 3
        # only three colors --> black and two grey
        push!(colors, [0,0,0])
        push!(colors, [90,90,90])
        push!(colors, [180,180,180])
    else
        # continuous sequence
        for i in 1:len
            c = round(Int64, i/len * 200)
            push!(colors, [c,c,c])
        end
    end
    # return the sequence
    return colors
end

















#-----------------------------------------------------------------------------------------------------------------------------
#
#   SVG PLOTTING OF 2D AND 3D LATTICES
#
#   Generates a SVG file that shows the lattice
#   SVG file can optionally be opened and converted to pdf
#
#   Parameters (necessary 2D and 3D):
#   -   lattice: The lattice object to plot
#
#   Parameters (optional 2D and 3D):
#   -   conversion: The factor of how many pixels represent a length 1 of real space distance (Default: 160)
#   -   border_percentage: percentage of the width which is allocated as border (Default: 0.1)
#   -   filename_output: The output filename, "AUTO" meaning it is generated automatically (Default: "AUTO")
#   -   site_radius: Radius of sites (Default: 25)
#   -   site_border_width_percentage: Border width of sites in percent of the diameter (Default: 0.2)
#   -   bond_thickness: The thickness of bonds (Default: 8)
#   -   visualize_periodic: Decide if periodic connections are also shown (Default: false)
#   -   colorcode_sites: Dictonary containing a colorcode for all sites matching the site index to a RGB color array
#   -   colorcode_bonds: Dictonary containing a colorcode for all bonds matching the interaction strength to a RGB color array
#   -   colorcode_bonds_automation: String specifying if bonds are colored automatically. Select either "COLOR" or "GREY" or "OFF"
#   -   openfile: Decide if the newly created SVG should already be opened by external viewer (Default: false)
#   -   export_pdf: Decide if svg should be converted to pdf (Default: true)
#
#   Parameters (optional but only for 3D):
#   -   bond_color_BG: The color of bonds that are in the farthest back of the plot (i.e. to which color to converge towards the back) 
#   -   site_color_BG: The color of sites that are in the farthest back of the plot (i.e. to which color to converge towards the back) 
#   -   lattice_rotation: Array of Floats indicating the rotation around the XYZ axis
#   -   camera_position_relative: Array of Floats indicating the position of the camera
#   -   DOF: Float indicating the strength of the Depth of Field in color gradient
#
#-----------------------------------------------------------------------------------------------------------------------------
function plotLattice2D(
		lattice::Lattice;
		conversion = 160,
		border_percentage=0.1,
		filename_output::String="AUTO",
		site_radius=25,
		site_border_width_percentage::Float64=0.2,
        site_labels="OFF",
		bond_thickness::Int64=8,
		visualize_periodic=false,
		colorcode_sites = Dict(0 => [255,255,255], 1 => [255,255,255]),
		colorcode_bonds = Dict("0" => [0,0,0], "1.0" => [0,0,0]),
        colorcode_bonds_automation::String = "OFF",
		openfile=false,
        export_pdf=true
		)

    # define the filename_output if it is set to AUTO
    if filename_output=="AUTO"
        filename_output = "$(lattice.filename[1:end-4])_plot.svg"
        filename_output = replace(filename_output, FOLDER_LATTICES, "")
    end

	# load positions and connections
	positions	= lattice.positions
	connections = lattice.connections

    # maybe overwrite dictonary
    if colorcode_bonds_automation == "GREY"
        # construct the list of interaction strengths
        cs_list = getConnectionStrengthList(lattice)
        # get the color code list
        cc_list = getGreySequence(size(cs_list, 1))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[string(cs_list[i])] = cc_list[i]
        end
    elseif colorcode_bonds_automation == "COLOR"
        # construct the list of interaction strengths
        cs_list = getConnectionStrengthList(lattice)
        # get the color code list
        cc_list = getColorSequence(size(cs_list, 1))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[string(cs_list[i])] = cc_list[i]
        end
    end

	# repair color dictonary
	colorcode_bonds["0"] = get(colorcode_bonds, "0", [0,0,0])
	colorcode_sites[0] = get(colorcode_sites, 0, [255,255,255])


	# define styles for the different sites, i.e. strings that are saved into the svg strings
	site_r 				= site_radius
	site_border			= "#000000"
	site_border_width	= "$(site_border_width_percentage*site_radius)px"

	# sites to plot
	pos_x = zeros(size(positions,1))
	pos_y = zeros(size(positions,1))
	for i in 1:size(positions,1)
		pos_x[i] = positions[i][1]
		pos_y[i] = positions[i][2]
	end
	sites_to_plot = positions
	indices_to_plot = lattice.positions_indices
	xvals = pos_x
	yvals = pos_y
	border 			= border_percentage * (maximum(pos_x) + maximum(pos_y) - minimum(pos_x) - minimum(pos_y))/4

	# connections to plot (all)
	connections_to_plot = copy(connections)
    # the neutral connection wrap, i.e. which wrap identifies a non-periodic connection
	neutral_connection_wrap = (0,0)
	if size(lattice.lattice_vectors,1) == 1
		neutral_connection_wrap = (0)
	end
	

	# define the width and height of the canvas
	width_UC 	= (maximum(xvals) + site_radius/conversion + border) - (minimum(xvals) - site_radius/conversion - border)
	height_UC	= (maximum(yvals) + site_radius/conversion + border) - (minimum(yvals) - site_radius/conversion - border)
	min_x	= minimum(xvals) - site_radius/conversion - border
	min_y	= minimum(yvals) - site_radius/conversion - border
	width	= conversion*(width_UC)
	height	= conversion*(height_UC)

	# define the conversion functions for coordinates
	function X(x)
		return conversion*(x - min_x)
	end
	function Y(y)
		return + conversion*(y - min_y)
	end

	# open the SVG file
	file = open(filename_output, "w")

	# write the headerstring
	write(file, getSVGHeaderString(round(Int64, width), round(Int64, height)))

	# write all connections
	for (i,c) in enumerate(connections_to_plot)
		if (c[4] == neutral_connection_wrap)
			connection_color = color_hex(get(colorcode_bonds, string(c[3]), colorcode_bonds["0"]))
			write(file, getSVGStringLine("path$(i)", [X(positions[Int(c[1])][1]), Y(positions[Int(c[1])][2])], [X(positions[Int(c[2])][1]), Y(positions[Int(c[2])][2])], connection_color, bond_thickness))
		elseif visualize_periodic
			connection_color = color_hex(get(colorcode_bonds, string(c[3]), colorcode_bonds["0"]))
			write(file, getSVGStringLine("path$(i)", [X(positions[Int(c[1])][1]), Y(positions[Int(c[1])][2])], [X(positions[Int(c[2])][1]), Y(positions[Int(c[2])][2])], connection_color, bond_thickness, dashed=true))
        end
	end

	# write all sites
	for (i,s) in enumerate(sites_to_plot)
        site_color_basic = get(colorcode_sites, indices_to_plot[i], colorcode_sites[0])
		site_color = color_hex(site_color_basic)
        label_color = color_hex([colorelment < 100 ? 255 : 0 for colorelment in site_color_basic])
        if site_labels == "POSITION INDEX"
		    write(file, getSVGStringEllipseStroked("el$(i)", X(s[1]), Y(s[2]), site_r, site_r, site_color, site_border, site_border_width, label=indices_to_plot[i], labelcolor=label_color))
        elseif site_labels == "LATTICE INDEX"
		    write(file, getSVGStringEllipseStroked("el$(i)", X(s[1]), Y(s[2]), site_r, site_r, site_color, site_border, site_border_width, label=i, labelcolor=label_color))
        else
		    write(file, getSVGStringEllipseStroked("el$(i)", X(s[1]), Y(s[2]), site_r, site_r, site_color, site_border, site_border_width))
        end
	end

	# write the footerstring to close the svg file
	write(file, getSVGFooterString())
	# close the file
	close(file)

	# convert to pdf
    if export_pdf
	    run(`inkscape $(filename_output) --export-pdf $(filename_output[1:end-4]).pdf`)
    end
	# if file shall be opened
	if openfile
        if is_linux()
		    run(`ristretto $(filename_output)`)
        elseif is_windows()
            try
                if export_pdf
                    run(`explorer $(filename_output[1:end-4]).pdf`)
                else
                    run(`explorer $(filename_output)`)
                end                
            end
        elseif is_apple()
            run(`open $(filename_output)`)
        else
            println("wanted to show file but operating system not supported for showing file")
        end
	end

	# return the output filename
	return filename_output

end
export plotLattice2D


function plotLattice3D(
		lattice::Lattice;
		border_percentage=0.1,
		filename_output::String="AUTO",
		site_radius=25,
		site_border_width_percentage::Float64=0.2,
        site_labels="OFF",
		bond_thickness::Int64=8,
		visualize_periodic=false,
		colorcode_sites = Dict(0 => [255,255,255], 1 => [255,255,255]),
		colorcode_bonds = Dict("0" => [0,0,0], "1.0" => [0,0,0]),
        colorcode_bonds_automation::String = "OFF",
		bond_color_BG = [100,100,100],
		site_color_BG = [125,125,125],
		lattice_rotation::Array{Float64}=[0.0,0.0,0.0],
		camera_position_relative::Array{Float64}=[0.0,0.0,5.0],
		DOF::Float64 = 0.8,
		openfile=false,
        export_pdf=true
		)
	
    # define the filename_output if it is set to AUTO
    if filename_output=="AUTO"
        filename_output = "$(lattice.filename[1:end-4])_plot.svg"
        filename_output = replace(filename_output, FOLDER_LATTICES, "")
    end

	# load positions and connections
	positions = copy(lattice.positions)
	connections = copy(lattice.connections)

    # maybe overwrite dictonary
    if colorcode_bonds_automation == "GREY"
        # construct the list of interaction strengths
        cs_list = getConnectionStrengthList(lattice)
        # get the color code list
        cc_list = getGreySequence(size(cs_list, 1))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[string(cs_list[i])] = cc_list[i]
        end
    elseif colorcode_bonds_automation == "COLOR"
        # construct the list of interaction strengths
        cs_list = getConnectionStrengthList(lattice)
        # get the color code list
        cc_list = getColorSequence(size(cs_list, 1))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[string(cs_list[i])] = cc_list[i]
        end
    end

	# repair color dictonary
	colorcode_bonds["0"] = get(colorcode_bonds, "0", [0,0,0])
	colorcode_sites[0] = get(colorcode_sites, 0, [255,255,255])

	# shift the positions so they are centered in (0,0,0)
	position_center = [0,0,0]
	for pos in positions
		position_center = position_center .+ pos
	end
	position_center = position_center ./ size(positions,1)
	for index in 1:size(positions,1)
		positions[index] = positions[index] .- position_center
	end

	# rotate around all axis
	for (index,angle_degrees) in enumerate(lattice_rotation)
		# get the angle in radians
		angle = angle_degrees *pi/180.0
		# define the rotation matrix
		R = eye(3) * cos(angle)
		R[index, index] = 1
		R[(index+1)%3+1,(index)%3+1] = sin(angle)
		R[(index)%3+1,(index+1)%3+1] = -sin(angle)
		#println(R)
		# rotate the positions
		for index in 1:size(positions,1)
			positions[index] = R*positions[index]
		end
	end


	# vector absolute values
	function vecabs(r)
		return sqrt(r[1]*r[1] + r[2]*r[2] + r[3]*r[3])
	end

	# determine the position on screen
	# with reference to https://en.wikipedia.org/wiki/3D_projection

	# get the arrays of screen positions
	b_x = zeros(size(positions, 1))
	b_y = zeros(size(positions, 1))
	# the array of distance to camera
	distances_to_camera = zeros(size(positions,1))

	# determine the constants to be used

	# camera position in xyz coordinates
	max_x = maximum([abs(pos[1]) for pos in positions])
	max_y = maximum([abs(pos[2]) for pos in positions])
	max_z = maximum([abs(pos[3]) for pos in positions])
	c_x = max_x * camera_position_relative[1]
	c_y = max_y * camera_position_relative[2]
	c_z = max_z * camera_position_relative[3]
	c_vec = [c_x, c_y, c_z]

	# camera angles
	theta_x = 2*pi * (0)
	theta_y = 2*pi * (0)
	theta_z = 2*pi * (0)

	# viewer position
	e_x = 0
	e_y = 0
	e_z = 1

	# calculate the projection of each points
	for index in 1:length(b_x)
		# determine the vector x
		x = positions[index][1] - c_x
		y = positions[index][2] - c_y
		z = positions[index][3] - c_z
		# calculate the distance to the camera of all positions
		distances_to_camera[index] = vecabs(positions[index] .- c_vec)
		# determine the displacements
		d_x = cos(theta_y)*(sin(theta_z)*y + cos(theta_z)*x)  -  sin(theta_y)*(z)
		d_y = sin(theta_x)*(cos(theta_y)*z + sin(theta_y)*(sin(theta_z)*y + cos(theta_z)*x))  +  cos(theta_x)*(cos(theta_z)*y - sin(theta_z)*x)
		d_z = cos(theta_x)*(cos(theta_y)*z + sin(theta_y)*(sin(theta_z)*y + cos(theta_z)*x))  -  sin(theta_x)*(cos(theta_z)*y - sin(theta_z)*x)
		# determine the positions on screen
		b_x[index] = (e_z/d_z) * d_x - e_x
		b_y[index] = (e_z/d_z) * d_y - e_y
	end


	min_i, max_i 	= -1,1
	min_j, max_j 	= -1,1
	border 			= border_percentage * (maximum(b_x) + maximum(b_y) - minimum(b_x) - minimum(b_y))/4


	# die farbe der sites
	function getProperSiteColor(distance, site_index)
		# get the color of that particular site
		site_color_bare = get(colorcode_sites, site_index, colorcode_sites[0])
		site_color_back = (site_color_bare .* (1-0.3*DOF)) .+ (site_color_BG .* 0.3*DOF)
		# the desired color of the site
		site_color_desired = [255, 255, 255]
		for i in 1:3
			# get the color max and min values
			site_color_max = site_color_back[i]
			site_color_min = site_color_bare[i]
			# parameters of the fit
			a = -(site_color_max-site_color_min)/((1/min_distance) - (1/max_distance))
			b = site_color_max - a * 1/max_distance
			# get the best grey value
			grey = (a / distance) + b
			grey = Int(floor(grey))
			# set the value
			site_color_desired[i] = grey
		end
		# return the hex version of the proper site color
		return site_color_desired
	end

	# proper radius
	radius_min = site_radius
	radius_max = Int(floor(radius_min/2.0))
	function getProperSiteRadius(distance)
		# overwrite the radius
		radius_max = radius_min*min_distance/max_distance
		# parameters of the fit
		a = -(radius_max-radius_min)/((1/min_distance) - (1/max_distance))
		b = radius_max - a * 1/max_distance
		# get the best radius
		radius = Int(floor((a / distance) + b))
		return radius
	end
	function getProperSiteBorderRadius(distance)
		# overwrite the radius
		radius_max = radius_min*min_distance/max_distance
		# parameters of the fit
		a = -(radius_max-radius_min)/((1/min_distance) - (1/max_distance))
		b = radius_max - a * 1/max_distance
		# get the best radius
		radius = ((a / distance) + b)*site_border_width_percentage
		return "$(radius)px"
	end

	# proper connection color
	function getProperConnectionColor(distance, strength)
		# parse to string
		strength = string(strength)
		# get the color of that particular bond
		bond_color_bare = get(colorcode_bonds, strength, colorcode_bonds["0"])
		bond_color_back = (bond_color_bare .* (1-DOF)) .+ (bond_color_BG .* DOF)
		# the desired color of the site
		bond_color_desired = [255, 255, 255]
		for i in 1:3
			# get the color max and min values
			bond_color_max = bond_color_back[i]
			bond_color_min = bond_color_bare[i]
			# parameters of the fit
			a = -(bond_color_max-bond_color_min)/((1/min_distance) - (1/max_distance))
			b = bond_color_max - a * 1/max_distance
			# get the best grey value
			grey = (a / distance) + b
			grey = Int(floor(grey))
			# set the value
			bond_color_desired[i] = grey
		end
		# return the hex version of the proper site color
		return color_hex(bond_color_desired)
	end
	function getProperSiteBorderColor(distance)
		# get the color of that 0 bond
		bond_color_bare = get(colorcode_bonds, "0", colorcode_bonds["0"])
		bond_color_back = (bond_color_bare .* (1-0.1*DOF)) .+ (bond_color_BG .* 0.1*DOF)
		# the desired color of the site
		bond_color_desired = [255, 255, 255]
		for i in 1:3
			# get the color max and min values
			bond_color_max = bond_color_back[i]
			bond_color_min = bond_color_bare[i]
			# parameters of the fit
			a = -(bond_color_max-bond_color_min)/((1/min_distance) - (1/max_distance))
			b = bond_color_max - a * 1/max_distance
			# get the best grey value
			grey = (a / distance) + b
			grey = Int(floor(grey))
			# set the value
			bond_color_desired[i] = grey
		end
		# return the hex version of the proper site color
		return color_hex(bond_color_desired)
	end

	# proper width
	width_min = bond_thickness
	width_max = width_min * 0.6
	function getProperConnectionWidth(distance)
		# parameters of the fit
		a = -(width_max-width_min)/((1/min_distance) - (1/max_distance))
		b = width_max - a * 1/max_distance
		# get the best radius
		width = (a / distance) + b
		return width
	end


	# DEFINE SITES THAT ARE PLOTTED

	# sites to plot
	sites_to_plot = Array[]
	for index in 1:length(b_x)
        if site_labels == "LATTICE INDEX"
		    push!(sites_to_plot, [b_x[index], b_y[index], distances_to_camera[index], sum(positions[index].*c_vec), lattice.positions_indices[index], index])
        else 
            push!(sites_to_plot, [b_x[index], b_y[index], distances_to_camera[index], sum(positions[index].*c_vec), lattice.positions_indices[index]])
        end
	end

	# sort the sites s.t. the z values arrange in descending order
	function sortfunctionSites(site)
		return -site[4]
	end
	sort!(sites_to_plot, by=sortfunctionSites)


	# DEFINE CONNECTIONS THAT ARE PLOTTED 

	# connections to plot (sort out periodic ones)
	connections_to_plot = Array[]
	neutral_connection_wrap = (0,0,0)
	if length(lattice.connections[1][4]) == 2
		neutral_connection_wrap = (0,0)
	elseif length(lattice.connections[1][4]) == 1
		neutral_connection_wrap = 0
	end
	for c in connections
		if c[4] != neutral_connection_wrap && !visualize_periodic
			continue
		end
		# startpositions
		p1 = positions[Int(c[1])]
		p2 = positions[Int(c[2])]
		pm = (p1.+p2) * 0.5
		# calculate the camera distance
		pm_distance = vecabs(pm.-c_vec)
		# push the connection
		#push!(connections_to_plot, [c; pm_distance; pm[3]-c_z])
		push!(connections_to_plot, [c; pm_distance; min(sum(p1.*c_vec), sum(p2.*c_vec)) - 0.01])
	end
	# sort the connections s.t. the z values arrange in descending order
	function sortfunctionConnections(con)
		return -con[6] #max(con[5],con[6])
	end
	sort!(connections_to_plot, by=sortfunctionConnections)





	# overwrite positions from here on with projected positions
	for index in 1:length(b_x)
		positions[index] = [b_x[index],b_y[index]]
	end

	# generate the xvals and yvals only for size
	xvals = b_x
	yvals = b_y
	# define the width and height
	width_UC 	= (maximum(xvals) + border) - (minimum(xvals) - border)
	height_UC	= (maximum(yvals) + border) - (minimum(yvals) - border)
	min_x	= minimum(xvals) - border
	min_y	= minimum(yvals) - border
	if width_UC > height_UC
		width	= 1200
		conversion = width/width_UC
		height	= conversion*(height_UC)
	else
		height	= 1200
		conversion = height/height_UC
		width	= conversion*(width_UC)
	end

	# max und min distance
	max_distance = maximum(distances_to_camera)
	min_distance = minimum(distances_to_camera)

	# define the conversion functions for coordinates
	function X(x)
		return conversion*(x - min_x)
	end
	function Y(y)
		return + conversion*(y - min_y)
	end

	# open the file
	file = open(filename_output, "w")

	# write the headerstring
	write(file, getSVGHeaderString(round(Int64, width), round(Int64, height)))


	# write all strings
	i = 0
	while size(connections_to_plot,1) > 0 || size(sites_to_plot,1) > 0
		# if no more connections, only plot sites
		if size(connections_to_plot,1) == 0
			# pop the first site
			site_to_plot = pop!(sites_to_plot)
            site_color_basic = getProperSiteColor(site_to_plot[3],site_to_plot[5])
            site_color = color_hex(site_color_basic)
            label_color = color_hex([colorelment < 100 ? 255 : 0 for colorelment in site_color_basic])
            # write 
            if site_labels == "POSITION INDEX" || site_labels == "LATTICE INDEX"
                write(file, getSVGStringEllipseStroked("el$(i)", 
                    X(site_to_plot[1]), Y(site_to_plot[2]),
                    getProperSiteRadius(site_to_plot[3]), getProperSiteRadius(site_to_plot[3]),
                    site_color,
                    getProperSiteBorderColor(site_to_plot[3]), getProperSiteBorderRadius(site_to_plot[3]),
                    label=Int(site_to_plot[end]), labelcolor=label_color))
            else
                write(file, getSVGStringEllipseStroked("el$(i)", 
                    X(site_to_plot[1]), Y(site_to_plot[2]),
                    getProperSiteRadius(site_to_plot[3]), getProperSiteRadius(site_to_plot[3]),
                    site_color,
                    getProperSiteBorderColor(site_to_plot[3]), getProperSiteBorderRadius(site_to_plot[3])))
            end
			# increase index
			i = i+1
		elseif size(sites_to_plot,1) == 0
			# pop the first connection
			connection_to_plot = pop!(connections_to_plot)
			c = connection_to_plot
			# write
			if visualize_periodic || c[4] == neutral_connection_wrap
                write(file,
                    getSVGStringLine("path$(i)", 
                    [X(positions[Int(c[1])][1]), Y(positions[Int(c[1])][2])],
                    [X(positions[Int(c[2])][1]), Y(positions[Int(c[2])][2])],
                    getProperConnectionColor(c[5], c[3]),
                    getProperConnectionWidth(c[5]), dashed=(c[4] != neutral_connection_wrap)))	
            end
			# increase index
			i = i+1
		else
			# decide which one to pop
			if sites_to_plot[end][4] < connections_to_plot[end][6] #min(connections_to_plot[end][5],connections_to_plot[end][6]) - 0.0001
				# pop the first site
				site_to_plot = pop!(sites_to_plot)
                site_color_basic = getProperSiteColor(site_to_plot[3],site_to_plot[5])
                site_color = color_hex(site_color_basic)
                label_color = color_hex([colorelment < 100 ? 255 : 0 for colorelment in site_color_basic])
                # write 
                if site_labels == "POSITION INDEX" || site_labels == "LATTICE INDEX"
                    write(file, getSVGStringEllipseStroked("el$(i)", 
                        X(site_to_plot[1]), Y(site_to_plot[2]),
                        getProperSiteRadius(site_to_plot[3]), getProperSiteRadius(site_to_plot[3]),
                        site_color,
                        getProperSiteBorderColor(site_to_plot[3]), getProperSiteBorderRadius(site_to_plot[3]),
                        label=Int(site_to_plot[end]), labelcolor=label_color))
                else
                    write(file, getSVGStringEllipseStroked("el$(i)", 
                        X(site_to_plot[1]), Y(site_to_plot[2]),
                        getProperSiteRadius(site_to_plot[3]), getProperSiteRadius(site_to_plot[3]),
                        site_color,
                        getProperSiteBorderColor(site_to_plot[3]), getProperSiteBorderRadius(site_to_plot[3])))
                end
				# increase index
				i = i+1
			else
				# pop the first connection
				connection_to_plot = pop!(connections_to_plot)
				c = connection_to_plot
				# write
				if visualize_periodic || c[4] == neutral_connection_wrap
					write(file,
						getSVGStringLine("path$(i)", 
						[X(positions[Int(c[1])][1]), Y(positions[Int(c[1])][2])],
						[X(positions[Int(c[2])][1]), Y(positions[Int(c[2])][2])],
						getProperConnectionColor(c[5], c[3]),
						getProperConnectionWidth(c[5]), dashed=(c[4] != neutral_connection_wrap)))	
				end	
			end
			# increase index
			i = i+1
		end
	end




	

	# write the footerstring to close the svg file
	write(file, getSVGFooterString())
	# close the file
	close(file)

	# convert to pdf
    if export_pdf
	    run(`inkscape $(filename_output) --export-pdf $(filename_output[1:end-4]).pdf`)
    end
	
	# if file shall be opened	
	if openfile
        if is_linux()
		    run(`ristretto $(filename_output)`)
        elseif is_windows()
            try
                if export_pdf
                    run(`explorer $(filename_output[1:end-4]).pdf`)
                else
                    run(`explorer $(filename_output)`)
                end                
            end 
        elseif is_apple()
            run(`open $(filename_output)`)
        else
            println("wanted to show file but operating system not supported for showing file")
        end
	end

	# return the output filename
	return filename_output

end
export plotLattice3D


function plotLattice(
		lattice::Lattice;
		conversion = 160,
		border_percentage=0.1,
		filename_output::String="AUTO",
		site_radius=25,
		site_border_width_percentage::Float64=0.2,
        site_labels="OFF",
		bond_thickness::Int64=8,
		visualize_periodic=false,
		colorcode_sites = Dict(0 => [255,255,255], 1 => [255,255,255]),
		colorcode_bonds = Dict("0" => [0,0,0], "1.0" => [0,0,0]),
        colorcode_bonds_automation::String = "OFF",
		bond_color_BG = [100,100,100],
		site_color_BG = [125,125,125],
		lattice_rotation::Array{Float64}=[0.0,0.0,0.0],
		camera_position_relative::Array{Float64}=[0.0,0.0,5.0],
		DOF::Float64 = 0.8,
		openfile=false,
        export_pdf=true
		)

    # check which dimension the lattice has
    dimension = length(lattice.positions[1])

    # determine which function to call, depending on dimension
    if dimension == 2
        return plotLattice2D(
            lattice,
		    conversion=conversion,
            border_percentage=border_percentage,
            filename_output=filename_output,
            site_radius=site_radius,
            site_border_width_percentage=site_border_width_percentage,
            site_labels = site_labels,
            bond_thickness=bond_thickness,
            visualize_periodic=visualize_periodic,
            colorcode_sites = colorcode_sites,
            colorcode_bonds = colorcode_bonds,
            colorcode_bonds_automation = colorcode_bonds_automation,
            openfile=openfile,
            export_pdf=export_pdf)
    elseif dimension == 3
        return plotLattice3D(
            lattice,
            border_percentage=border_percentage,
            filename_output=filename_output,
            site_radius=site_radius,
            site_border_width_percentage=site_border_width_percentage,
            site_labels = site_labels,
            bond_thickness=bond_thickness,
            visualize_periodic=visualize_periodic,
            colorcode_sites = colorcode_sites,
            colorcode_bonds = colorcode_bonds,
            colorcode_bonds_automation = colorcode_bonds_automation,
            bond_color_BG = bond_color_BG,
            site_color_BG = site_color_BG,
            lattice_rotation=lattice_rotation,
            camera_position_relative=camera_position_relative,
            DOF=DOF,
            openfile=openfile,
            export_pdf=export_pdf)
    else
        println("plotting to SVG not implemented for dimension $(dimension) of the lattice")
        return
    end
end
export plotLattice






function plotPlaquettes2D(
		lattice::Lattice,
        plaquettes,
        plaquette_values;
		conversion = 160,
		border_percentage=0.1,
		filename_output::String="AUTO",
		site_radius=25,
		site_border_width_percentage::Float64=0.2,
        site_labels="OFF",
		bond_thickness::Int64=8,
		visualize_periodic=false,
		colorcode_sites = Dict(0 => [255,255,255], 1 => [255,255,255]),
		colorcode_bonds = Dict("0" => [0,0,0], "1.0" => [0,0,0]),
		colorcode_plaquettes = Dict(0.0 => [0,0,0], 1.0 => [255,0,0], -1.0 => [0,100,255]),
        colorcode_bonds_automation::String = "OFF",
		openfile=false,
        export_pdf=true
		)

    # define the filename_output if it is set to AUTO
    if filename_output=="AUTO"
        filename_output = "$(lattice.filename[1:end-4])_plq_plot.svg"
        filename_output = replace(filename_output, FOLDER_LATTICES, "")
    end

	# load positions and connections
	positions	= lattice.positions
	connections = lattice.connections

    # maybe overwrite dictonary
    if colorcode_bonds_automation == "GREY"
        # construct the list of interaction strengths
        cs_list = getConnectionStrengthList(lattice)
        # get the color code list
        cc_list = getGreySequence(size(cs_list, 1))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[string(cs_list[i])] = cc_list[i]
        end
    elseif colorcode_bonds_automation == "COLOR"
        # construct the list of interaction strengths
        cs_list = getConnectionStrengthList(lattice)
        # get the color code list
        cc_list = getColorSequence(size(cs_list, 1))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[string(cs_list[i])] = cc_list[i]
        end
    end

	# repair color dictonary
	colorcode_bonds["0"] = get(colorcode_bonds, "0", [0,0,0])
	colorcode_sites[0] = get(colorcode_sites, 0, [255,255,255])
    colorcode_plaquettes[0.0] = get(colorcode_plaquettes, 0.0, [0,0,0])

	# define styles for the different sites, i.e. strings that are saved into the svg strings
	site_r 				= site_radius
	site_border			= "#000000"
	site_border_width	= "$(site_border_width_percentage*site_radius)px"

	# sites to plot
	pos_x = zeros(size(positions,1))
	pos_y = zeros(size(positions,1))
	for i in 1:size(positions,1)
		pos_x[i] = positions[i][1]
		pos_y[i] = positions[i][2]
	end
	sites_to_plot = positions
	indices_to_plot = lattice.positions_indices
	xvals = pos_x
	yvals = pos_y
	border 			= border_percentage * (maximum(pos_x) + maximum(pos_y) - minimum(pos_x) - minimum(pos_y))/4

	# connections to plot (all)
	connections_to_plot = copy(connections)
    # the neutral connection wrap, i.e. which wrap identifies a non-periodic connection
	neutral_connection_wrap = (0,0)
	if size(lattice.lattice_vectors,1) == 1
		neutral_connection_wrap = (0)
	end
	

	# define the width and height of the canvas
	width_UC 	= (maximum(xvals) + site_radius/conversion + border) - (minimum(xvals) - site_radius/conversion - border)
	height_UC	= (maximum(yvals) + site_radius/conversion + border) - (minimum(yvals) - site_radius/conversion - border)
	min_x	= minimum(xvals) - site_radius/conversion - border
	min_y	= minimum(yvals) - site_radius/conversion - border
	width	= conversion*(width_UC)
	height	= conversion*(height_UC)

	# define the conversion functions for coordinates
	function X(x)
		return conversion*(x - min_x)
	end
	function Y(y)
		return + conversion*(y - min_y)
	end

	# open the SVG file
	file = open(filename_output, "w")

	# write the headerstring
	write(file, getSVGHeaderString(round(Int64, width), round(Int64, height)))

    # write all plaquettes
    for (i,p) in enumerate(plaquettes)
        # get the color
        color_plaquette = color_hex(get(colorcode_plaquettes, plaquette_values[i], colorcode_plaquettes[0.0]))
        # get the plaquette point list
        points = Array[]
        for index in p
            push!(points, [X(positions[index][1]), Y(positions[index][2])])
        end
        # plot plaquette
        write(file, getSVGStringPlaquette("plaq$(i)", points, color_plaquette))
    end
	# write all connections
	for (i,c) in enumerate(connections_to_plot)
		if (c[4] == neutral_connection_wrap)
			connection_color = color_hex(get(colorcode_bonds, string(c[3]), colorcode_bonds["0"]))
			write(file, getSVGStringLine("path$(i)", [X(positions[Int(c[1])][1]), Y(positions[Int(c[1])][2])], [X(positions[Int(c[2])][1]), Y(positions[Int(c[2])][2])], connection_color, bond_thickness))
		elseif visualize_periodic
			connection_color = color_hex(get(colorcode_bonds, string(c[3]), colorcode_bonds["0"]))
			write(file, getSVGStringLine("path$(i)", [X(positions[Int(c[1])][1]), Y(positions[Int(c[1])][2])], [X(positions[Int(c[2])][1]), Y(positions[Int(c[2])][2])], connection_color, bond_thickness, dashed=true))
        end
	end

	# write all sites
	for (i,s) in enumerate(sites_to_plot)
        site_color_basic = get(colorcode_sites, indices_to_plot[i], colorcode_sites[0])
		site_color = color_hex(site_color_basic)
        label_color = color_hex([colorelment < 100 ? 255 : 0 for colorelment in site_color_basic])
        if site_labels == "POSITION INDEX"
		    write(file, getSVGStringEllipseStroked("el$(i)", X(s[1]), Y(s[2]), site_r, site_r, site_color, site_border, site_border_width, label=indices_to_plot[i], labelcolor=label_color))
        elseif site_labels == "LATTICE INDEX"
		    write(file, getSVGStringEllipseStroked("el$(i)", X(s[1]), Y(s[2]), site_r, site_r, site_color, site_border, site_border_width, label=i, labelcolor=label_color))
        else
		    write(file, getSVGStringEllipseStroked("el$(i)", X(s[1]), Y(s[2]), site_r, site_r, site_color, site_border, site_border_width))
        end
	end

	# write the footerstring to close the svg file
	write(file, getSVGFooterString())
	# close the file
	close(file)

	# convert to pdf
    if export_pdf
	    run(`inkscape $(filename_output) --export-pdf $(filename_output[1:end-4]).pdf`)
    end
	# if file shall be opened
	if openfile
        if is_linux()
		    run(`ristretto $(filename_output)`)
        elseif is_windows()
            try
                if export_pdf
                    run(`explorer $(filename_output[1:end-4]).pdf`)
                else
                    run(`explorer $(filename_output)`)
                end                
            end
        elseif is_apple()
            run(`open $(filename_output)`)
        else
            println("wanted to show file but operating system not supported for showing file")
        end
	end

	# return the output filename
	return filename_output

end
export plotPlaquettes2D




# function to show the lattice
function showLattice(
		lattice::Lattice;
		conversion = 160,
		site_radius=25,
        site_labels="OFF",
		bond_thickness::Int64=8,
		visualize_periodic=false,
		colorcode_sites = Dict(0 => [255,255,255], 1 => [255,255,255]),
		colorcode_bonds = Dict("0" => [0,0,0], "1.0" => [0,0,0]),
        colorcode_bonds_automation::String = "OFF",
        background_color = [200,200,200],
        forcePyPlot=false
    )
    # if forcing PyPlot, do this and exit
    if forcePyPlot
        showLatticePyPlot(
            lattice,
            conversion = conversion,
            site_radius=site_radius,
            site_labels=site_labels,
            bond_thickness=bond_thickness,
            visualize_periodic=visualize_periodic,
            colorcode_sites = colorcode_sites,
            colorcode_bonds = colorcode_bonds,
            colorcode_bonds_automation=colorcode_bonds_automation,
            background_color = background_color
        )
        return
    end
    # CHECK IF MAYAVI IS PRESENT
    try
        @pyimport mayavi.mlab as mlab
        MAYAVI_AVAILABLE = true
        try
        showLatticeMayavi(
            lattice,
            conversion = conversion,
            site_radius=site_radius,
            site_labels=site_labels,
            bond_thickness=bond_thickness,
            visualize_periodic=visualize_periodic,
            colorcode_sites = colorcode_sites,
            colorcode_bonds = colorcode_bonds,
            colorcode_bonds_automation=colorcode_bonds_automation,
            background_color = background_color
        )
        catch x
        end
        return
    catch error
        if isa(error, PyCall.PyError)
            #println("PyError occured when importing:")
            #println(error)
            MAYAVI_AVAILABLE = false
            #println("Using PyPlot now")
            try
                showLatticePyPlot(
                    lattice,
                    conversion = conversion,
                    site_radius=site_radius,
                    site_labels=site_labels,
                    bond_thickness=bond_thickness,
                    visualize_periodic=visualize_periodic,
                    colorcode_sites = colorcode_sites,
                    colorcode_bonds = colorcode_bonds,
                    colorcode_bonds_automation=colorcode_bonds_automation,
                    background_color = background_color
                )
            catch error2
                println("Error occured when plotting:")
                println(error2)
            end
        else
            println("Error occured when loading:")
            println(error)
        end
    end

end
export showLattice


# show Lattice function to plot in pyplot
function showLatticePyPlot(
		lattice::Lattice;
		conversion = 160,
		site_radius=25,
        site_labels="OFF",
		bond_thickness::Int64=8,
		visualize_periodic=false,
		colorcode_sites = Dict(0 => [255,255,255], 1 => [255,255,255]),
		colorcode_bonds = Dict("0" => [0,0,0], "1.0" => [0,0,0]),
        colorcode_bonds_automation::String = "OFF",
        background_color = [200,200,200]
    )

    # maybe overwrite dictonary
    if colorcode_bonds_automation == "GREY"
        # construct the list of interaction strengths
        cs_list = getConnectionStrengthList(lattice)
        # get the color code list
        cc_list = getGreySequence(size(cs_list, 1))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[string(cs_list[i])] = cc_list[i]
        end
    elseif colorcode_bonds_automation == "COLOR"
        # construct the list of interaction strengths
        cs_list = getConnectionStrengthList(lattice)
        # get the color code list
        cc_list = getColorSequence(size(cs_list, 1))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[string(cs_list[i])] = cc_list[i]
        end
    end

	# repair color dictonary
	colorcode_bonds["0"] = get(colorcode_bonds, "0", [0,0,0])
	colorcode_sites[0] = get(colorcode_sites, 0, [255,255,255])

    # PLOTTING WITH PYPLOT
    # create a figure
    fig = figure(figsize=(10,10), facecolor=background_color ./ 255.0)
    # add a 3d projection subplot
    ax = fig[:add_subplot](111, projection="3d", axisbg=background_color ./ 255.0)
    # define a function to plot a sphere for every site
    function plotSite(site,radius,color; detail=10)
        x = site[1]
        y = site[2]
        if length(site) == 3
            z = site[3]
        else
            z = 0.0
        end
        # Make data
        n = detail
        u = linspace(0,2*pi,n);
        v = linspace(0,pi,n);
        px = cos(u) * sin(v)';
        py = sin(u) * sin(v)';
        pz = ones(n) * cos(v)';
        px = (px.*radius) .+ x
        py = (py.*radius) .+ y
        pz = (pz.*radius) .+ z
        # Plot the surface
        surf(px,py,pz, rstride=1, cstride=1, linewidth=0, antialiased=true, color=color)
    end
    # define a function to plot a bond as a tube
    function plotBond(from, to, radius, color; detail_ring=8, detail_length=3)
        # check to be 3d
        if length(from) == 2
            from = [from[1], from[2], 0.0]
        end        
        if length(to) == 2
            to = [to[1], to[2], 0.0]
        end        
        # https://de.mathworks.com/matlabcentral/answers/151235-tube-plot-with-x-y-coordinates-and-radius?requestedDomain=www.mathworks.com
        # plot bond
        vec_u = to .- from;
        t = nullspace(vec_u')';
        vec_v = t[1,:]
        vec_w = t[2,:]
        # normalize
        vec_v = radius.* vec_v ./ norm(vec_v)
        vec_w = radius.* vec_w ./ norm(vec_w)
        #println("v=$(vec_v)    w=$(vec_w)")
        #[S,T] = linspace(0,1,detail)*linspace(0,2*pi,detail)';
        S = linspace(0,1,detail_length)
        T = linspace(0,2*pi,detail_ring)
        px = zeros(detail_length, detail_ring)
        py = zeros(detail_length, detail_ring)
        pz = zeros(detail_length, detail_ring)
        for i in 1:detail_length
        for j in 1:detail_ring        
            #P = repmat(from',m*n,1) .+ S.*vec_u .+ cos(T).*vec_v .+ sin(T).*vec_w;
            px[i,j] = from[1] + S[i]*vec_u[1] + cos(T[j])*vec_v[1] + sin(T[j])*vec_w[1];
            py[i,j] = from[2] + S[i]*vec_u[2] + cos(T[j])*vec_v[2] + sin(T[j])*vec_w[2];
            pz[i,j] = from[3] + S[i]*vec_u[3] + cos(T[j])*vec_v[3] + sin(T[j])*vec_w[3];
        end
        end
        surf(px,py,pz, rstride=1, cstride=1, linewidth=0, antialiased=true, color=color)
    end
   
    # plot all bonds
    for c in lattice.connections
        # skip if wrong directioin
        if c[1] < c[2] || (!visualize_periodic && sum([abs(el) for el in c[4]]) != 0)
            continue
        end
        # get data
        from = lattice.positions[Int(c[1])] .* conversion
        to = lattice.positions[Int(c[2])] .* conversion
        color = get(colorcode_bonds, string(c[3]), colorcode_bonds["0"]) ./ 255.0
        # plot
        plotBond(from, to, bond_thickness, color)
    end
    # plot all sites
    for index in 1:size(lattice.positions,1)
        p = lattice.positions[index]
        if length(p) == 2
            p = [p[1], p[2], 0.0]
        end
        p = p.*conversion
        plotSite(p, site_radius, get(colorcode_sites, lattice.positions_indices[index], colorcode_sites[0])./255.0)
    end
    # plot all site labels
    for index in 1:size(lattice.positions,1)
        p = lattice.positions[index]
        if length(p) == 2
            p = [p[1], p[2], 0.0]
        end
        p = p.*conversion
        # add label
        if site_labels == "LATTICE INDEX"
            ax[:text](p[1]+site_radius, p[2]+site_radius, p[3]+site_radius,  "$(index)", zorder=1, size=site_radius, color="k") 
        elseif site_labels == "POSITION INDEX"
            ax[:text](p[1]+site_radius, p[2]+site_radius, p[3]+site_radius,  "$(lattice.positions_indices[index])", zorder=1, size=site_radius, color="k") 
        end
    end


    # turn off everthing but the plot
    ax[:set_axis_off]()
    # get the current limits    
    x_limits = ax[:get_xlim3d]()
    y_limits = ax[:get_ylim3d]()
    z_limits = ax[:get_zlim3d]()
    # build the current ranges
    x_range = abs(x_limits[2] - x_limits[1])
    x_middle = (x_limits[2] + x_limits[1]) / 2.0
    y_range = abs(y_limits[2] - y_limits[1])
    y_middle = (y_limits[2] + y_limits[1]) / 2.0
    z_range = abs(z_limits[2] - z_limits[1])
    z_middle = (z_limits[2] + z_limits[1]) / 2.0
    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*maximum([x_range, y_range, z_range])
    ax[:set_xlim3d](x_middle - plot_radius, x_middle + plot_radius)
    ax[:set_ylim3d](y_middle - plot_radius, y_middle + plot_radius)
    ax[:set_zlim3d](z_middle - plot_radius, z_middle + plot_radius)
    # tighten layout
    tight_layout()
    # show
    show()
end
# show Lattice function to plot in mayavi
function showLatticeMayavi(
		lattice::Lattice;
		conversion = 160,
		site_radius=25,
        site_labels="OFF",
		bond_thickness::Int64=8,
		visualize_periodic=false,
		colorcode_sites = Dict(0 => [255,255,255], 1 => [255,255,255]),
		colorcode_bonds = Dict("0" => [0,0,0], "1.0" => [0,0,0]),
        colorcode_bonds_automation::String = "OFF",
        background_color = [200,200,200]
    )

    # CHECK IF MAYAVI IS PRESENT
    @pyimport mayavi.mlab as mlab

    # maybe overwrite dictonary
    if colorcode_bonds_automation == "GREY"
        # construct the list of interaction strengths
        cs_list = getConnectionStrengthList(lattice)
        # get the color code list
        cc_list = getGreySequence(size(cs_list, 1))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[string(cs_list[i])] = cc_list[i]
        end
    elseif colorcode_bonds_automation == "COLOR"
        # construct the list of interaction strengths
        cs_list = getConnectionStrengthList(lattice)
        # get the color code list
        cc_list = getColorSequence(size(cs_list, 1))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[string(cs_list[i])] = cc_list[i]
        end
    end

	# repair color dictonary
	colorcode_bonds["0"] = get(colorcode_bonds, "0", [0,0,0])
	colorcode_sites[0] = get(colorcode_sites, 0, [255,255,255])

    
    # PLOTTING WITH MAYAVI
    # create a figure
    mlab.figure(bgcolor=(background_color[1], background_color[2], background_color[3]))
    # define a function to plot a sphere for every site
    function plotSite(site,radius,color; detail=10)
        x = site[1]
        y = site[2]
        if length(site) == 3
            z = site[3]
        else
            z = 0.0
        end
        mlab.points3d([x],[y],[z],scale_factor=radius, color=(color[1], color[2], color[3]))
    end
    # define a function to plot a bond as a tube
    function plotBond(from, to, radius, color; detail_ring=8, detail_length=3)
        # check to be 3d
        if length(from) == 2
            from = [from[1], from[2], 0.0]
        end        
        if length(to) == 2
            to = [to[1], to[2], 0.0]
        end        
        mlab.plot3d([from[1], to[1]], [from[2], to[2]], [from[3], to[3]], color=(color[1], color[2], color[3]), tube_radius=radius)
    end
   
    # plot all bonds
    for c in lattice.connections
        # skip if wrong directioin
        if c[1] < c[2] || (!visualize_periodic && sum([abs(el) for el in c[4]]) != 0)
            continue
        end
        # get data
        from = lattice.positions[Int(c[1])] .* conversion
        to = lattice.positions[Int(c[2])] .* conversion
        color = get(colorcode_bonds, string(c[3]), colorcode_bonds["0"]) ./ 255.0
        # plot
        plotBond(from, to, bond_thickness, color)
    end
    # plot all sites
    for index in 1:size(lattice.positions,1)
        p = lattice.positions[index]
        if length(p) == 2
            p = [p[1], p[2], 0.0]
        end
        p = p.*conversion
        plotSite(p, site_radius, get(colorcode_sites, lattice.positions_indices[index], colorcode_sites[0])./255.0)
    end
    # plot all site labels
    for index in 1:size(lattice.positions,1)
        p = lattice.positions[index]
        if length(p) == 2
            p = [p[1], p[2], 0.0]
        end
        p = p.*conversion
        # add label
        if site_labels == "LATTICE INDEX"
            mlab.text3d(p[1]+site_radius, p[2]+site_radius, p[3]+site_radius, "$(index)", scale=site_radius, color=(0,0,0)) 
        elseif site_labels == "POSITION INDEX"
            mlab.text3d(p[1]+site_radius, p[2]+site_radius, p[3]+site_radius, "$(lattice.positions_indices[index])", scale=site_radius, color=(0,0,0)) 
        end
    end


    # show plot
    mlab.show()


    # return the figure
    return fig
end













#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#
#   CALCULATIONS OF BAND STRUCTURE AND SPIN GROUNDSTATES
#
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------






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









#-----------------------------------------------------------------------------------------------------------------------------
#
#   METHODS FOR CALCULATING THE BAND STRUCTURE OF A MATRIX IN K SPACE
#
#-----------------------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------------------
#
#   BAND STRUCTURE ALONG A PATH IN K SPACE
#
#   path has to be of following format
#   Array[
#       ["name1", [coordinates1]],
#       ["name2", [coordinates2]],
#       ...
#   ]
#   path does not close! to close, insert first point again
#
#   Parameters (necessary):
#   - lattice: The lattice object of which the band structure should be calculated
#   - path: The path along which the band structure should be calculated (in format given above)
#
#   Parameters (optional):
#   - reduceLattice: If the lattice should be reduced to a 1x1(x1) lattice of the original unitcell (i.e. for purposes of Luttinger Tisza, etc.)
#   - percentages: Either array of floats giving the individual percentages of the segments or the String "EQUAL" for equal length in the plot
#   - resolution: How many k-points to calculate in total
#   - enforce_hermitian: If the matrix should be enforced to be hermitian
#   - limits_energy: The y-axis (energy axis) limits of the plot
#   - plot_title: A title for the whole plot ("AUTO" for automatic title)
#   - plot_color: The color of bands
#   - figsize: the figure size as given to PyPlot
#
#-----------------------------------------------------------------------------------------------------------------------------
# BAND STRUCTURE ALONG PATH
function calculateBandStructureAlongPath(
        interaction_matrix::Array{Complex,2},
        path;
        percentages="EQUAL",
        resolution=1000,
        enforce_hermitian=false,
        limits_energy="AUTO",
        plot_title="",
        plot_color="b",
        figsize=(6,4),
        showPlot=true
            )
    
    # normalize percentages
    if percentages == "EQUAL"
        percentages = ones(size(path,1)-1)
    end
    percentages = percentages ./ sum(percentages)
    # build up the segment list
    segments = Array[]
    for i in 1:size(path,1)-1
        segment = [i, i+1, percentages[i]]
        push!(segments, segment)
    end
    # if LT is checked
    #if check_LT
    #    LT_k = Array[]
    #end
    # segment data, i.e. the bandstructure over the segments
    segments_data = Array[]
    resolution_actual = 0
    hlines = []
    # iterate over all segments
    for segment in segments
        # get the grid in between two points
        resolution_local = Int(floor(segment[3]*resolution))
        multipliers = linspace(0, 1, resolution_local)
        resolution_actual += resolution_local
        push!(hlines, resolution_actual+0)
        #println(segment)
        k1 = convert(Array{Float64,1}, path[Int(segment[1])][2:end])
        k2 = convert(Array{Float64,1}, path[Int(segment[2])][2:end])
        #println(k1)
        #println(k2)
        # insert bands
        bands = Array[]
        for b in 1:size(interaction_matrix,1)
            push!(bands, zeros(resolution_local))
        end
        # calculate all energies
        for i in 1:resolution_local
            # get the current k
            k = k2 .* multipliers[i] .+ k1 .* (1-multipliers[i])
            # if LT is checked, push current k
            #if check_LT
            #    push!(LT_k, k)
            #end
            # get the interaction matrix for this k
            matrix = copy(interaction_matrix)
            # diagonalize the matrix
            eigenvalues = eigvals(matrix)
            # save all the eigenvalues
            for b in 1:size(bands, 1)
                if imag(eigenvalues[b]) > 0
                    if imag(eigenvalues[b]) > 1e-15
                        println(imag(eigenvalues[b]))
                        println(matrix)
                        bands[b][i] = eigenvalues[b]
                    else
                        bands[b][i] = real(eigenvalues[b])
                    end                    
                else
                    bands[b][i] = eigenvalues[b]
                end
            end
        end
        # push the obtained back structure into the data array
        push!(segments_data, bands)
    end
    # generate the complete band structure
    bandstructure = Array[zeros(resolution_actual) for b in segments_data[1]]
    index = 1
    for i in 1:size(segments_data,1)
        segment = segments[i]
        data = segments_data[i]
        for b in 1:size(bandstructure,1)
            bandstructure[b][index:hlines[i]] = data[b]
        end
        index = hlines[i]+1
    end
    # if LT is checked, give the results
    #if check_LT
    #    LT_v = checkLuttingerTisza(lattice, LT_k, only_GS=false)
    #    println("$(100.0*sum(LT_v)/length(LT_v)) % of all eigenvalues are valid in LT")
    #end
    # plot the eigenvalues
    rc("font", family="serif")
    fig = figure(figsize=figsize)
    if plot_title == "AUTO"
        title("energy spectrum along path of interaction matrix")
    elseif plot_title == ""
        # do nothing title related
    else
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
    figurename = "interaction_matrix"
    figurename1 = "$(FOLDER_SPECTRA)bandstructure_path_$(figurename[1:end-4]).pdf"
    figurename2 = "$(FOLDER_SPECTRA)bandstructure_path_$(figurename[1:end-4]).png"
    buildFolderSpectra()
    savefig(figurename1)
    savefig(figurename2)
    if showPlot
        show()
        #print("Continue? ")
        #readline()
    end
    return fig
end
function calculateBandStructureAlongPath(
        lattice::Lattice,
        path;
        reduceLattice=true,
        percentages="EQUAL",
        resolution=1000,
        enforce_hermitian=false,
        limits_energy="AUTO",
        plot_title="",
        plot_color="b",
        figsize=(6,4),
        showPlot=true,
        majorana=false
            )
    
    # check if to reduce the lattice
    if reduceLattice && lattice.unitcell.filename != UNITCELL_DUMMY_FILENAME
        lattice = getLatticePeriodic(lattice.unitcell, ones(Int64, size(lattice.unitcell.lattice_vectors,1)), save=false)
    end
    # normalize percentages
    if percentages == "EQUAL"
        percentages = ones(size(path,1)-1)
    end
    percentages = percentages ./ sum(percentages)
    # build up the segment list
    segments = Array[]
    for i in 1:size(path,1)-1
        segment = [i, i+1, percentages[i]]
        push!(segments, segment)
    end
    # if LT is checked
    #if check_LT
    #    LT_k = Array[]
    #end
    # segment data, i.e. the bandstructure over the segments
    segments_data = Array[]
    resolution_actual = 0
    hlines = []
    # iterate over all segments
    for segment in segments
        # get the grid in between two points
        resolution_local = Int(floor(segment[3]*resolution))
        multipliers = linspace(0, 1, resolution_local)
        resolution_actual += resolution_local
        push!(hlines, resolution_actual+0)
        #println(segment)
        k1 = convert(Array{Float64,1}, path[Int(segment[1])][2:end])
        k2 = convert(Array{Float64,1}, path[Int(segment[2])][2:end])
        #println(k1)
        #println(k2)
        # insert bands
        bands = Array[]
        for b in 1:size(lattice.positions,1)
            push!(bands, zeros(resolution_local))
        end
        # calculate all energies
        for i in 1:resolution_local
            # get the current k
            k = k2 .* multipliers[i] .+ k1 .* (1-multipliers[i])
            # if LT is checked, push current k
            #if check_LT
            #    push!(LT_k, k)
            #end
            # get the interaction matrix for this k
            matrix = getInteractionMatrixKSpace(lattice, k, enforce_hermitian=enforce_hermitian, majorana=majorana)
            # diagonalize the matrix
            eigenvalues = eigvals(matrix)
            # save all the eigenvalues
            for b in 1:size(bands, 1)
                if imag(eigenvalues[b]) > 0
                    if imag(eigenvalues[b]) > 1e-15
                        println(imag(eigenvalues[b]))
                        println(matrix)
                        bands[b][i] = eigenvalues[b]
                    else
                        bands[b][i] = real(eigenvalues[b])
                    end                    
                else
                    bands[b][i] = eigenvalues[b]
                end
            end
        end
        # push the obtained back structure into the data array
        push!(segments_data, bands)
    end
    # generate the complete band structure
    bandstructure = Array[zeros(resolution_actual) for b in segments_data[1]]
    index = 1
    for i in 1:size(segments_data,1)
        segment = segments[i]
        data = segments_data[i]
        for b in 1:size(bandstructure,1)
            bandstructure[b][index:hlines[i]] = data[b]
        end
        index = hlines[i]+1
    end
    # if LT is checked, give the results
    #if check_LT
    #    LT_v = checkLuttingerTisza(lattice, LT_k, only_GS=false)
    #    println("$(100.0*sum(LT_v)/length(LT_v)) % of all eigenvalues are valid in LT")
    #end
    # plot the eigenvalues
    rc("font", family="serif")
    fig = figure(figsize=figsize)
    if plot_title == "AUTO"
        if majorana
            title("majorana energy spectrum along path of lattice \"$(lattice.filename)\"")
        else
            title("energy spectrum along path of lattice \"$(lattice.filename)\"")
        end
    elseif plot_title == ""
        # do nothing title related
    else
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
    if majorana
        figurename1 = "$(FOLDER_SPECTRA)majorana_bandstructure_path_$(figurename[1:end-4]).pdf"
        figurename2 = "$(FOLDER_SPECTRA)majorana_bandstructure_path_$(figurename[1:end-4]).png"
    else
        figurename1 = "$(FOLDER_SPECTRA)bandstructure_path_$(figurename[1:end-4]).pdf"
        figurename2 = "$(FOLDER_SPECTRA)bandstructure_path_$(figurename[1:end-4]).png"
    end
    buildFolderSpectra()
    savefig(figurename1)
    savefig(figurename2)
    if showPlot
        show()
        #print("Continue? ")
        #readline()
    end
    return fig
end
function calculateBandStructureAlongPath(
        unitcell::Unitcell,
        path;
        percentages="EQUAL",
        resolution=1000,
        enforce_hermitian=false,
        limits_energy="AUTO",
        plot_title="",
        plot_color="b",
        figsize=(6,4),
        showPlot=true,
        majorana=false
            )
    
    # make a lattice from the unitcell
    lattice = getLatticePeriodic(unitcell, ones(Int64, size(unitcell.lattice_vectors,1)), save=false)
    lattice.filename = replace(identity(unitcell.filename), FOLDER_UNITCELLS, FOLDER_LATTICES)
    # push to lattice based method and return the result
    return calculateBandStructureAlongPath(
        lattice,
        path;
        reduceLattice=false,
        percentages=percentages,
        resolution=resolution,
        enforce_hermitian=enforce_hermitian,
        limits_energy=limits_energy,
        plot_title=plot_title,
        plot_color=plot_color,
        figsize=figsize,
        showPlot=showPlot,
        majorana=majorana
            )
end
export calculateBandStructureAlongPath

# BAND STRUCTURE OF STRIP (1D periodic)
function calculateBandStructureOfStrip(
        unitcell::Unitcell,
        periodic_direction::Int64,
        finite_direction_N::Int64;
        percentages="EQUAL",
        resolution=1000,
        enforce_hermitian=false,
        limits_energy="AUTO",
        plot_title="",
        plot_color="b",
        figsize=(6,4),
        showPlot=true,
        majorana=false
    )

    # obtain the lattice
    dimensions = ones(Int64, size(unitcell.lattice_vectors, 1)) .* finite_direction_N
    dimensions[periodic_direction] = -1
    lattice = getLattice(unitcell, dimensions)

    # obtain the path
    last_vector = lattice.lattice_vectors[1]
    last_vector = last_vector
    a = sqrt(sum(last_vector.*last_vector))
    gamma_zero  = last_vector .* 0.0
    gamma_minus = last_vector .* (-2*pi) / (a*a)
    gamma_plus  = last_vector .* ( 2*pi) / (a*a)
    K_minus     = last_vector .* (-1*pi) / (a*a)
    K_plus      = last_vector .* ( 1*pi) / (a*a)

    path = Array[
        ["Gamma (-1)";  gamma_minus],
        ["-pi";         K_minus],
        ["Gamma (0)";   gamma_zero],
        ["pi";          K_plus],
        ["Gamma (+1)";  gamma_plus]
    ]

    # calculate the band structure
    calculateBandStructureAlongPath(
        lattice,
        path;
        reduceLattice=false,
        percentages=percentages,
        resolution=resolution,
        enforce_hermitian=enforce_hermitian,
        limits_energy=limits_energy,
        plot_title=plot_title,
        plot_color=plot_color,
        figsize=figsize,
        showPlot=showPlot,
        majorana=majorana
    )
end
export calculateBandStructureOfStrip


# FULL 2D Bandstructure
function calculateBandStructure2D(
        unitcell::Unitcell,
        kx, ky;
        enforce_hermitian=false,
        limits_energy="AUTO",
        plot_title="",
        plot_color="b",
        figsize=(6,4),
        showPlot=true,
        majorana=false
            )
    
    
    # insert bands
    bands = Array[]
    kx_vals = zeros(length(kx), length(ky))
    ky_vals = zeros(length(kx), length(ky))
    for b in 1:size(unitcell.basis,1)
        push!(bands, zeros(length(kx), length(ky)))
    end
    # calculate all energies
    for i in 1:length(kx)
    for j in 1:length(ky)
        # get the current k
        k = [kx[i], ky[j]]
        kx_vals[i,j] = k[1]
        ky_vals[i,j] = k[2]
        # get the interaction matrix for this k
        matrix = getInteractionMatrixKSpace(unitcell, k, enforce_hermitian=enforce_hermitian, majorana=majorana)
        # diagonalize the matrix
        eigenvalues = eigvals(matrix)
        # save all the eigenvalues
        for b in 1:size(bands, 1)
            if imag(eigenvalues[b]) > 0
                if imag(eigenvalues[b]) > 1e-15
                    println(imag(eigenvalues[b]))
                    println(matrix)
                    bands[b][i,j] = eigenvalues[b]
                else
                    bands[b][i,j] = real(eigenvalues[b])
                end                    
            else
                bands[b][i,j] = eigenvalues[b]
            end
        end
    end
    end
    # generate the complete band structure
    bandstructure = bands
    # if LT is checked, give the results
    #if check_LT
    #    LT_v = checkLuttingerTisza(lattice, LT_k, only_GS=false)
    #    println("$(100.0*sum(LT_v)/length(LT_v)) % of all eigenvalues are valid in LT")
    #end
    # plot the eigenvalues
    rc("font", family="serif")
    fig = figure(figsize=figsize)
    ax = fig[:add_subplot](111, projection="3d")
    if plot_title == "AUTO"
        if majorana
            title("majorana energy spectrum along path of unitcell \"$(unitcell.filename)\"")
        else
            title("energy spectrum along path of unitcell \"$(unitcell.filename)\"")
        end
    elseif plot_title == ""
        # do nothing title related
    else
        title(plot_title)
    end
    xlabel("momentum kx")
    ylabel("momentum ky")
    zlabel("energy")
    for b in bandstructure
        ax[:plot_surface](kx_vals,ky_vals,b, rstride=1, cstride=1, cmap="coolwarm", linewidth=0)
    end
    axx = ax[:get_xaxis]()
    # check if specific boundaries are desired
    if !(limits_energy == "AUTO")
        ylim(limits_energy[1], limits_energy[2])
    end
    # tighten the layout
    tight_layout()
    # save the plot
    figurename = split(unitcell.filename, FOLDER_SPECTRA[end])[end]
    if majorana
        figurename1 = "$(FOLDER_SPECTRA)majorana_bandstructure_path_$(figurename[1:end-4]).pdf"
        figurename2 = "$(FOLDER_SPECTRA)majorana_bandstructure_path_$(figurename[1:end-4]).png"
    else
        figurename1 = "$(FOLDER_SPECTRA)bandstructure_path_$(figurename[1:end-4]).pdf"
        figurename2 = "$(FOLDER_SPECTRA)bandstructure_path_$(figurename[1:end-4]).png"
    end
    buildFolderSpectra()
    savefig(figurename1)
    savefig(figurename2)
    if showPlot
        show()
        #print("Continue? ")
        #readline()
    end
    return fig   
end
export calculateBandStructure2D






# SOME DEFAULT PATHS
DEFAULT_PATH_FCC = Array[
    ["gamma"; [0,0,0]],
    ["X"; [2*pi, 0, 0]],
    ["W"; [2*pi, pi, 0]],
    ["L"; [pi, pi, pi]],
    ["gamma"; [0,0,0]],
    ["K"; [3*pi/2, 3*pi/2, 0]],
    ["X"; [2*pi, 0, 0]]
]
export DEFAULT_PATH_FCC

DEFAULT_PATH_TRIANGULAR = Array[
    ["gamma"; [0,0]],
    ["K"; [2*pi/sqrt(3.0), 2*pi/3]],
    ["M"; [2*pi/sqrt(3.0), 0]],
    ["gamma"; [0,0]]
]
export DEFAULT_PATH_TRIANGULAR

DEFAULT_PATH_SQUAREOCTAGON_2 = Array[
    ["gamma"; [0,0]],
    ["K"; [2,  0].*(pi / (1.0 + 1.0/sqrt(2.0)))],
    ["M"; [1, -1].*(pi / (1.0 + 1.0/sqrt(2.0)))],
    ["gamma"; [0,0]]
]
export DEFAULT_PATH_SQUAREOCTAGON_2


DEFAULT_PATH_SQUARE_LONG = Array[
    ["M";     [pi,   0.0]],
    ["Gamma"; [0.0,  0.0]],
    ["K'";    [pi,   -pi]],
    ["M";     [pi,   0.0]],
    ["M'";    [0.0,   pi]],
    ["K";     [pi,    pi]],
    ["Gamma"; [0.0,  0.0]]
]
export DEFAULT_PATH_SQUARE_LONG

DEFAULT_PATH_SQUARE_SHORT = Array[
    ["Gamma"; [0.0,  0.0]],
    ["M";     [pi,   0.0]],
    ["K";     [pi,    pi]],
    ["Gamma"; [0.0,  0.0]]
]
export DEFAULT_PATH_SQUARE_SHORT

DEFAULT_PATH_SQUARE = DEFAULT_PATH_SQUARE_SHORT
export DEFAULT_PATH_SQUARE

























#-----------------------------------------------------------------------------------------------------------------------------
#
#   METHODS FOR JOB SCHEDULING
#
#-----------------------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------------------
#
#   FUNCTION TO CREATE A JOB DIRECTORY FOR SCHEDULING
#
#   USAGE: Give a jobdirectory name and a main program to run as well as a default input file and the function creates
#   a job directory containing all input files as well as two scripts to run the program in parallel or single core to
#   calculate all input files.
#   Usage should be in the following steps:
#   1) provide an input file template of the form: <parameter name>\t<default value>
#   2) provide a main program which can be called by a command similar to "PROGRAM_NAME INPUTFILE.txt"
#   3) call this function to create job directory and all input files
#   4) Enter job directory and run one of the two julia execution scripts
#       - Single core execution one after another: "julia script_single_execution.jl"
#       - MPI based multi-core execution: "mpirun -x PATH --host <host1, host2, host3,...> --map-by node -np <number of processes> nice -n19 julia script_MPI_execution.jl"
#   5) collect output files
#
#
#   Parameters (necessary):
#   - directory_path:       The name of the directory to be created
#   - input_file_template:  The input file template from which parameters are extracted
#   - main_program_file:    The file which contains the main program to be run
#
#   Parameters (optional):
#   - parameters:           A dictonary of parameter values to be used
#   - paramters_use_ARGS:   Bool to determine if the ARGS string is searched for parameters
#   - input_file_seperator: The seperator used for input file seperation
#   - main_program_execution_command:   The command to launch the program PROGRAM with a given input file INPUTFILE. Default: julia PROGRAM INPUTFILE
#   - files_to_copy:        A list of file that are copied to the destination folder as well
#
#-----------------------------------------------------------------------------------------------------------------------------
function createJobDirectory(
    directory_path::String,
    input_file_template_name::String,
    main_program_file::String
    ;
    parameters::Dict=Dict(),
    paramters_use_ARGS::Bool=true,
    input_file_seperator::String=":\t",
    main_program_execution_command::String="julia PROGRAM INPUTFILE",
    files_to_copy::Array{String,1}=String[]
        )
    

    # -----------------
    # START INFORMATION
    # -----------------

    # correct the directory input
    if !endswith(directory_path, "/")
        directory_path = "$(directory_path)/"
    end
    # Start: Print the directory path
    println("Job directory \"$(directory_path)\" will be created")
    # Print the input template filename
    println("Input template is located at \"$(input_file_template_name)\"\n")

    


    # ---------------------------------------------------
    # OBTAIN THE INPUT FILE TEMPLATE AND EXTRACT KEYWORDS
    # ---------------------------------------------------

    # Open the template file and read all lines
    input_file_template         = open(input_file_template_name)
    input_file_template_lines   = readlines(input_file_template)
    close(input_file_template)

    # Create a new list for all parameters
    parameters_RAW              = Array[]
    parameters_passed_values    = Array[]

    # push in all relevant lines
    for line in input_file_template_lines
        line = strip(line)
        if startswith(line, "#") || line==""
            continue
        end
        push!(parameters_RAW, split(line, "\t"))
    end

    # print all parameters that are found
    println("Found $(size(parameters_RAW, 1)) parameters in template")
    for p in parameters_RAW
        if length(p) == 1
            println("- \"$(p[1])\" without default value")
        elseif length(p) == 2
            println("- \"$(p[1])\" with default value $(p[2])")
        end
    end
    println("")




    # ---------------------------------------------------
    # OBTAIN THE CORRECT VALUES TO BE USED
    # ---------------------------------------------------

    # Iterate over all accessible parameters
    for p in parameters_RAW

        # list of values to use
        values_to_use   = []

        # OPTION 1, check in Command line arguments
        if paramters_use_ARGS
            # look through all command line arguments
            for arg in ARGS
                # check if it fulfils the structure
                if startswith(arg, "$(p[1])=")
                    # get the string representation
                    arg_value = replace(arg, "$(p[1])=", "")
                    # push to the list of values to use
                    push!(values_to_use, arg_value)
                    # break the loop
                    break
                end
            end
        end
        
        # OPTION 2, check in dictonary
        if length(values_to_use)==0 && p[1] in keys(parameters)
            # get the string representation
            dict_value = parameters(p[1])
            # push to the list of values to use
            push!(values_to_use, dict_value)
        end

        # OPTION 3, check default values of template
        if length(values_to_use)==0 && length(p) >= 2
            # get the string representation
            default_value = p[2]
            # push to the list of values to use
            push!(values_to_use, default_value)
        end


        # PARSING OF VALUES

        # create a list with parsed entries
        values_to_use_parsed = []
        
        # iterate over all values to be used
        for v in values_to_use
            # determine which type of value
            if typeof(v) != String && typeof(v) != SubString{String}
                # already parsed, push into array
                push!(values_to_use_parsed, v)
            elseif startswith(v, "[") && endswith(v, "]")
                # list of format [el1, el2, el3, ...]
                v = split(strip(v[2:end-1]), ",")
                for i in 1:length(v)
                    v[i] = strip(v[i])
                    push!(values_to_use_parsed, parse(Float64,v[i]))
                end
            elseif startswith(v, "linspace(") && endswith(v, ")")
                # list via linspace(start, end, steps)
                v = split(v[10:end-1], ",")
                for i in 1:length(v)
                    v[i] = strip(v[i])
                end
                v = linspace(parse(Float64,v[1]),parse(Float64,v[2]),parse(Float64,v[3]))
                for i in 1:length(v)
                    push!(values_to_use_parsed, v[i])
                end
            else
                # it is probably either a number or a string, try parsing to float
                if isnull(tryparse(Float64,v))
                    # not a float, just push in list as it is
                    push!(values_to_use_parsed, v)
                else
                    # push parsed into list
                    push!(values_to_use_parsed, parse(Float64,v))
                end
            end
        end


        # append the list of paramters to the global list
        push!(parameters_passed_values, values_to_use_parsed)

    end

    # print all values that are found
    println("Found values for parameters:")
    for (p_index,p) in enumerate(parameters_RAW)
        println("- \"$(p[1])\", $(length(parameters_passed_values[p_index])) values: $(parameters_passed_values[p_index])")
    end
    println("")






    # -----------
    # CREATE JOBS
    # -----------

    # create a new list of dictonaries which are the jobs
    joblist = []
    # add a dummy job without information
    push!(joblist, Dict())

    # iterate over all parameters
    for p in 1:length(parameters_RAW)

        # find out if single value or multiple values
        if length(parameters_passed_values[p]) == 1
            # only one value, all jobs get this value
            for j in joblist
                j["$(parameters_RAW[p][1])"] = parameters_passed_values[p][1]
            end
        elseif length(parameters_passed_values[p]) >= 1
            # more than one value, every job gets copied for every value
            joblist_backup = joblist
            joblist = []
            for j in joblist_backup
                for i in 1:length(parameters_passed_values[p])
                    j_tmp = copy(j)
                    j_tmp["$(parameters_RAW[p][1])"] = parameters_passed_values[p][i]
                    push!(joblist, j_tmp)
                end
            end
        else
            # PROBLEM: NO VALUE FOR PARAMETER, exlude from job
            println("DANGER: Excluding parameter $(parameters_RAW[p][1]) from jobs because it has no value")
        end
        
    end

    # print the success
    println("Job building is done, $(length(joblist)) jobs are constructed\n")

    # check if only one job without parameters
    if length(joblist) == 1 && length(keys(joblist[1])) == 0
        println("Some mistake occured, only job has no parameters. Aborting now...\n")
        return
    end





    # --------------------
    # CREATE THE DIRECTORY
    # --------------------

    # check if the directory exists
    if isdir(directory_path)
        # Directory exists
        println("Directory \"$(directory_path)\" already exists")
        print("Delete and reconstruct? [y/n]: ")
        answer = strip(readline()[end-1:end])
        while !in(answer, ["y", "n"])
            print("Delete and reconstruct? [y/n]: ")
            answer = strip(readline()[end-1:end])
        end
        # depending on answer
        if answer == "y"
            # ask for recursive deletion
            print("Delete recursively? [y/n]: ")
            answer = strip(readline()[end-1:end])
            while !in(answer, ["y", "n"])
                print("Delete and reconstruct? [y/n]: ")
                answer = strip(readline()[end-1:end])
            end
            if answer == "y"
                # delete the complete directory
                rm(directory_path, recursive=true)
            else
                # delete the complete directory
                rm(directory_path, recursive=false)
            end
            # print
            println("Deleted already existing directory $(directory_path)")
            # create the directory
            mkdir(directory_path)
            println("Recreated directory $(directory_path)\n")

        else
            # abort the program
            println("Not deleting directory, aborting program now...")
            return
        end
    else
        # Directory does not exist, create it
        mkdir(directory_path)
        println("Created directory $(directory_path)\n")
    end




    # ----------------------
    # CREATE THE INPUT FILES
    # ----------------------

    # iterate over all jobs
    for i in 1:length(joblist)
        # dump job i to file input_i.txt
        jobfilename = "input_$(i).txt"
        f = open("$(directory_path)$(jobfilename)", "w")
        # write entries of the dictonary into the file
        for k in keys(joblist[i])
            write(f, "$(k)$(input_file_seperator)$(joblist[i][k])\n")
        end
        # close the file
        close(f)
        # set the correct input filename in the job
        joblist[i]["input_filename"] = jobfilename
    end
    # print the success
    println("Created all $(length(joblist)) jobfiles\n")




    # ------------------------------
    # COPY / CREATE THE SCRIPT FILES
    # ------------------------------    

    # Copy main program file
    cp(main_program_file, "$(directory_path)$(main_program_file)")
    println("Copied main program file \"$(main_program_file)\"")

    # Create single execution julia script
    script_single_exe = open("$(directory_path)script_single_execution.jl", "w")
    # write all lines
    for i in 1:length(joblist)
        # compose a line
        line = "run(`$(main_program_execution_command)`)\n"
        # change the program part
        line = replace(line, "PROGRAM", main_program_file)
        # change the input part
        line = replace(line, "INPUTFILE", joblist[i]["input_filename"])
        # set the correct input filename in the job
        write(script_single_exe, line)
    end   
    # close the script
    close(script_single_exe)
    println("Created single execution script \"script_single_execution.jl\"")


    # Create single execution julia script
    script_MPI_exe = open("$(directory_path)script_MPI_execution.jl", "w")


    # write all lines, starting with imports
    write(script_MPI_exe, "# use with:\n")
    write(script_MPI_exe, "# mpirun -x PATH --host <host1,host2,...> --map-by node -np <number of processes> nice -n19 julia script_MPI_execution.jl\n")
    write(script_MPI_exe, "using MPI\n")
    write(script_MPI_exe, "MPI.Init()\n")
    # get the number of the MPI process
    write(script_MPI_exe, "commworld = MPI.COMM_WORLD\n")
    write(script_MPI_exe, "commrank  = MPI.Comm_rank(commworld)\n")
    write(script_MPI_exe, "commsize  = MPI.Comm_size(commworld)\n")
    write(script_MPI_exe, "max_i     = $(length(joblist))\n")
    write(script_MPI_exe, "init_i    = commrank+1\n")
    write(script_MPI_exe, "for i in init_i:commsize:max_i\n")
    # compose main line    
    line = "\trun(pipeline(`$(main_program_execution_command)`, stdout=\"bash_stdout_INPUTFILE\", stderr=\"bash_stderr_INPUTFILE\"))\n\tprintln(\"Done input file \\\"INPUTFILE\\\"\")\n"
    # change the program part
    line = replace(line, "PROGRAM", main_program_file)
    # change the input part
    line = replace(line, "INPUTFILE", "input_\$(i).txt")
    # set the correct input filename in the job
    write(script_MPI_exe, line)
    write(script_MPI_exe, "end\n")
    # Finalize MPI script
    write(script_MPI_exe, "MPI.Finalize()\n")
    # close the script
    close(script_MPI_exe)
    println("Created MPI execution script \"script_MPI_execution.jl\"")

    # Copy additional files
    for fn in files_to_copy
        cp(fn, "$(directory_path)$(fn)")
        println("Copied additional file \"$(fn)\"")
    end

end
export createJobDirectory





# MODULE END
end