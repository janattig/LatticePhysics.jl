################################################################################
#
#   The unified Module for all different sorts of lattices
#   and lattice based calculations.
#
#   Module structure (mostly included in subfiles):
#
#   1) TYPE DEFINITION
#       - Definition of UNITCELL TYPE
#       - Definition of LATTICE TYPE
#       - saving / loading types
#       - information on types
#       - connectivity of tpyes
#
#   2) PRE-IMPLEMENTED UNITCELLS (various 2D and 3D stuff)
#
#   3) CONSTRUCTION OF UNITCELLS FROM SITES / CONSTRUCTION HELPER
#
#   4) CONSTRUCTION OF LATTICES
#       - Building Lattices as (periodic/open/semiperiodic) repeated patterns of Unitcells
#       - Building Lattices by bond distance to an origin site
#       - Building Lattices in a shape around an origin site
#
#   5) MODIFICATION OF LATTICES
#       - Adding / Removing Connections
#       - Adding / Removing sites
#       - Interaction Strengths
#
#   6) PLAQUETTE RELATED FUNCTIONS
#       - obtaining plaquettes of lattice
#       - printing plaquette statistics
#
#   7) INTERACTION MATRICES (REAL AND MOMENTUM SPACE)
#
#   8) PATHS (IN MOMENTUM SPACE)
#       - type definition
#       - modification of paths
#       - pre-implemented paths
#
#   9) BRILLOUIN ZONES
#       - type definition
#       - pre-implemented BZs
#       - plotting of BZ
#
#   10) SVG PLOTTING
#       - 2D lattices
#       - TODO 3D lattices
#       - 2D plaquettes
#
#  11) BLENDER PLOTTING
#
#  12) BANDSTRUCTURES (CALCULATION AND PLOTTING)
#
#  13) FERMI SURFACES (CALCULATION AND PLOTTING)
#       - 2D unitcells
#       - 3D unitcells
#
#  14) LUTTINGER TISZA (CALCULATION AND PLOTTING)
#       - LTBandstructure type
#       - Spin interaction matrices
#       - calculation and plotting of LT Bandstructures (along path in k space)
#       - calculation and plotting of LT Groundstates (in k space)
#
################################################################################

# Start of module
module LatticePhysics








################################################################################
#
#   IMPORTED LIBRARIES THAT ARE USED IN LATTICEPHYSICS.JL
#
################################################################################

# JLD is used for all save / load applications
using JLD

# PyPlot is used for all plotting of band structures (plotting of lattices is done via SVG file creation)
using PyCall
using PyPlot

# Optim is used for minimizing the band structure to find the ground state energy of the system
using Optim






################################################################################
#
#   DEFAULT FILENAMES FOR FOLDERS AND OBJECTS
#
################################################################################


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






################################################################################
#
#   FUNCTIONS TO BUILD DEFAULT FOLDERS
#
################################################################################
function buildFolderUnitcells()
    mkpath(FOLDER_UNITCELLS)
end
export buildFolderUnitcells
function buildFolderLattices()
    mkpath(FOLDER_LATTICES)
end
export buildFolderLattices
function buildFolderSpectra()
    mkpath(FOLDER_SPECTRA)
end
export buildFolderSpectra


function buildFolders()
    buildFolderUnitcells()
    buildFolderLattices()
    buildFolderSpectra()
end
export buildFolders










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
#       - Testing a Unitcell object
#       - Print Information (printInfo(..))
#       - obtain the connection strengths of an object
#
#   5) OBTAINING CONNECTIVITY INFORMATION
#       - getConnectivityList(..)
#       - getConnectionList(..)
#
################################################################################

# included in subfile (1153 lines)
include("LatticePhysics_type_definitions.jl")






################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT UNITCELLS in 2D and 3D
#
#   STRUCTURE OF THE FILE
#
#   2D UNITCELLS
#   - SQUARE / RECTANGLE
#       - getUnitcellSquare
#       - getUnitcellExtendedSquare
#       - getUnitcellCheckerboard
#       - getUnitcellShastrySutherland
#       - getUnitcellAdvancedSquare
#       - getUnitcellSquareOctagon
#       - getUnitcellBCC2D
#       - getUnitcellFullyConnectedSquare
#   - TRIANGULAR
#       - getUnitcellTriangular
#       - getUnitcellHoneycomb
#       - getUnitcellKagome
#       - getUnitcellKagomeMinus
#       - getUnitcellHoneycombXXX
#
#   3D UNITCELLS
#   - CUBIC / FCC
#       - getUnitcellDiamond
#       - getUnitcellFCC
#       - getUnitcellBCC
#       - getUnitcellPyrochlore
#   - (X,3)y FAMILY
#       - getUnitcell_8_3_a
#       - getUnitcell_8_3_b
#       - getUnitcell_8_3_c
#       - getUnitcell_8_3_n
#       - getUnitcell_9_3_a
#       - getUnitcell_10_3_a / getUnitcellHyperoctagon
#       - getUnitcell_10_3_b / getUnitcellHyperhoneycomb
#       - getUnitcell_10_3_c
#       - getUnitcell_10_3_d
#
################################################################################

# included in subfile (4235 lines)
include("LatticePhysics_unitcell_implementations.jl")




################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT UNITCELL CONSTRUCTION RELATED STUFF
#
#   STRUCTURE OF THE FILE
#
#   1) UNITCELL FROM SITE COLLECTION (automatically determine the connections)
#       - 2D (not exported)
#       - 3D (not exported)
#       - independent of the dimension
#
#   2) UNITCELL GENERATING CODE
#
################################################################################

# included in subfile (554 lines)
include("LatticePhysics_unitcell_construction.jl")







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

# included in subfile (2079 lines)
include("LatticePhysics_lattice_construction.jl")



################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT LATTICE MODIFICATION FUNCTIONS
#   (MOST OF THEM ALSO WORK FOR UNITCELLS)
#
#   STRUCTURE OF THE FILE
#
#   1) CONNECTIONS ADDING REMOVING
#       - Add a new connection
#       - remove connections (based on indices)
#       - remove connections (based on strength)
#
#   2) SITES ADDING REMOVING
#       - add a new site
#       - remove a site (based on index)
#
#   3) CONNECTION STRENGTH MODIFICATION
#       - all connections
#       - mapping of strengths
#       - evaluate strengths
#       - optimize connections
#
################################################################################

# included in subfile (748 lines)
include("LatticePhysics_lattice_modification.jl")








################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT PLAQUETTE RELATED FUNCTIONS
#
#   STRUCTURE OF THE FILE
#
#   1) FINDING PLAQUETTES IN LATTICE OBJECTS
#   - of a site
#   - of complete lattice
#
#   2) PRINTING PLAQUETTE STATISTICS
#   - of a site
#   - of complete lattice
#
#
#   BUG: Plaquettes dont implement periodic boundaries so far
#
################################################################################

# included in subfile (333 lines)
include("LatticePhysics_plaquette_operations.jl")



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

# included in subfile (246 lines)
include("LatticePhysics_interaction_matrices.jl")





################################################################################
#
#   METHODS FOR CONSTRUCTION OF PATHS INSIDE THE BZ OF A UNITCELL
#
#   STRUCTURE OF THE FILE
#
#   1) TYPE PATH
#       - type definition
#       - printInfo function
#       - path String function (e.g. X-G-M-X)
#
#   2) CONSTRUCTION FUNCTIONS PATH
#       - function to add points
#       - TODO function to remove points
#       - function to scale the total resolution
#       - function to set the total resolution
#
#   3) DEFAULT PATHS
#       - getDefaultPathTriangular
#       - getDefaultPathSquare
#       - getDefaultPathSquareOctagon
#       - getDefaultPathFCC
#
#   4) TODO FUNCTION TO CREATE A DEFAULT PATH BASED ON A UNITCELL OBJECT
#
#   5) FUNCTION TO PLOT A PATH (in both 2D and 3D)
#
################################################################################

# included in subfile (799 lines)
include("LatticePhysics_paths.jl")




################################################################################
#
#   METHODS FOR CONSTRUCTION AND PLOTTING OF BRILLOUIN ZONES (BZ) OF A UNITCELL
#
#   STRUCTURE OF THE FILE
#
#   1) TYPE BRILLOUINZONE
#       - type definition
#       - TODO printInfo function
#
#   2) DEFAULT BZ
#       - TODO getDefaultBZSquare
#       - getDefaultBZFCC
#
#   3) TODO FUNCTION TO CREATE A DEFAULT BZ BASED ON A UNITCELL OBJECT
#
#   4) PLOTTING OF A BZ
#
################################################################################

# included in subfile (511 lines)
include("LatticePhysics_brillouin_zones.jl")








################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT LATTICE PLOTTING FUNCTIONS FOR
#   PLOTTING TO SVG IMAGES
#
#   STRUCTURE OF THE FILE
#
#   1) HELPER SVG FUNCTIONS (NOTE: not exported)
#       - Header / Footer String
#       - Ellipse (stroked) String
#       - Line String
#       - Plaquette String
#       - Color conversion
#       - Color collections / sequences
#
#   2) PLOTTING LATTICES
#       - plot in 2D
#       - TODO plot in 3D
#       - TODO plot independent of dimension
#
#   3) PLOTTING PLAQUETTES (only 2D)
#
################################################################################

# included in subfile (2094 lines)
include("LatticePhysics_SVG_plotting.jl")



################################################################################
#
#   METHODS FOR DUMPING A LATTICE TO BLENDER INPUT FILE
#   (input file can be used in Blender with the provided AddOn)
#
#   STRUCTURE OF THE FILE
#
#   1) SAVING LATTICE DATA TO FILE
#
################################################################################

# included in subfile (450 lines)
include("LatticePhysics_blender.jl")




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

# included in subfile (470 lines)
include("LatticePhysics_bandstructures.jl")






################################################################################
#
#   METHODS FOR CONSTRUCTION OF FERMI SURFACES (2D & 3D)
#
#   STRUCTURE OF THE FILE
#
#   1) CALCULATION OF FERMI SURFACE
#
#   2) PLOTTING OF FERMI SURFACE
#       - plotting from points
#       - plotting from unitcell
#
################################################################################

# included in subfile (773 lines)
include("LatticePhysics_fermi_surfaces.jl")





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
#      (- LT constraint and deviation functions (NOT EXPORTED) )
#       - spin interaction matrices
#       - calculation of band structures
#
#   3) PLOTTING OF LT BAND STRUCTURES
#       - plotting of Bandstructure objects
#       - plotting of bandstructures of unitcells along paths
#
#   4) CALCULATION OF LT GROUND STATES (k space manifold)
#
#   5) PLOTTING OF LT GROUND STATES (k space manifold)
#       - plotting from points
#       - plotting from unitcell
#
################################################################################

# included in subfile (1849 lines)
include("LatticePhysics_luttinger_tisza.jl")









# MODULE END
# total lines: 581 + 1153 + 4235 + 554 + 2079 + 748 + 333 + 246 + 799 + 511 + 2094 + 450 + 470 + 773 + 1849
# = 16875 lines
end
