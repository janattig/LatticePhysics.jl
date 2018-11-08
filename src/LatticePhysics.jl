################################################################################
#
#   The unified Module for all different sorts of lattices
#   and lattice based calculations.
#
#   Module structure mostly included in subfiles see below
#
################################################################################


# Start of module
# module LatticePhysics

# TODO LIST
# --> TODO TYPES ALWAYS FIRST ARGUMENTS IN FUNCTIONS



################################################################################
#
#   module LatticePhysics_Base
#
#   --> ABSTRACT TYPE DEFINITION
#           - AbstractSite{L,D}
#           - AbstractBond{L,N}
#           - AbstractUnitcell{S,B}
#           - AbstractLattice{S,B,U}
#
#   --> INTERFACING ABSTRACT TYPES
#           - AbstractSite{L,D}
#               - positions
#               - label
#           - AbstractBond{L,N}
#               - connecting identifiers
#               - label
#               - wrap
#           - AbstractUnitcell{S,B}
#               - Bravais lattice vectors
#               - sites
#               - bonds
#           - AbstractLattice{S,B,U}
#               - Unitcell
#               - Bravais lattice vectors
#               - sites
#               - bonds
#
#   --> NAIVE STRUCT DEFINITIONS
#           - Site{L,D}
#           - Bond{L,N}
#           - Unitcell{S,B}
#           - Lattice{S,B,U}
#
#   --> IMPLEMENTING JULIA.BASE FUNCTIONS
#           - show
#           - print
#
################################################################################

# include the relevant subfile
include("LatticePhysics_Base.jl")



################################################################################
#
#   module LatticePhysics_IO
#   -> LatticePhysics_Base
#   -> JLD2
#
#   --> IO FOR UNITCELLS / LATTICES
#           - save to xml ??
#           - JLD ??
#
################################################################################



################################################################################
#
#   module LatticePhysics_ConnectionRepresentation
#   -> LatticePhysics_Base
#   -> LinearAlgebra
#
#   --> FUNCTIONS FOR IDENTIFYING NEIGHBOURS
#
#   --> STATIC LOOKUP TABLES
#           - connectivity information
#           - lists of identifiers for every site
#           - list of booleans
#
#   --> INTERACTION MATRICES
#           - construction based on functions and dictonaries
#           - real space / momentum space?
#
################################################################################



################################################################################
#
#   module LatticePhysics_UnitcellDefinitions
#   -> LatticePhysics_Base
#   -> LinearAlgebra
#
#   --> PRE-IMPLEMENTED UNITCELLS 2D
#
#   --> PRE-IMPLEMENTED UNITCELLS 3D
#
################################################################################

# include the relevant subfile
include("LatticePhysics_UnitcellDefinitions.jl")



################################################################################
#
#   module LatticePhysics_UnitcellConstruction
#   -> LatticePhysics_Base
#   -> LinearAlgebra
#
#   --> CONSTRUCTION OF UNITCELLS FROM SITES
#
#   --> CONSTRUCTION HELPER
#           - print the generating code for a Unitcell
#
################################################################################

# include the relevant subfile
include("LatticePhysics_UnitcellConstruction.jl")



################################################################################
#
#   module LatticePhysics_LatticeConstruction
#   -> LatticePhysics_Base
#   -> LinearAlgebra
#
#   --> NAIVE CONSTRUCTION OF LATTICES
#           - periodic pattern of unitcells
#           - semi-periodic pattern of unitcells
#           - open pattern of unitcells
#           - by bond distance to an origin site
#           - in a shape around an origin site
#
#   --> ADDING NNN TO LATTICE
#
#   --> LATTICE CONSTRUCTION ALGORITHM FROM SPIN SPIRALS PAPER
#
################################################################################

# include the relevant subfile
include("LatticePhysics_LatticeConstruction.jl")



################################################################################
#
#   module LatticePhysics_LatticeModification
#   -> LatticePhysics_Base
#   -> LinearAlgebra
#
#   --> MODIFICATION OF LATTICES COMPONENTS
#       - Adding / Removing Connections
#       - Adding / Removing sites
#
#   --> RELABELING BONDS / SITES
#       - relabeling bonds
#       - relabeling sites
#
#   --> SPATIAL MODIFICATIONS
#       - Scaling  in space
#       - Rotating in space
#       - Shifting in space
#       - Lattice vector optimization (relabeling unitcells)
#
################################################################################



################################################################################
#
#   module LatticePhysics_Plaquettes
#   -> LatticePhysics_Base
#   -> LinearAlgebra
#
#   --> TYPE DEFINITION OF PLAQUETTE TYPE
#
#   --> PLAQUETTE OBTAINING FUNCTIONS
#           - get plaquettes of a site
#           - get plaquettes of a unitcell / lattice
#
#   --> PLAQUETTE INFORMATION FUNCTIONS
#           - obtain / print plaquette statistics
#
################################################################################



################################################################################
#
#   module LatticePhysics_ReciprocalSpace
#   -> LatticePhysics_Base
#   -> LinearAlgebra
#
#   --> TYPE DEFINITIONS
#           - Path
#           - BrillouinZone
#
#   --> PATHS IN MOMENTUM SPACE
#           - modification of paths
#           - pre-implemented paths
#           - paths based on Unitcells
#
#   --> BRILLOUIN ZONES
#           - type definition
#           - pre-implemented BZs
#           - construction of default BZs (2D)
#           - construction of default BZs (3D)
#
################################################################################



################################################################################
#
#   module LatticePhysics_ReciprocalSpace_Plotting
#   -> LatticePhysics_Base
#   -> LatticePhysics_ReciprocalSpace
#   -> LinearAlgebra
#
#   --> PLOTTING OF PATHS
#           - 2D paths
#           - 3D paths
#
#   --> PLOTTING OF BRILLOUIN ZONES
#           - 2D Brillouin zones
#           - 3D Brillouin zones
#
################################################################################



################################################################################
#
#   module LatticePhysics_SVGBase
#
#   --> SVG STRINGS FOR FILES
#           - svg header
#           - svg footer
#
#   --> SVG STRINGS FOR GEOMETRIC OBJECTS
#           - ellipsoids
#           - rectangles
#           - paths
#           - text
#
################################################################################



################################################################################
#
#   module LatticePhysics_Plotting_SVG
#   -> LatticePhysics_Base
#   -> LatticePhysics_SVGBase
#   -> LinearAlgebra
#
#   --> PLOTTING OF LATTICES TO SVG FILES
#           - 2D lattices
#           - 3D lattices
#
################################################################################



################################################################################
#
#   module LatticePhysics_PlaquettePlotting_SVG
#   -> LatticePhysics_Base
#   -> LatticePhysics_Plaquettes
#   -> LatticePhysics_SVGBase
#   -> LatticePhysics_Plotting_SVG
#   -> LinearAlgebra
#
#   --> PLOTTING OF PLAQUETTES AND LATTICES TO SVG FILES
#           - 2D plaquettes and lattices
#
################################################################################



################################################################################
#
#   module LatticePhysics_Plotting_Blender
#   -> LatticePhysics_Base
#   -> LinearAlgebra
#
#   --> PLOTTING OF LATTICES TO TEXT FILES TO BE READ BY BLENDER
#           - 2D lattices
#           - 3D lattices
#
################################################################################



################################################################################
#
#   module LatticePhysics_Plotting_PyPlot
#   -> PyPlot
#   -> LatticePhysics_Base
#   -> LinearAlgebra
#
#   --> PLOTTING OF LATTICES WITH PYPLOT
#           - 2D lattices
#           - 3D lattices
#
################################################################################



################################################################################
#
#   module LatticePhysics_FreeFermion_Calculations
#   -> LatticePhysics_Base
#   -> LatticePhysics_ReciprocalSpace
#   -> LinearAlgebra
#
#   --> TYPE DEFINITIONS
#           - FermionBandstructure
#           - FermiSurface
#
#   --> IO OF TYPES
#
#   --> CALCULATION OF BANDSTRUCTURES AND EIGENVALUES
#           - in real space (pure eigenvalues)
#           - in momentum space (along cardinal directions)
#           - in momentum space (along reciprocal lattice directions)
#           - in momentum space (along path)
#
#   --> CALCULATION OF FERMI SURFACES
#           - calculation for 2D unitcells
#           - calculation for 3D unitcells
#
################################################################################



################################################################################
#
#   module LatticePhysics_FreeFermion_Plotting
#   -> PyPlot
#   -> LatticePhysics_Base
#   -> LatticePhysics_ReciprocalSpace
#   -> LatticePhysics_FreeFermion_Calculations
#   -> LinearAlgebra
#
#   --> PLOTTING OF BANDSTRUCTURES AND EIGENVALUES
#           - in real space (pure eigenvalues)
#           - in momentum space (along cardinal directions)
#           - in momentum space (along reciprocal lattice directions)
#           - in momentum space (along path)
#
#   --> PLOTTING OF FERMI SURFACES
#           - for 1D fermi surfaces in 2D
#           - for 2D fermi surfaces in 3D
#
################################################################################



################################################################################
#
#   module LatticePhysics_LuttingerTisza_Calculations
#   -> LatticePhysics_Base
#   -> LatticePhysics_ReciprocalSpace
#   -> LinearAlgebra
#
#   --> TYPE DEFINITIONS
#           - LTBandstructure
#           - LTGroundstateManifold
#
#   --> IO OF TYPES
#
#   --> SPIN INTERACTION MATRICES
#
#   --> CALCULATION OF LT CONSTRAINTS
#
#   --> CALCULATION OF LT BANDSTRUCTURES AND EIGENVALUES
#           - in real space (pure eigenvalues)
#           - in momentum space (along cardinal directions)
#           - in momentum space (along reciprocal lattice directions)
#           - in momentum space (along path)
#
#   --> CALCULATION OF LT GROUNDSTATE MANIFOLDS
#           - calculation for 2D unitcells
#           - calculation for 3D unitcells
#
################################################################################



################################################################################
#
#   module LatticePhysics_LuttingerTisza_Plotting
#   -> PyPlot
#   -> LatticePhysics_Base
#   -> LatticePhysics_ReciprocalSpace
#   -> LatticePhysics_LuttingerTisza_Calculations
#   -> LinearAlgebra
#
#   --> PLOTTING OF LT BANDSTRUCTURES AND EIGENVALUES
#           - in real space (pure eigenvalues)
#           - in momentum space (along cardinal directions)
#           - in momentum space (along reciprocal lattice directions)
#           - in momentum space (along path)
#
#   --> PLOTTING OF LT GROUNDSTATE MANIFOLDS
#           - plotting of manifolds in 2D
#           - plotting of manifolds in 3D
#
################################################################################






# End of module
# end
