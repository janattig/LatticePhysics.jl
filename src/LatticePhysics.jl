################################################################################
#
#   The unified Module for all different sorts of lattices
#   and lattice based calculations.
#
#   Module structure (mostly included in subfiles see below)
#
################################################################################
################################################################################
#
#   module LatticePhysics_Base
#
#   --> ABSTRACT TYPE DEFINITION
#           - AbstractSite{D,L}
#           - AbstractBond{N,L}
#           - AbstractUnitcell{D,N,L}
#           - AbstractLattice{D,N,L}
#
#   --> INTERFACING ABSTRACT TYPES
#           - AbstractSite{D,L}
#               - positions
#               - label
#           - AbstractBond{N,L}
#               - connecting identifiers
#               - label
#               - wrap
#           - AbstractUnitcell{D,N,L}
#               - Bravais lattice vectors
#               - sites
#               - bonds
#           - AbstractLattice{D,N,L}
#               - Unitcell
#               - Bravais lattice vectors
#               - sites
#               - bonds
#
#   --> STRUCT DEFINITIONS
#           - Site{D,L}
#           - Bond{N,L}
#           - Unitcell{D,N,L}
#           - Lattice{D,N,L}
#
#   --> IMPLEMENTING JULIA.BASE FUNCTIONS
#           - show
#           - print
#
################################################################################
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
#           - TODO create paths based on unitcell
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






#
#   8) PATHS (IN MOMENTUM SPACE)
#       - type definition
#       - modification of paths
#       - pre-implemented paths
#       - TODO create paths based on unitcell
#
#   9) BRILLOUIN ZONES
#       - type definition
#       - pre-implemented BZs
#       - construction of default BZs (2D)
#       - construction of default BZs (3D)
#       - plotting of BZ
#
#  10) SVG PLOTTING
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



# End of module
end
