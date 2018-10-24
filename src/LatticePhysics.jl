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
#
#   --> IO FOR UNITCELLS / LATTICES
#           - save to xml ??
#           - JLD ??
#
################################################################################
################################################################################
#
#   module LatticePhysics_ConnectionRepresentation
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
#
#   --> PRE-IMPLEMENTED UNITCELLS 2D
#
#   --> PRE-IMPLEMENTED UNITCELLS 3D
#
################################################################################
################################################################################
#
#   module LatticePhysics_UnitcellConstruction
#
#   --> CONSTRUCTION OF UNITCELLS FROM SITES
#
#   --> CONSTRUCTION HELPER
#
################################################################################
################################################################################
#
#   module LatticePhysics_LatticeConstruction
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
#   6) PLAQUETTE RELATED FUNCTIONS
#       - obtaining plaquettes of lattice
#       - printing plaquette statistics
#
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
