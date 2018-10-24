################################################################################
#
#   The unified Module for all different sorts of lattices
#   and lattice based calculations.
#
#   Module structure (mostly included in subfiles):
#
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
#   --> INTERACING JULIA.BASE FUNCTIONS
#           - show
#           - print
#
################################################################################
#
#   module LatticePhysics_IO
#
#   --> IO FOR UNITCELLS / LATTICES
#           - save to xml ??
#           - JLD ??
#
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


#   --> 2) CONSTRUCTION OF LATTICES
#           - Building Lattices as (periodic/open/semiperiodic) repeated patterns of Unitcells
#           - Building Lattices by bond distance to an origin site
#           - Building Lattices in a shape around an origin site



#   2) PRE-IMPLEMENTED UNITCELLS (various 2D and 3D stuff)
#
#   3) CONSTRUCTION OF UNITCELLS FROM SITES / CONSTRUCTION HELPER
#
#
#   5) MODIFICATION OF LATTICES
#       - Adding / Removing Connections
#       - Adding / Removing sites
#       - Interaction Strengths
#       - Scaling / Rotating / Shifting in real space
#       - Lattice vector optimization
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
