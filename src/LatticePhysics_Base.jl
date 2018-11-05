################################################################################
#
#   module LatticePhysics_Base
#
#   --> ABSTRACT TYPE DEFINITION
#           - AbstractSite{D,L}
#           - AbstractBond{N,L}
#           - AbstractUnitcell{D,N,L, S,B}
#           - AbstractLattice{D,N,L, S,B}
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
#   --> NAIVE STRUCT DEFINITIONS
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

# module start
# module LatticePhysics_Base




# ABSTRACT TYPES INCLUDING INTERFACE DEFINITION

# Sites
include("LatticePhysics_Base/abstract_site.jl")
# Bonds
include("LatticePhysics_Base/abstract_bond.jl")
# Unitcells
include("LatticePhysics_Base/abstract_unitcell.jl")
# Lattices
include("LatticePhysics_Base/abstract_lattice.jl")



# CONCRETE TYPES INCLUDING INTERFACE IMPLEMENTATION

# Sites
include("LatticePhysics_Base/concrete_site.jl")
# Bonds
include("LatticePhysics_Base/concrete_bond.jl")
# Unitcells
include("LatticePhysics_Base/concrete_unitcell.jl")
# Lattices
include("LatticePhysics_Base/concrete_lattice.jl")



# module end
# end
