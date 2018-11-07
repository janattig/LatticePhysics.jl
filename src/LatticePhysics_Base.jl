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



# CUSTOM NICE PRINTING
include("LatticePhysics_Base/show.jl")


# module end
# end
