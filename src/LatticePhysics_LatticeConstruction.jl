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

# module start
# module LatticePhysics_LatticeConstruction


# CONSTRUCTION OF LATTICES FROM UNITCELLS
# BY PUTTING TOGETHER UNITCELLS
include("LatticePhysics_LatticeConstruction/unitcell_patterns_periodic.jl")



# module end
# end
