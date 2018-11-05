################################################################################
#
#	ABSTRACT TYPE
#
#   Lattice{D,N,L,S,B,U}
#   --> D is the dimension of embedding space
#   --> N is the dimension of the Bravais lattice that the bond is located in
#   --> L is the label type
#   --> S is the site type (<: AbstractSite{L,D})
#   --> B is the bond type (<: AbstractBond{L,N})
#
#   FILE CONTAINS
#       - abstract type definition
#       - interface definition
#       - interface testing
#
################################################################################


################################################################################
#
#   ABSTRACT TYPE DEFINITION
#
################################################################################
abstract type AbstractLattice{D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N}} end
