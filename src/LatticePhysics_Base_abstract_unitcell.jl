################################################################################
#
#	ABSTRACT TYPE
#
#   Unitcell{D,N,L,B,S}
#   --> D is the dimension of embedding space
#   --> N is the dimension of the Bravais lattice that the bond is located in
#   --> L is the label type
#   --> B is the bond type (<: AbstractBond{L,N})
#   --> S is the site type (<: AbstractSite{L,D})
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
abstract type AbstractUnitcell{D,N,L,B<:AbstractBond{L,N},S<:AbstractSite{L,D}} end
