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
abstract type AbstractUnitcell{D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N}} end






################################################################################
#
#	INTERFACING / ACCESSING UNITCELLS
#	(functions have to be overwritten by concrete types)
#
################################################################################


# default constructor interface
# used for creation of new unitcells
function newUnitcell(
            lattice_vectors :: Vector{Vector{<:Real}},
            sites           :: Vector{S},
            bonds           :: Vector{B},
            ::Type{U}
        ) :: U where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N},U<:AbstractUnitcell{D,N,L,S,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'newUnitcell' for concrete unitcell type " *
            string(U) * " with bond type " * string(B) *
            " and site type " * string(S)   )
end





# accessing a list of lattice vectors
function latticeVectors(
            unitcell :: U
        ) :: Vector{Vector{Float64}} where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N},U<:AbstractUnitcell{D,N,L,S,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'latticeVectors' for concrete unitcell type " * string(U) )
end


# accessing a list of sites
function sites(
            unitcell :: U
        ) :: Vector{S} where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N},U<:AbstractUnitcell{D,N,L,S,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'sites' for concrete unitcell type " *
            string(U) * " with site type " * string(S)   )
end


# accessing a list of bonds
function bonds(
            unitcell :: U
        ) :: Vector{B} where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N},U<:AbstractUnitcell{D,N,L,S,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'bonds' for concrete unitcell type " *
            string(U) * " with bond type " * string(B)   )
end
