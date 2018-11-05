################################################################################
#
#	ABSTRACT TYPE
#
#   Lattice{D,N,L,S,B}
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






################################################################################
#
#	INTERFACING / ACCESSING LATTICES
#	(functions have to be overwritten by concrete types)
#
################################################################################


# default constructor interface
# used for creation of new lattices
function newLattice(
            lattice_vectors :: Vector{<:Vector{<:Real}},
            sites           :: Vector{S},
            bonds           :: Vector{B},
            ::Type{LA}
        ) :: LA where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N},LA<:AbstractLattice{D,N,L,S,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'newLattice' for concrete lattice type " *
            string(LA) * " with bond type " * string(B) *
            " and site type " * string(S)   )
end




# accessing a list of lattice vectors
function latticeVectors(
            lattice :: LA
        ) :: Vector{Vector{Float64}} where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N},LA<:AbstractLattice{D,N,L,S,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'latticeVectors' for concrete lattice type " * string(LA) )
end


# accessing a list of sites
function sites(
            lattice :: LA
        ) :: Vector{S} where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N},LA<:AbstractLattice{D,N,L,S,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'sites' for concrete lattice type " *
            string(LA) * " with site type " * string(S)   )
end


# accessing a list of bonds
function bonds(
            lattice :: LA
        ) :: Vector{B} where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N},LA<:AbstractLattice{D,N,L,S,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'bonds' for concrete lattice type " *
            string(LA) * " with bond type " * string(B)   )
end
