################################################################################
#
#	ABSTRACT TYPE
#
#   Unitcell{S,B}
#   --> S is the site type (<: AbstractSite{LS,D})
#       --> D is the dimension of embedding space
#       --> LS is the label type of sites
#   --> B is the bond type (<: AbstractBond{LB,N})
#       --> N is the dimension of the Bravais lattice that the bond is located in
#       --> LB is the label type of bonds
#
#   FILE CONTAINS
#       - abstract type definition
#       - interface definition
#       - TODO interface testing
#
################################################################################

################################################################################
#
#   ABSTRACT TYPE DEFINITION
#
################################################################################
abstract type AbstractUnitcell{
        S <: AbstractSite{LS,D} where {LS,D},
        B <: AbstractBond{LB,N} where {LB,N}
    } end






################################################################################
#
#	INTERFACING / ACCESSING UNITCELLS
#	(functions have to be overwritten by concrete types)
#
################################################################################


# default constructor interface
# used for creation of new unitcells
function newUnitcell(
            lattice_vectors :: Vector{<:Vector{<:Real}},
            sites           :: Vector{S},
            bonds           :: Vector{B},
            ::Type{U}
        ) :: U where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'newUnitcell' for concrete unitcell type " *
            string(U) * " with bond type " * string(B) *
            " and site type " * string(S)   )
end





# accessing a list of lattice vectors
function latticeVectors(
            unitcell :: U
        ) :: Vector{Vector{Float64}} where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'latticeVectors' for concrete unitcell type " * string(U) )
end


# accessing a list of sites
function sites(
            unitcell :: U
        ) :: Vector{S} where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'sites' for concrete unitcell type " *
            string(U) * " with site type " * string(S)   )
end


# accessing a list of bonds
function bonds(
            unitcell :: U
        ) :: Vector{B} where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'bonds' for concrete unitcell type " *
            string(U) * " with bond type " * string(B)   )
end








# some more beauty interface
# builds on interface defined above but can also be overwritten

# number of sites
function numSites(
            unitcell :: U
        ) :: Int64 where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # return length of the site array that the unitcell type implements
    return length(sites(unitcell))
end

# number of bonds
function numBonds(
            unitcell :: U
        ) :: Int64 where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # return length of the bond array that the unitcell type implements
    return length(bonds(unitcell))
end




# access a specific site or bond

# site
function site(
            unitcell :: U,
            index    :: Int64
        ) :: S where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # return the respective site
    return sites(unitcell)[index]
end

# bond
function bond(
            unitcell :: U,
            index    :: Int64
        ) :: B where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # return the respective bond
    return bonds(unitcell)[index]
end






# specific Bravais lattices

# a1
function a1(
            unitcell :: U
        ) :: Vector{Float64} where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # return the first entry in the lattice vectors list
    return latticeVectors(unitcell)[1]
end

# a2
function a2(
            unitcell :: U
        ) :: Vector{Float64} where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # return the second entry in the lattice vectors list
    return latticeVectors(unitcell)[2]
end

# a3
function a3(
            unitcell :: U
        ) :: Vector{Float64} where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # return the third entry in the lattice vectors list
    return latticeVectors(unitcell)[3]
end
