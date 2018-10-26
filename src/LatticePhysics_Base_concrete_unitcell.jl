################################################################################
#
#	CONCRETE TYPE
#
#   Unitcell{D,N,L,S,B} <: AbstractUnitcell{D,N,L,S,B}
#   --> D is the dimension of embedding space
#   --> N is the dimension of the Bravais lattice that the bond is located in
#   --> L is the label type
#   --> B is the bond type (<: AbstractBond{L,N})
#   --> S is the site type (<: AbstractSite{L,D})
#
#   FILE CONTAINS
#       - concrete struct definition
#       - interface implementation
#
################################################################################



################################################################################
#
#   CONCRETE STRUCT DEFINITION
#
################################################################################
mutable struct Unitcell{D,N,L,S,B} <: AbstractUnitcell{D,N,L,S,B}

    # basis vectors of the Bravais lattice
    lattice_vectors	:: Vector{Vector{Float64}}

    # basis sites within the unitcell
    sites			:: Vector{S}

    # list of bonds
    bonds			:: Vector{B}

end






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
            ::Type{Unitcell{D,N,L,S,B}}
        ) :: Unitcell{D,N,L,S,B} where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N}}

    # return a newly created object
    Unitcell{D,N,L,S,B}(lattice_vectors, sites, bonds)
end





# accessing a list of lattice vectors
function latticeVectors(
            unitcell :: Unitcell{D,N,L,S,B}
        ) :: Vector{Vector{Float64}} where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N}}

    # return the list of lattice vectors
    return unitcell.lattice_vectors
end


# accessing a list of sites
function sites(
            unitcell :: Unitcell{D,N,L,S,B}
        ) :: Vector{S} where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N}}

    # return the list of sites
    return unitcell.sites
end


# accessing a list of bonds
function bonds(
            unitcell :: Unitcell{D,N,L,S,B}
        ) :: Vector{B} where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N}}

    # return the list of bonds
    return unitcell.bonds
end
