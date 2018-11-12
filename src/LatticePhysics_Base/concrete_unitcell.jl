################################################################################
#
#	CONCRETE TYPE
#
#   Unitcell{D,N,L,S,B} <: AbstractUnitcell{S,B}
#   --> S is the site type (<: AbstractSite{LS,D})
#       --> D is the dimension of embedding space
#       --> LS is the label type of sites
#   --> B is the bond type (<: AbstractBond{LB,N})
#       --> N is the dimension of the Bravais lattice that the bond is located in
#       --> LB is the label type of bonds
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
mutable struct Unitcell{S,B} <: AbstractUnitcell{S,B}

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
            :: Type{Unitcell{S,B}},
            lattice_vectors :: Vector{<:Vector{<:Real}},
            sites           :: Vector{S},
            bonds           :: Vector{B}
        ) :: Unitcell{S,B} where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # return a newly created object
    Unitcell{S,B}(lattice_vectors, sites, bonds)
end





# accessing a list of lattice vectors
function latticeVectors(
            unitcell :: Unitcell{S,B}
        ) :: Vector{Vector{Float64}} where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # return the list of lattice vectors
    return unitcell.lattice_vectors
end
# setting a list of lattice vectors
function latticeVectors!(
            unitcell        :: Unitcell{S,B},
            lattice_vectors :: Vector{<:Vector{<:Real}}
        ) where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # set the list of lattice vectors
    unitcell.lattice_vectors = lattice_vectors
end


# accessing a list of sites
function sites(
            unitcell :: Unitcell{S,B}
        ) :: Vector{S} where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # return the list of sites
    return unitcell.sites
end
# setting a list of sites
function sites!(
            unitcell :: Unitcell{S,B},
            sites    :: Vector{S}
        ) where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # set the list of sites
    unitcell.sites = sites
end


# accessing a list of bonds
function bonds(
            unitcell :: Unitcell{S,B}
        ) :: Vector{B} where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # return the list of bonds
    return unitcell.bonds
end
# setting a list of bonds
function bonds!(
            unitcell :: Unitcell{S,B},
            bonds    :: Vector{B}
        ) where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # set the list of bonds
    unitcell.bonds = bonds
end
