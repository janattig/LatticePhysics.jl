################################################################################
#
#	CONCRETE TYPE
#
#   Lattice{S,B,U} <: AbstractLattice{S,B,U}
#   --> S is the site type (<: AbstractSite{LS,D})
#       --> D is the dimension of embedding space
#       --> LS is the label type of sites
#   --> B is the bond type (<: AbstractBond{LB,N})
#       --> N is the dimension of the Bravais lattice that the bond is located in
#       --> LB is the label type of bonds
#   --> U is the unitcell type (<: AbstractUnitcell{SU,BU})
#       --> SU is the site type of the unitcell (<: AbstractSite{LUS,DU})
#           --> DU is the dimension of embedding space (of the UC)
#           --> LUS is the label type of sites (of the UC)
#       --> BU is the bond type of the unitcell (<: AbstractBond{LUB,NU})
#           --> NU is the dimension of the Bravais lattice that the bond is located in (of the UC)
#           --> LUB is the label type of bonds (of the UC)
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
mutable struct Lattice{S,B,U} <: AbstractLattice{S,B,U}

    # basis vectors of the Bravais lattice
    lattice_vectors	:: Vector{Vector{Float64}}

    # all sites within the lattice
    sites			:: Vector{S}

    # list of bonds
    bonds			:: Vector{B}

    # unitcell
    unitcell        :: U

end






################################################################################
#
#	INTERFACING / ACCESSING LATTICES
#	(functions have to be overwritten by concrete types)
#
################################################################################


# default constructor interface
# used for creation of new lattices
function newLattice(
            ::Type{Lattice{S,B,U}},
            lattice_vectors :: Vector{<:Vector{<:Real}},
            sites           :: Vector{S},
            bonds           :: Vector{B},
            unitcell        :: U
        ) :: Lattice{S,B,U} where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # return a newly created object
    Lattice{S,B,U}(lattice_vectors, sites, bonds, unitcell)
end





# accessing a list of lattice vectors
function latticeVectors(
            lattice :: Lattice{S,B,U}
        ) :: Vector{Vector{Float64}} where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # return the list of lattice vectors
    return lattice.lattice_vectors
end
# setting a list of lattice vectors
function latticeVectors!(
            lattice         :: Lattice{S,B,U},
            lattice_vectors :: Vector{<:Vector{<:Real}}
        ) where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # set the list of lattice vectors
    lattice.lattice_vectors = lattice_vectors
end


# accessing a list of sites
function sites(
            lattice :: Lattice{S,B,U}
        ) :: Vector{S} where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # return the list of sites
    return lattice.sites
end
# setting a list of sites
function sites!(
            lattice :: Lattice{S,B,U},
            sites   :: Vector{S}
        ) where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # set the list of sites
    lattice.sites = sites
end

# accessing a list of bonds
function bonds(
            lattice :: Lattice{S,B,U}
        ) :: Vector{B} where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # return the list of bonds
    return lattice.bonds
end
# setting a list of bonds
function bonds!(
            lattice :: Lattice{S,B,U},
            bonds   :: Vector{B}
        ) where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # set the list of bonds
    lattice.bonds = bonds
end


# accessing the unitcell
function unitcell(
            lattice :: Lattice{S,B,U}
        ) :: U where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # return the unitcell
    return lattice.unitcell
end
# setting the unitcell
function unitcell!(
            lattice  :: Lattice{S,B,U},
            unitcell :: U
        ) :: U where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N}}

    # set the unitcell
    lattice.unitcell = unitcell
end
