################################################################################
#
#	CONCRETE TYPE
#
#   Lattice{D,N,L,S,B} <: AbstractLattice{D,N,L,S,B}
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
mutable struct Lattice{D,N,L,S,B} <: AbstractLattice{D,N,L,S,B}

    # basis vectors of the Bravais lattice
    lattice_vectors	:: Vector{Vector{Float64}}

    # all sites within the lattice
    sites			:: Vector{S}

    # list of bonds
    bonds			:: Vector{B}

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
            lattice_vectors :: Vector{<:Vector{<:Real}},
            sites           :: Vector{S},
            bonds           :: Vector{B},
            ::Type{Lattice{D,N,L,S,B}}
        ) :: Lattice{D,N,L,S,B} where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N}}

    # return a newly created object
    Lattice{D,N,L,S,B}(lattice_vectors, sites, bonds)
end





# accessing a list of lattice vectors
function latticeVectors(
            lattice :: Lattice{D,N,L,S,B}
        ) :: Vector{Vector{Float64}} where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N}}

    # return the list of lattice vectors
    return lattice.lattice_vectors
end


# accessing a list of sites
function sites(
            lattice :: Lattice{D,N,L,S,B}
        ) :: Vector{S} where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N}}

    # return the list of sites
    return lattice.sites
end


# accessing a list of bonds
function bonds(
            lattice :: Lattice{D,N,L,S,B}
        ) :: Vector{B} where {D,N,L,S<:AbstractSite{L,D},B<:AbstractBond{L,N}}

    # return the list of bonds
    return lattice.bonds
end
