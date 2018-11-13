################################################################################
#
#	ABSTRACT TYPE
#
#   Lattice{S,B,U}
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
abstract type AbstractLattice{
        S <: AbstractSite{LS,D} where {LS,D},
        B <: AbstractBond{LB,N} where {LB,N},
        U <: AbstractUnitcell{SU,BU} where {
            SU <: AbstractSite{LUS,DU} where {LUS,DU},
            BU <: AbstractBond{LUB,NU} where {LUB,NU},
        }
    } end






################################################################################
#
#	INTERFACING / ACCESSING LATTICES
#	(functions have to be overwritten by concrete types)
#
################################################################################


# default constructor interface
# used for creation of new lattices
function newLattice(
            ::Type{L},
            lattice_vectors :: Vector{<:Vector{<:Real}},
            sites           :: Vector{S},
            bonds           :: Vector{B},
            unitcell        :: U
        ) :: L where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'newLattice' for concrete lattice type " *
            string(L) * " with bond type " * string(B) *
            " and site type " * string(S)   )
end




# accessing a list of lattice vectors
function latticeVectors(
            lattice :: L
        ) :: Vector{Vector{Float64}} where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'latticeVectors' for concrete lattice type " * string(L) )
end
# setting a list of lattice vectors
function latticeVectors!(
            lattice         :: L,
            lattice_vectors :: Vector{<:Vector{<:Real}}
        ) where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'latticeVectors!' for concrete lattice type " * string(L) )
end


# accessing a list of sites
function sites(
            lattice :: L
        ) :: Vector{S} where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'sites' for concrete lattice type " *
            string(L) * " with site type " * string(S)   )
end
# setting a list of sites
function sites!(
            lattice :: L,
            sites   :: Vector{S}
        ) where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'sites!' for concrete lattice type " *
            string(L) * " with site type " * string(S)   )
end

# accessing a list of bonds
function bonds(
            lattice :: L
        ) :: Vector{B} where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'bonds' for concrete lattice type " *
            string(L) * " with bond type " * string(B)   )
end
# setting a list of bonds
function bonds!(
            lattice :: L,
            bonds   :: Vector{B}
        ) where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'bonds!' for concrete lattice type " *
            string(L) * " with bond type " * string(B)   )
end


# accessing the unitcell
function unitcell(
            lattice :: L
        ) :: U where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'unitcell' for concrete lattice type " *
            string(L) * " with unitcell type " * string(U)   )
end
# setting the unitcell
function unitcell!(
            lattice  :: L,
            unitcell :: U
        ) where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'unitcell!' for concrete lattice type " *
            string(L) * " with unitcell type " * string(U)   )
end





# SIMILAR FUNCTION (can be overwritten but does not have to be overwritten)

# without new parameters
function similar(
            l :: L
        ) where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # return a new lattice object
    return newLattice(
        L,
        deepcopy(latticeVectors(l)),
        deepcopy(sites(l)),
        deepcopy(bonds(l)),
        deepcopy(unitcell(l))
    )
end
# with new parameters
function similar(
            l :: L,
            lattice_vectors :: Vector{<:Vector{<:Real}},
            sites           :: Vector{S},
            bonds           :: Vector{B},
            unitcell        :: U
        ) where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # create a new lattice object
    l_new = similar(l)
    # set parameters
    latticeVectors!(l_new, lattice_vectors)
    sites!(l_new, sites)
    bonds!(l_new, bonds)
    unitcell!(l_new, unitcell)
    # return the new object
    return l_new
end





# some more beauty interface
# builds on interface defined above but can also be overwritten

# number of sites
function numSites(
            lattice :: L
        ) :: Int64 where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # return length of the site array that the lattice type implements
    return length(sites(lattice))
end

# number of bonds
function numBonds(
            lattice :: L
        ) :: Int64 where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # return length of the bond array that the lattice type implements
    return length(bonds(lattice))
end





# access a specific site or bond

# site
function site(
            lattice :: L,
            index   :: Int64
        ) :: S where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # return the respective site
    return sites(lattice)[index]
end

# bond
function bond(
            lattice :: L,
            index   :: Int64
        ) :: B where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # return the respective bond
    return bonds(lattice)[index]
end







# specific Bravais lattices

# a1
function a1(
            lattice :: L
        ) :: Vector{Float64} where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # return the first entry in the lattice vectors list
    return latticeVectors(lattice)[1]
end

# a2
function a2(
            lattice :: L
        ) :: Vector{Float64} where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # return the second entry in the lattice vectors list
    return latticeVectors(lattice)[2]
end

# a3
function a3(
            lattice :: L
        ) :: Vector{Float64} where {D,N,LS,LB,U,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},L<:AbstractLattice{S,B,U}}

    # return the third entry in the lattice vectors list
    return latticeVectors(lattice)[3]
end
