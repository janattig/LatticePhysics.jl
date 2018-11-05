# CUSTOM NICE PRINTING FOR VARIOUS TYPES

# import the relevant function
import Base.show


# ABSTRACT TYPES

# single SITE
function Base.show(io::IO, s::S) where {L,D,S<:AbstractSite{L,D}}
    print(io, S, " @", point(s), ": ", label(s))
end

# single BOND
function Base.show(io::IO, b::B) where {L,N,B<:AbstractBond{L,N}}
    print(io, B, " ", from(b), "-->", to(b), " @", wrap(b), ": ", label(b))
end

# single UNITCELL
function Base.show(io::IO, u::U) where {D,L,N,S<:AbstractSite{L,D},B<:AbstractBond{L,N},U<:AbstractUnitcell{D,N,L,S,B}}
    print(io, "Unitcell object\n--> type ", U, "\n--> ", length(sites(u)), " sites of type ", S, "\n--> ", length(bonds(u)), " bonds of type ", B)
end

# single LATTICE
function Base.show(io::IO, la::LA) where {D,L,N,S<:AbstractSite{L,D},B<:AbstractBond{L,N},LA<:AbstractLattice{D,N,L,S,B}}
    print(io, "Lattice object\n--> type ", U, "\n--> ", length(sites(la)), " sites of type ", S, "\n--> ", length(bonds(la)), " bonds of type ", B)
end
