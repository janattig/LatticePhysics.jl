# SQUARE LATTICE

# Implementation with unknown label type --> raise error
function getUnitcellSquare(
            label_type :: Type{L};
            version :: Int64 = 1,
            unitcell_type :: Type{U} = Unitcell{2,2,L,Site{L,2},Bond{L,2}}
        ) :: U where {L,S<:AbstractSite{L,2},B<:AbstractBond{L,2},U<:AbstractUnitcell{2,2,L,S,B}}

    # give an error
    error("Label type " * string(L) * " of square lattice unitcell not implemented yet")
end

# Implementation with numbered label types (DEFAULT: Int64)
function getUnitcellSquare(
            label_type :: Type{L} = Int64;
            version :: Int64 = 1,
            unitcell_type :: Type{U} = Unitcell{2,2,L,Site{L,2},Bond{L,2}}
        ) :: U where {L<:Number,S<:AbstractSite{L,2},B<:AbstractBond{L,2},U<:AbstractUnitcell{2,2,L,S,B}}
    # return a different unitcell depending on version
    if version == 1
        # return a new Unitcell
        return newUnitcell(
            # lattice vectors
            Vector{Float64}[
                Float64[1, 0],
                Float64[0, 1]
            ],
            # sites
            S[
                newSite(Float64[0,0], L(1), S)
            ],
            # bonds
            B[
                newBond(1,1, L(1), (+1,0), B),
                newBond(1,1, L(1), (-1,0), B),
                newBond(1,1, L(1), (0,+1), B),
                newBond(1,1, L(1), (0,-1), B)
            ],
            # give the unitcell type
            U
        )
    else
        # give an error
        error("Version " * string(version) * " of square lattice unitcell for label type " * string(L) * " not implemented yet")
    end
end
