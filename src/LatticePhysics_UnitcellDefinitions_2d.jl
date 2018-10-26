# SQUARE LATTICE

# Implementation with unknown label type
function getUnitcellSquare(
            label_type :: Type{L};
            version :: Int64 = 1,
            unitcell_type :: Type{U} = Unitcell{2,2,L,Site{L,2},Bond{L,2}}
        ) :: U where {L,S<:AbstractSite{L,2},B<:AbstractBond{L,2},U<:AbstractUnitcell{2,2,L,S,B}}

    # give an error
    error("Label type " * string(L) * " of square lattice unitcell not implemented yet")
end

# Implementation with Int64 label type
function getUnitcellSquare(
            label_type :: Type{Int64};
            version :: Int64 = 1,
            unitcell_type :: Type{U} = Unitcell{2,2,Int64,Site{Int64,2},Bond{Int64,2}}
        ) :: U where {S<:AbstractSite{Int64,2},B<:AbstractBond{Int64,2},U<:AbstractUnitcell{2,2,Int64,S,B}}
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
                newSite(Float64[0,0], 1, S)
            ],
            # bonds
            B[
                newBond(1,1, 1, (+1,0), B),
                newBond(1,1, 1, (-1,0), B),
                newBond(1,1, 1, (0,+1), B),
                newBond(1,1, 1, (0,-1), B)
            ],
            # give the unitcell type
            U
        )
    else
        # give an error
        error("Version " * string(version) * " of square lattice unitcell for label type Int64 not implemented yet")
    end
end
