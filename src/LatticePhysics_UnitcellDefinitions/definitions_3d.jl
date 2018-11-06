# CUBIC (SIMPLE CUBIC / SC) LATTICE


# Implementation with any other label type (Fallback)
function getUnitcellCubic(
            unitcell_type :: Type{U};
            version :: Int64 = 1,
        ) :: U where {L,S<:AbstractSite{L,3},B<:AbstractBond{L,3},U<:AbstractUnitcell{3,3,L,S,B}}

    # give an error
    error("Label type " * string(L) * " of cubic lattice unitcell not implemented yet")
end

# Implementation with String label types
function getUnitcellCubic(
            unitcell_type :: Type{U};
            version :: Int64 = 1,
        ) :: U where {L<:AbstractString,S<:AbstractSite{L,3},B<:AbstractBond{L,3},U<:AbstractUnitcell{3,3,L,S,B}}
    # return a different unitcell depending on version
    if version == 1
        # return a new Unitcell
        return newUnitcell(
            # lattice vectors
            Vector{Float64}[
                Float64[1, 0, 0],
                Float64[0, 1, 0],
                Float64[0, 0, 1]
            ],
            # sites
            S[
                newSite(Float64[0,0,0], L("1"), S)
            ],
            # bonds
            B[
                newBond(1,1, L("1"), (+1,0,0), B),
                newBond(1,1, L("1"), (-1,0,0), B),
                newBond(1,1, L("1"), (0,+1,0), B),
                newBond(1,1, L("1"), (0,-1,0), B),
                newBond(1,1, L("1"), (0,0,+1), B),
                newBond(1,1, L("1"), (0,0,-1), B)
            ],
            # give the unitcell type
            U
        )
    else
        # give an error
        error("Version " * string(version) * " of cubic lattice unitcell for label type " * string(L) * " not implemented yet")
    end
end

# Implementation with Number label types
function getUnitcellCubic(
            unitcell_type :: Type{U};
            version :: Int64 = 1,
        ) :: U where {L<:Number,S<:AbstractSite{L,3},B<:AbstractBond{L,3},U<:AbstractUnitcell{3,3,L,S,B}}
    # return a different unitcell depending on version
    if version == 1
        # return a new Unitcell
        return newUnitcell(
            # lattice vectors
            Vector{Float64}[
                Float64[1, 0, 0],
                Float64[0, 1, 0],
                Float64[0, 0, 1]
            ],
            # sites
            S[
                newSite(Float64[0,0,0], L(1), S)
            ],
            # bonds
            B[
                newBond(1,1, L(1), (+1,0,0), B),
                newBond(1,1, L(1), (-1,0,0), B),
                newBond(1,1, L(1), (0,+1,0), B),
                newBond(1,1, L(1), (0,-1,0), B),
                newBond(1,1, L(1), (0,0,+1), B),
                newBond(1,1, L(1), (0,0,-1), B)
            ],
            # give the unitcell type
            U
        )
    else
        # give an error
        error("Version " * string(version) * " of cubic lattice unitcell for label type " * string(L) * " not implemented yet")
    end
end


# wrapper function for only passing the label type (DEFAULT: Int64)
function getUnitcellCubic(
            label_type  :: Type{L}  = Int64;
            version     :: Int64    = 1,
        ) :: Unitcell{3,3,L,Site{L,3},Bond{L,3}} where L
    # create a suitable unitcell of the given type
    return getUnitcellCubic(Unitcell{3,3,L,Site{L,3},Bond{L,3}}, version=version)
end
