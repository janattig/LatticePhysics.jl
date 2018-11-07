# SQUARE LATTICE


# Implementation with any other label type (Fallback)
function getUnitcellSquare(
            unitcell_type :: Type{U};
            version :: Int64 = 1,
        ) :: U where {LS,LB,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2},U<:AbstractUnitcell{S,B}}

    # give an error
    error("Label types " * string(LS) * " / " * string(LB) * " for square lattice not implemented yet")
end

# Implementation with String label types
function getUnitcellSquare(
            unitcell_type :: Type{U};
            version :: Int64 = 1,
        ) :: U where {LS<:AbstractString,LB<:AbstractString,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2},U<:AbstractUnitcell{S,B}}
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
                newSite(Float64[0,0], LS("1"), S)
            ],
            # bonds
            B[
                newBond(1,1, LB("1"), (+1,0), B),
                newBond(1,1, LB("1"), (-1,0), B),
                newBond(1,1, LB("1"), (0,+1), B),
                newBond(1,1, LB("1"), (0,-1), B)
            ],
            # give the unitcell type
            U
        )
    else
        # give an error
        error("Version " * string(version) * " of square lattice unitcell for label type " * string(LS) * " / " * string(LB) * " not implemented yet")
    end
end

# Implementation with Number label types
function getUnitcellSquare(
            unitcell_type :: Type{U};
            version :: Int64 = 1,
        ) :: U where {LS<:Number,LB<:Number,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2},U<:AbstractUnitcell{S,B}}
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
                newSite(Float64[0,0], LS(1), S)
            ],
            # bonds
            B[
                newBond(1,1, LB(1), (+1,0), B),
                newBond(1,1, LB(1), (-1,0), B),
                newBond(1,1, LB(1), (0,+1), B),
                newBond(1,1, LB(1), (0,-1), B)
            ],
            # give the unitcell type
            U
        )
    else
        # give an error
        error("Version " * string(version) * " of square lattice unitcell for label type " * string(LS) * " / " * string(LB) * " not implemented yet")
    end
end


# wrapper function for only passing the label type (DEFAULT: Int64)
function getUnitcellSquare(
            label_type  :: Type{L}  = Int64;
            version     :: Int64    = 1,
        ) :: Unitcell{Site{L,2},Bond{L,2}} where L
    # create a suitable unitcell of the given type
    return getUnitcellSquare(Unitcell{Site{L,2},Bond{L,2}}, version=version)
end
