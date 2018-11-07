################################################################################
#
#   CUBIC (SC / SIMPLE CUBIC) LATTICE
#
################################################################################

# Implementation with any unitcell type and version (Fallback)
function getUnitcellCubic(
            unitcell_type   :: Type{U};
            version         :: Int64 = 1,
        ) :: U where {LS,LB,S<:AbstractSite{LS,3},B<:AbstractBond{LB,3},U<:AbstractUnitcell{S,B}}

    # this might fail due to missing implementation
    try
        # try to call the respective subfunction
        return getUnitcellCubic(unitcell_type, Val(version))
    catch error_caught
        # check and possibliy give the error directly
        if isa(error_caught, MethodError)
            # print that there is a missing function
            error("Version " * string(version) * " of cubic lattice unitcell (label types " * string(LS) * " / " * string(LB) * ") not implemented yet")
        else
            # throw the error again
            throw(error_caught)
        end
    end
end

# wrapper function for passing the label type (DEFAULT: Int64)
function getUnitcellCubic(
            label_type  :: Type{L}  = Int64;
            version     :: Int64    = 1,
        ) :: Unitcell{Site{L,3},Bond{L,3}} where L
    # create a suitable unitcell of the given type
    return getUnitcellCubic(Unitcell{Site{L,3},Bond{L,3}}, version=version)
end





################################################################################
#
#   VERSION IMPLEMENTATIONS FROM HERE ON
#
################################################################################

# Implementation
# - version 1
# - labels <: Any
# --> Fallback (raises error)
function getUnitcellCubic(
            unitcell_type :: Type{U},
            version :: Val{1}
        ) :: U where {LS,LB,S<:AbstractSite{LS,3},B<:AbstractBond{LB,3},U<:AbstractUnitcell{S,B}}

    # return a new Unitcell
    error("Version 1 of cubic lattice unitcell has no implementation for label types " * string(LS) * " / " * string(LB) * " yet")
end

# Implementation
# - version 1
# - labels <: AbstractString
function getUnitcellCubic(
            unitcell_type :: Type{U},
            version :: Val{1}
        ) :: U where {LS<:AbstractString,LB<:AbstractString,S<:AbstractSite{LS,3},B<:AbstractBond{LB,3},U<:AbstractUnitcell{S,B}}

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
            newSite(Float64[0,0,0], LS("1"), S)
        ],
        # bonds
        B[
            newBond(1,1, LB("1"), (+1,0,0), B),
            newBond(1,1, LB("1"), (-1,0,0), B),
            newBond(1,1, LB("1"), (0,+1,0), B),
            newBond(1,1, LB("1"), (0,-1,0), B),
            newBond(1,1, LB("1"), (0,0,+1), B),
            newBond(1,1, LB("1"), (0,0,-1), B)
        ],
        # give the unitcell type
        U
    )
end

# Implementation
# - version 1
# - labels <: Number
function getUnitcellCubic(
            unitcell_type :: Type{U},
            version :: Val{1}
        ) :: U where {LS<:Number,LB<:Number,S<:AbstractSite{LS,3},B<:AbstractBond{LB,3},U<:AbstractUnitcell{S,B}}

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
            newSite(Float64[0,0,0], LS(1), S)
        ],
        # bonds
        B[
            newBond(1,1, LB(1), (+1,0,0), B),
            newBond(1,1, LB(1), (-1,0,0), B),
            newBond(1,1, LB(1), (0,+1,0), B),
            newBond(1,1, LB(1), (0,-1,0), B),
            newBond(1,1, LB(1), (0,0,+1), B),
            newBond(1,1, LB(1), (0,0,-1), B)
        ],
        # give the unitcell type
        U
    )
end
