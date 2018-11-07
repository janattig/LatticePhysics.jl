################################################################################
#
#   TRIANGULAR LATTICE
#
################################################################################

# Implementation with any unitcell type and version (Fallback)
function getUnitcellTriangular(
            unitcell_type   :: Type{U};
            version         :: Int64 = 1,
        ) :: U where {LS,LB,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2},U<:AbstractUnitcell{S,B}}

    # this might fail due to missing implementation
    try
        # try to call the respective subfunction
        return getUnitcellTriangular(unitcell_type, Val(version))
    catch error_caught
        # check and possibliy give the error directly
        if isa(error_caught, MethodError)
            # print that there is a missing function
            error("Version " * string(version) * " of triangular unitcell (label types " * string(LS) * " / " * string(LB) * ") not implemented yet")
        else
            # throw the error again
            throw(error_caught)
        end
    end
end

# wrapper function for passing a common label type (DEFAULT: Int64)
function getUnitcellTriangular(
            label_type  :: Type{L}  = Int64;
            version     :: Int64    = 1,
        ) :: Unitcell{Site{L,2},Bond{L,2}} where L
    # create a suitable unitcell of the given type
    return getUnitcellTriangular(Unitcell{Site{L,2},Bond{L,2}}, version=version)
end

# wrapper function for passing both label types
function getUnitcellTriangular(
            label_type_site :: Type{LS},
            label_type_bond :: Type{LB};
            version         :: Int64 = 1,
        ) :: Unitcell{Site{LS,2},Bond{LB,2}} where {LS,LB}
    # create a suitable unitcell of the given type
    return getUnitcellTriangular(Unitcell{Site{LS,2},Bond{LB,2}}, version=version)
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
function getUnitcellTriangular(
            unitcell_type :: Type{U},
            version :: Val{1}
        ) :: U where {LS,LB,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2},U<:AbstractUnitcell{S,B}}

    # return a new Unitcell
    error("Version 1 of triangular unitcell has no implementation for label types " * string(LS) * " / " * string(LB) * " yet")
end

# Implementation
# - version 1
# - labels <: AbstractString
function getUnitcellTriangular(
            unitcell_type :: Type{U},
            version :: Val{1}
        ) :: U where {LS<:AbstractString,LB<:AbstractString,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2},U<:AbstractUnitcell{S,B}}

    # return a new Unitcell
    return newUnitcell(
        # lattice vectors
        Vector{Float64}[
            Float64[sqrt(3.0)/2, -0.5],
            Float64[sqrt(3.0)/2, +0.5]
        ],
        # sites
        S[
            newSite(Float64[0,0], LS("1"), S)
        ],
        # bonds
        B[
            newBond(1,1, LB("1"), (+1, 0), B),
            newBond(1,1, LB("1"), (-1, 0), B),
            newBond(1,1, LB("1"), ( 0,+1), B),
            newBond(1,1, LB("1"), ( 0,-1), B),
            newBond(1,1, LB("1"), (+1,-1), B),
            newBond(1,1, LB("1"), (-1,+1), B)
        ],
        # give the unitcell type
        U
    )
end

# Implementation
# - version 1
# - labels <: Number
function getUnitcellTriangular(
            unitcell_type :: Type{U},
            version :: Val{1}
        ) :: U where {LS<:Number,LB<:Number,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2},U<:AbstractUnitcell{S,B}}

    # return a new Unitcell
    return newUnitcell(
        # lattice vectors
        Vector{Float64}[
            Float64[sqrt(3.0)/2, -0.5],
            Float64[sqrt(3.0)/2, +0.5]
        ],
        # sites
        S[
            newSite(Float64[0,0], LS(1), S)
        ],
        # bonds
        B[
            newBond(1,1, LB(1), (+1, 0), B),
            newBond(1,1, LB(1), (-1, 0), B),
            newBond(1,1, LB(1), ( 0,+1), B),
            newBond(1,1, LB(1), ( 0,-1), B),
            newBond(1,1, LB(1), (+1,-1), B),
            newBond(1,1, LB(1), (-1,+1), B)
        ],
        # give the unitcell type
        U
    )
end
