################################################################################
#
#   HONEYCOMB LATTICE
#
################################################################################



# REFERENCE / FALLBACK (for generic AbstractUnitcell type)

# Referencing to individual functions by wrapping in Val()
function getUnitcellHoneycomb(
            unitcell_type   :: Type{U},
            version         :: Int64 = 1
        ) :: U where {LS,LB,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2},U<:AbstractUnitcell{S,B}}

    # call the respective subfunction by converting to val type
    return getUnitcellHoneycomb(unitcell_type, Val(version))
end

# Fallback for all implementations (if Val{V} is not found)
function getUnitcellHoneycomb(
            unitcell_type   :: Type{U},
            version         :: Val{V}
        ) :: U where {LS,LB,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2},U<:AbstractUnitcell{S,B},V}

    # fallback / fail due to missing implementation
    error("Version " * string(V) * " of honeycomb unitcell (label types " * string(LS) * " / " * string(LB) * ") not implemented yet")
end



# WRAPPER FUNCTIONS (for concrete Unitcell type) call general function

# wrapper function for passing no label types (and version) (DEFAULT)
function getUnitcellHoneycomb(
            version :: Int64 = 1
        ) :: Unitcell{Site{Int64,2},Bond{Int64,2}}

    # create a suitable unitcell of the Unitcell type
    return getUnitcellHoneycomb(Unitcell{Site{Int64,2},Bond{Int64,2}}, version)
end

# wrapper function for passing common label type (and version)
function getUnitcellHoneycomb(
            label_type :: Type{L},
            version    :: Int64 = 1
        ) :: Unitcell{Site{L,2},Bond{L,2}} where L

    # create a suitable unitcell of the Unitcell type
    return getUnitcellHoneycomb(Unitcell{Site{L,2},Bond{L,2}}, version)
end

# wrapper function for passing site / bond label types (and version)
function getUnitcellHoneycomb(
            label_type_site :: Type{LS},
            label_type_bond :: Type{LB},
            version         :: Int64 = 1
        ) :: Unitcell{Site{LS,2},Bond{LB,2}} where {LS,LB}

    # create a suitable unitcell of the Unitcell type
    return getUnitcellHoneycomb(Unitcell{Site{LS,2},Bond{LB,2}}, version)
end



################################################################################
#
#   VERSION IMPLEMENTATIONS FROM HERE ON
#
################################################################################


# Implementation
# - version 1
# - labels <: Any
# --> FALLBACK (raises error)
function getUnitcellHoneycomb(
            unitcell_type :: Type{U},
            version       :: Val{1}
        ) :: U where {LS,LB,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2}, U<:AbstractUnitcell{S,B}}

    # error since this version has no implementation yet
    error("Version 1 of honeycomb unitcell has no implementation for label types " * string(LS) * " / " * string(LB) * " yet")
end

# Implementation
# - version 1
# - labels <: Number / Number
function getUnitcellHoneycomb(
            unitcell_type :: Type{U},
            version       :: Val{1}
        ) :: U where {LS<:Number,LB<:Number,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2}, U<:AbstractUnitcell{S,B}}

    # return a new Unitcell
    return newUnitcell(
        # Type of the unitcell
        U,
        # Bravais lattice vectors
        Vector{Float64}[
            Float64[sqrt(3.0)/2, -0.5],
            Float64[sqrt(3.0)/2, +0.5]
        ],
        # Sites
        S[
            newSite(S, Float64[0.0, 0.0],         LS(1)),
            newSite(S, Float64[1/sqrt(3.0), 0.0], LS(2))
        ],
        # Bonds
        B[
            newBond(B, 1, 2, LB(1), (0, 0)),
            newBond(B, 1, 2, LB(1), (-1, 0)),
            newBond(B, 1, 2, LB(1), (0, -1)),
            newBond(B, 2, 1, LB(1), (0, 0)),
            newBond(B, 2, 1, LB(1), (1, 0)),
            newBond(B, 2, 1, LB(1), (0, 1))
        ]
    )
end

# Implementation
# - version 1
# - labels <: AbstractString / AbstractString
function getUnitcellHoneycomb(
            unitcell_type :: Type{U},
            version       :: Val{1}
        ) :: U where {LS<:AbstractString,LB<:AbstractString,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2}, U<:AbstractUnitcell{S,B}}

    # return a new Unitcell
    return newUnitcell(
        # Type of the unitcell
        U,
        # Bravais lattice vectors
        Vector{Float64}[
            Float64[sqrt(3.0)/2, -0.5],
            Float64[sqrt(3.0)/2, +0.5]
        ],
        # Sites
        S[
            newSite(S, Float64[0.0, 0.0],         LS("1")),
            newSite(S, Float64[1/sqrt(3.0), 0.0], LS("2"))
        ],
        # Bonds
        B[
            newBond(B, 1, 2, LB("1"), (0, 0)),
            newBond(B, 1, 2, LB("1"), (-1, 0)),
            newBond(B, 1, 2, LB("1"), (0, -1)),
            newBond(B, 2, 1, LB("1"), (0, 0)),
            newBond(B, 2, 1, LB("1"), (1, 0)),
            newBond(B, 2, 1, LB("1"), (0, 1))
        ]
    )
end



# Implementation
# - version 4
# - labels <: Any
# --> FALLBACK (raises error)
function getUnitcellHoneycomb(
            unitcell_type :: Type{U},
            version       :: Val{4}
        ) :: U where {LS,LB,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2}, U<:AbstractUnitcell{S,B}}

    # error since this version has no implementation yet
    error("Version 4 of honeycomb unitcell has no implementation for label types " * string(LS) * " / " * string(LB) * " yet")
end

# Implementation
# - version 4
# - labels <: Number
function getUnitcellHoneycomb(
            unitcell_type :: Type{U},
            version       :: Val{4}
        ) :: U where {LS<:Number,LB<:Number,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2}, U<:AbstractUnitcell{S,B}}

    # return a new Unitcell
    return newUnitcell(
        # Type of the unitcell
        U,
        # Bravais lattice vectors
        Vector{Float64}[
            Float64[sqrt(3.0)/2, -0.5],
            Float64[sqrt(3.0)/2, +0.5]
        ],
        # Sites
        S[
            newSite(S, Float64[0.0, 0.0],         LS(1)),
            newSite(S, Float64[1/sqrt(3.0), 0.0], LS(2))
        ],
        # Bonds
        B[
            newBond(B, 1, 2, LB(1), (0, 0)),
            newBond(B, 1, 2, LB(2), (-1, 0)),
            newBond(B, 1, 2, LB(3), (0, -1)),
            newBond(B, 2, 1, LB(1), (0, 0)),
            newBond(B, 2, 1, LB(2), (1, 0)),
            newBond(B, 2, 1, LB(3), (0, 1))
        ]
    )
end

# Implementation
# - version 4
# - labels <: AbstractString
function getUnitcellHoneycomb(
            unitcell_type :: Type{U},
            version       :: Val{4}
        ) :: U where {LS<:AbstractString,LB<:AbstractString,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2}, U<:AbstractUnitcell{S,B}}

    # return a new Unitcell
    return newUnitcell(
        # Type of the unitcell
        U,
        # Bravais lattice vectors
        Vector{Float64}[
            Float64[sqrt(3.0)/2, -0.5],
            Float64[sqrt(3.0)/2, +0.5]
        ],
        # Sites
        S[
            newSite(S, Float64[0.0, 0.0],         LS("1")),
            newSite(S, Float64[1/sqrt(3.0), 0.0], LS("2"))
        ],
        # Bonds
        B[
            newBond(B, 1, 2, LB("x"), (0, 0)),
            newBond(B, 1, 2, LB("y"), (-1, 0)),
            newBond(B, 1, 2, LB("z"), (0, -1)),
            newBond(B, 2, 1, LB("x"), (0, 0)),
            newBond(B, 2, 1, LB("y"), (1, 0)),
            newBond(B, 2, 1, LB("z"), (0, 1))
        ]
    )
end
