# PRINT  THE CODE THAT GENERATES A UNITCELL
function printUnitcellGeneratingCode(
            unitcell    :: U,
            name        :: String,
            version     :: Int64 = 1;
            labeltype_site  :: DataType = Nothing,
            labeltype_bond  :: DataType = Nothing,
            print_fallback  :: Bool     = true
        ) where {D,LS,LB,N, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}
    # set the correct label type for sites and bonds
    sitelabeltype = (labeltype_site == Nothing ? LS : labeltype_site) :: DataType
    bondlabeltype = (labeltype_bond == Nothing ? LB : labeltype_bond) :: DataType
    # generate the complete generating code as a string
    generating_code = "function getUnitcell" * uppercase(name[1]) * name[2:end] * "(\n" *
        "            unitcell_type :: Type{U},\n" *
        "            version :: Val{1}\n" *
        "        ) :: U where {LS<:" *
        string(sitelabeltype) * ",LB<:" * string(bondlabeltype) *
        ",S<:AbstractSite{LS," * string(D) * "},B<:AbstractBond{LB," * string(N) *
        "},U<:AbstractUnitcell{S,B}}"
    # print the generating code
    print(generating_code)
end

function getUnitcellSquare(
            unitcell_type :: Type{U},
            version :: Val{1}
        ) :: U where {LS<:AbstractString,LB<:AbstractString,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2},U<:AbstractUnitcell{S,B}}

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
end
