# PRINT  THE CODE THAT GENERATES A UNITCELL

# code for a specific unitcell version
function getUnitcellGeneratingCode(
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

    # build the string of lattice vectors
    lattice_vector_string = ""
    for l in latticeVectors(unitcell)
        lattice_vector_string = lattice_vector_string * "            Float64" * string(l) * ",\n"
    end
    if length(lattice_vector_string) > 2
        lattice_vector_string = lattice_vector_string[1:end-2] * "\n"
    end

    # build a string for all sites
    site_string = ""
    for s in sites(unitcell)
        if sitelabeltype <: AbstractString
            site_string = site_string * "            newSite(Float64" * string(point(s)) * ", LS(\"" * string(label(s)) *  "\"), S),\n"
        else
            site_string = site_string * "            newSite(Float64" * string(point(s)) * ", LS(" * string(label(s)) *  "), S),\n"
        end
    end
    if length(site_string) > 2
        site_string = site_string[1:end-2] * "\n"
    end

    # build a string for all bonds
    bond_string = ""
    for b in bonds(unitcell)
        if bondlabeltype <: AbstractString #newBond(1,1, LB("1"), (+1,0), B)
            bond_string = bond_string * "            newBond(" *
                string(from(b)) * ", " * string(to(b)) *
                ", LB(\"" * string(label(b)) *  "\"), " * string(wrap(b)) * ", B),\n"
        else
            bond_string = bond_string * "            newBond(" *
                string(from(b)) * ", " * string(to(b)) *
                ", LB(" * string(label(b)) *  "), " * string(wrap(b)) * ", B),\n"
        end
    end
    if length(bond_string) > 2
        bond_string = bond_string[1:end-2] * "\n"
    end

    # generate the complete generating code as a string
    generating_code = "function getUnitcell" * uppercase(name[1]) * name[2:end] * "(\n" *
        "            unitcell_type :: Type{U},\n" *
        "            version       :: Val{" * string(version) * "}\n" *
        "        ) :: U where {LS<:" *
        string(sitelabeltype) * ",LB<:" * string(bondlabeltype) *
        ",S<:AbstractSite{LS," * string(D) * "},B<:AbstractBond{LB," * string(N) * "}, " *
        "U<:AbstractUnitcell{S,B}}\n" *
        "    \n" *
        "    # return a new Unitcell\n" *
        "    return newUnitcell(\n" *
        "        # Bravais lattice vectors\n" *
        "        Vector{Float64}[\n" *
        lattice_vector_string *
        "        ],\n" *
        "        # Sites\n" *
        "        S[\n" *
        site_string *
        "        ],\n" *
        "        # Bonds\n" *
        "        B[\n" *
        bond_string *
        "        ],\n" *
        "        # Type of the unitcell\n" *
        "        U\n" *
        "    )\n" *
        "end"

    # print the generating code
    return generating_code
end

# print the code
function printUnitcellGeneratingCode(
            io          :: IO,
            unitcell    :: U,
            name        :: String,
            version     :: Int64 = 1;
            labeltype_site  :: DataType = Nothing,
            labeltype_bond  :: DataType = Nothing,
            print_fallback  :: Bool     = true
        ) where {D,LS,LB,N, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    # get the code
    code = getUnitcellGeneratingCode(unitcell, name, version, labeltype_site=labeltype_site, labeltype_bond=labeltype_bond, print_fallback=print_fallback)

    # print the code
    print(io, code)
end
function printUnitcellGeneratingCode(
            unitcell    :: U,
            name        :: String,
            version     :: Int64 = 1;
            labeltype_site  :: DataType = Nothing,
            labeltype_bond  :: DataType = Nothing,
            print_fallback  :: Bool     = true
        ) where {D,LS,LB,N, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    # get the code
    code = getUnitcellGeneratingCode(unitcell, name, version, labeltype_site=labeltype_site, labeltype_bond=labeltype_bond, print_fallback=print_fallback)

    # print the code
    print(code)
end
