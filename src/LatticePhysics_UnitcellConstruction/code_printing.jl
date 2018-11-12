# PRINT  THE CODE THAT GENERATES A UNITCELL

# code for a specific unitcell version
function getCodeUnitcellVersion(
            unitcell    :: U,
            name        :: String = "myunitcell",
            version     :: Int64  = 1;
            labeltype_site  :: DataType = Nothing,
            labeltype_bond  :: DataType = Nothing
        ) :: String where {D,LS,LB,N, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

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
            site_string = site_string * "            newSite(S, Float64" * string(point(s)) * ", LS(\"" * string(label(s)) *  "\")),\n"
        else
            site_string = site_string * "            newSite(S, Float64" * string(point(s)) * ", LS(" * string(label(s)) *  ")),\n"
        end
    end
    if length(site_string) > 2
        site_string = site_string[1:end-2] * "\n"
    end

    # build a string for all bonds
    bond_string = ""
    for b in bonds(unitcell)
        if bondlabeltype <: AbstractString #newBond(B, 1,1, LB("1"), (+1,0))
            bond_string = bond_string * "            newBond(B, " *
                string(from(b)) * ", " * string(to(b)) *
                ", LB(\"" * string(label(b)) *  "\"), " * string(wrap(b)) * "),\n"
        else
            bond_string = bond_string * "            newBond(B, " *
                string(from(b)) * ", " * string(to(b)) *
                ", LB(" * string(label(b)) *  "), " * string(wrap(b)) * "),\n"
        end
    end
    if length(bond_string) > 2
        bond_string = bond_string[1:end-2] * "\n"
    end

    # generate the complete generating code as a string
    generating_code = "" *
        "# Implementation\n" *
        "# - version " * string(version) * "\n" *
        "# - labels <: " * string(sitelabeltype) * " / " * string(bondlabeltype) * "\n" *
        "function getUnitcell" * uppercase(name[1]) * name[2:end] * "(\n" *
        "            unitcell_type :: Type{U},\n" *
        "            version       :: Val{" * string(version) * "}\n" *
        "        ) :: U where {LS<:" *
        string(sitelabeltype) * ",LB<:" * string(bondlabeltype) *
        ",S<:AbstractSite{LS," * string(D) * "},B<:AbstractBond{LB," * string(N) * "}, " *
        "U<:AbstractUnitcell{S,B}}\n" *
        "    \n" *
        "    # return a new Unitcell\n" *
        "    return newUnitcell(\n" *
        "        # Type of the unitcell\n" *
        "        U,\n" *
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
        "        ]\n" *
        "    )\n" *
        "end"

    # print the generating code
    return generating_code
end
# code for a specific unitcell version but simplified
function getCodeUnitcellSimplifiedVersion(
            unitcell    :: U,
            name        :: String = "myunitcell",
            version     :: Int64  = 1;
            labeltype_site  :: DataType = Nothing,
            labeltype_bond  :: DataType = Nothing
        ) :: String where {D,LS,LB,N, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

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
    for (i,s) in enumerate(sites(unitcell))
        if sitelabeltype <: AbstractString
            site_string = site_string * "            newSite(S, Float64" * string(point(s)) * ", LS(\"" * string(i) *  "\")),\n"
        else
            site_string = site_string * "            newSite(S, Float64" * string(point(s)) * ", LS(" * string(i) *  ")),\n"
        end
    end
    if length(site_string) > 2
        site_string = site_string[1:end-2] * "\n"
    end

    # build a string for all bonds
    bond_string = ""
    for b in bonds(unitcell)
        if bondlabeltype <: AbstractString #newBond(B, 1,1, LB("1"), (+1,0))
            bond_string = bond_string * "            newBond(B, " *
                string(from(b)) * ", " * string(to(b)) *
                ", LB(\"1\"), " * string(wrap(b)) * "),\n"
        else
            bond_string = bond_string * "            newBond(B, " *
                string(from(b)) * ", " * string(to(b)) *
                ", LB(1), " * string(wrap(b)) * "),\n"
        end
    end
    if length(bond_string) > 2
        bond_string = bond_string[1:end-2] * "\n"
    end

    # generate the complete generating code as a string
    generating_code = "" *
        "# Implementation\n" *
        "# - version " * string(version) * "\n" *
        "# - labels <: " * string(sitelabeltype) * " / " * string(bondlabeltype) * "\n" *
        "function getUnitcell" * uppercase(name[1]) * name[2:end] * "(\n" *
        "            unitcell_type :: Type{U},\n" *
        "            version       :: Val{" * string(version) * "}\n" *
        "        ) :: U where {LS<:" *
        string(sitelabeltype) * ",LB<:" * string(bondlabeltype) *
        ",S<:AbstractSite{LS," * string(D) * "},B<:AbstractBond{LB," * string(N) * "}, " *
        "U<:AbstractUnitcell{S,B}}\n" *
        "    \n" *
        "    # return a new Unitcell\n" *
        "    return newUnitcell(\n" *
        "        # Type of the unitcell\n" *
        "        U,\n" *
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
        "        ]\n" *
        "    )\n" *
        "end"

    # print the generating code
    return generating_code
end

# code for unitcell templates
function getCodeUnitcellTemplate() :: String

    # build the string of lattice vectors
    lattice_vector_string = ""
    lattice_vector_string = lattice_vector_string * "            <A1>,\n"
    lattice_vector_string = lattice_vector_string * "            <A2>,\n"
    lattice_vector_string = lattice_vector_string * "            ...\n"

    # build a string for all sites
    site_string = ""
    site_string = site_string * "            newSite(S, <POSITION>, LS(<LABEL>)),\n"
    site_string = site_string * "            newSite(S, <POSITION>, LS(<LABEL>)),\n"
    site_string = site_string * "            ...\n"

    # build a string for all bonds
    bond_string = ""
    bond_string = bond_string * "            newBond(B, <FROM>, <TO>, LB(<LABEL>), <WRAP>),\n"
    bond_string = bond_string * "            newBond(B, <FROM>, <TO>, LB(<LABEL>), <WRAP>),\n"
    bond_string = bond_string * "            ...\n"


    # generate the complete generating code as a string
    generating_code = "" *
        "# Implementation\n" *
        "# - version <VERSION>\n" *
        "# - labels <: <SITELABELTYPE> / <BONDLABELTYPE>\n" *
        "function getUnitcell<UNITCELLNAME>(\n" *
        "            unitcell_type :: Type{U},\n" *
        "            version       :: Val{<VERSION>}\n" *
        "        ) :: U where {LS<:<SITELABELTYPE>,LB<:<BONDLABELTYPE>, " *
        "S<:AbstractSite{LS,<SPACEDIMENSION>},B<:AbstractBond{LB,<LATTICEDIMENSION>}, " *
        "U<:AbstractUnitcell{S,B}}\n" *
        "    \n" *
        "    # return a new Unitcell\n" *
        "    return newUnitcell(\n" *
        "        # Type of the unitcell\n" *
        "        U,\n" *
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
        "        ]\n" *
        "    )\n" *
        "end"

    # print the generating code
    return generating_code
end

# code for unitcell version fallback
function getCodeUnitcellVersionFallback(
            unitcell    :: U,
            name        :: String = "myunitcell",
            version     :: Int64  = 1
        ) :: String where {D,LS,LB,N, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    # generate the complete generating code as a string
    generating_code = "" *
        "# Implementation\n" *
        "# - version " * string(version) * "\n" *
        "# - labels <: Any\n" *
        "# --> FALLBACK (raises error)\n" *
        "function getUnitcell" * uppercase(name[1]) * name[2:end] * "(\n" *
        "            unitcell_type :: Type{U},\n" *
        "            version       :: Val{" * string(version) * "}\n" *
        "        ) :: U where {LS,LB,S<:AbstractSite{LS," * string(D) * "},B<:AbstractBond{LB," * string(N) * "}, " *
        "U<:AbstractUnitcell{S,B}}\n" *
        "    \n" *
        "    # error since this version has no implementation yet\n" *
        "    error(\"Version " * string(version) * " of " * name * " unitcell has no implementation for label types \" * string(LS) * \" / \" * string(LB) * \" yet\")\n" *
        "end"

    # print the generating code
    return generating_code
end



# get the code for Unitcell file header
function getCodeUnitcellFileHeader(
            name :: String,
            D, N
        )

    # compile the string
    complete_code = ""

    # add the main comment to the top
    complete_code *= "################################################################################\n"
    complete_code *= "#\n"
    complete_code *= "#   " * uppercase(name) * " LATTICE\n"
    complete_code *= "#\n"
    complete_code *= "################################################################################\n\n\n\n"

    # FIRST big headline (fallback / reference)
    complete_code *= "# REFERENCE / FALLBACK (for generic AbstractUnitcell type)\n\n"

    # referencing function for abstract type
    complete_code *= "# Referencing to individual functions by wrapping in Val()\n"
    complete_code *= "function getUnitcell" * uppercase(name[1]) * name[2:end] * "(\n"
    complete_code *= "            unitcell_type   :: Type{U},\n"
    complete_code *= "            version         :: Int64 = 1\n"
    complete_code *= "        ) :: U where {LS,LB,S<:AbstractSite{LS," * string(D) * "},B<:AbstractBond{LB," * string(N) * "},U<:AbstractUnitcell{S,B}}\n"
    complete_code *= "    \n"
    complete_code *= "    # call the respective subfunction by converting to val type\n"
    complete_code *= "    return getUnitcell" * uppercase(name[1]) * name[2:end] * "(unitcell_type, Val(version))\n"
    complete_code *= "end\n\n"

    # fallback function for abstract type
    complete_code *= "# Fallback for all implementations (if Val{V} is not found)\n"
    complete_code *= "function getUnitcell" * uppercase(name[1]) * name[2:end] * "(\n"
    complete_code *= "            unitcell_type   :: Type{U},\n"
    complete_code *= "            version         :: Val{V}\n"
    complete_code *= "        ) :: U where {LS,LB,S<:AbstractSite{LS," * string(D) * "},B<:AbstractBond{LB," * string(N) * "},U<:AbstractUnitcell{S,B},V}\n"
    complete_code *= "    \n"
    complete_code *= "    # fallback / fail due to missing implementation\n"
    complete_code *= "    error(\"Version \" * string(V) * \" of " * name * " unitcell (label types \" * string(LS) * \" / \" * string(LB) * \") not implemented yet\")\n"
    complete_code *= "end\n\n\n\n"


    # SECOND big headline (wrapper for concrete type)
    complete_code *= "# WRAPPER FUNCTIONS (for concrete Unitcell type) call general function\n\n"

    # wrapper with only version (DEFAULT)
    complete_code *= "# wrapper function for passing no label types (and version) (DEFAULT)\n"
    complete_code *= "function getUnitcell" * uppercase(name[1]) * name[2:end] * "(\n"
    complete_code *= "            version :: Int64 = 1\n"
    complete_code *= "        ) :: Unitcell{Site{Int64," * string(D) * "},Bond{Int64," * string(N) * "}}\n"
    complete_code *= "    \n"
    complete_code *= "    # create a suitable unitcell of the Unitcell type\n"
    complete_code *= "    return getUnitcell" * uppercase(name[1]) * name[2:end] * "(Unitcell{Site{Int64," * string(D) * "},Bond{Int64," * string(N) * "}}, version)\n"
    complete_code *= "end\n\n"

    # wrapper with common label type and version
    complete_code *= "# wrapper function for passing common label type (and version)\n"
    complete_code *= "function getUnitcell" * uppercase(name[1]) * name[2:end] * "(\n"
    complete_code *= "            label_type :: Type{L},\n"
    complete_code *= "            version    :: Int64 = 1\n"
    complete_code *= "        ) :: Unitcell{Site{L," * string(D) * "},Bond{L," * string(N) * "}} where L\n"
    complete_code *= "    \n"
    complete_code *= "    # create a suitable unitcell of the Unitcell type\n"
    complete_code *= "    return getUnitcell" * uppercase(name[1]) * name[2:end] * "(Unitcell{Site{L," * string(D) * "},Bond{L," * string(N) * "}}, version)\n"
    complete_code *= "end\n\n"

    # wrapper with common label type and version
    complete_code *= "# wrapper function for passing site / bond label types (and version)\n"
    complete_code *= "function getUnitcell" * uppercase(name[1]) * name[2:end] * "(\n"
    complete_code *= "            label_type_site :: Type{LS},\n"
    complete_code *= "            label_type_bond :: Type{LB},\n"
    complete_code *= "            version         :: Int64 = 1\n"
    complete_code *= "        ) :: Unitcell{Site{LS," * string(D) * "},Bond{LB," * string(N) * "}} where {LS,LB}\n"
    complete_code *= "    \n"
    complete_code *= "    # create a suitable unitcell of the Unitcell type\n"
    complete_code *= "    return getUnitcell" * uppercase(name[1]) * name[2:end] * "(Unitcell{Site{LS," * string(D) * "},Bond{LB," * string(N) * "}}, version)\n"
    complete_code *= "end\n\n\n\n"



    # add the distinction to version code
    complete_code *= "################################################################################\n"
    complete_code *= "#\n"
    complete_code *= "#   VERSION IMPLEMENTATIONS FROM HERE ON\n"
    complete_code *= "#\n"
    complete_code *= "################################################################################"

    # return the code
    return complete_code
end




# print the unitcell version code
function printCodeUnitcellVersion(
            io          :: IO,
            unitcell    :: U,
            name        :: String = "myunitcell",
            version     :: Int64  = 1;
            labeltype_site  :: DataType = Nothing,
            labeltype_bond  :: DataType = Nothing
        ) where {D,LS,LB,N, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    # get the code
    code = getCodeUnitcellVersion(unitcell, name, version, labeltype_site=labeltype_site, labeltype_bond=labeltype_bond)

    # print the code
    print(io, code)
end
function printCodeUnitcellVersion(
            unitcell    :: U,
            name        :: String = "myunitcell",
            version     :: Int64  = 1;
            labeltype_site  :: DataType = Nothing,
            labeltype_bond  :: DataType = Nothing
        ) where {D,LS,LB,N, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    # get the code
    code = getCodeUnitcellVersion(unitcell, name, version, labeltype_site=labeltype_site, labeltype_bond=labeltype_bond)

    # print the code
    print(code)
end


# write a complete file for a unitcell
function writeUnitcellFile(
            unitcell    :: U,
            uc_name     :: String = "myunitcell",
            folder      :: String = "./",
            version     :: Int64  = 1
            ;
            generate_fallback :: Bool = true,
            generate_default  :: Bool = true
        ) where {D,LS,LB,N, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    # get the name for the function
    name = filter(c -> isletter(c), uc_name)
    # get the filename
    filename = folder * (folder[end] == "/" ? "" : "/") * lowercase(filter(c -> isletter(c) || c=='-', uc_name)) * ".jl"

    # mode for file opening
    openmode = "a"


    # compile the string
    complete_code = ""

    # check if the file already exists (then, leave out all header stuff)
    if !isfile(filename)

        # other mode of opening the file
        openmode = "w"

        # add the main comment to the top
        complete_code *= getCodeUnitcellFileHeader(name, D, N)

        # add space
        complete_code *= "\n\n\n"

        # DEFAULT VERSION

        # check if defaults should be generated
        if generate_default
            # check if default needs a fallback
            if generate_fallback
                # add the fallback for version 1 to the code
                complete_code *= getCodeUnitcellVersionFallback(unitcell, name, 1) * "\n\n"
            end
            # generate a default version for Number
            complete_code *= getCodeUnitcellSimplifiedVersion(unitcell, name, 1, labeltype_bond=Number, labeltype_site=Number) * "\n\n"
            # generate a default version for AbstractString
            complete_code *= getCodeUnitcellSimplifiedVersion(unitcell, name, 1, labeltype_bond=AbstractString, labeltype_site=AbstractString) * "\n\n"
        end

    end


    # CUSTOM VERSION

    # maybe some space
    if version != 1
        complete_code *= "\n\n"
    end

    # check if custom version needs a fallback
    if generate_fallback && version != 1
        # add the fallback for version 1 to the code
        complete_code *= getCodeUnitcellVersionFallback(unitcell, name, version) * "\n\n"
    end
    # generate a normal version code
    complete_code *= getCodeUnitcellVersion(unitcell, name, version) * "\n\n"





    # open the respective file and write the code into that file
    f = open(filename, openmode)
    write(f, complete_code)
    close(f)

end
