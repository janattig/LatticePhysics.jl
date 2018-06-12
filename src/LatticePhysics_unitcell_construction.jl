################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT UNITCELL CONSTRUCTION RELATED STUFF
#
#   STRUCTURE OF THE FILE
#
#   1) UNITCELL FROM SITE COLLECTION (automatically determine the connections)
#       - 2D
#       - 3D
#       - independent of the dimension
#
#   2) UNITCELL GENERATING CODE
#
################################################################################



#-----------------------------------------------------------------------------------------------------------------------------
#
#   GET A UNITCELL OBJECT BY ONLY PASSING SITES AND LATTICE VECTORS
#
#   Parameters (required)
#   - sites: an array indicating the position of sites
#   - lattice_vectors: an array indicating the lattice vectors of the lattice
#
#   Parameters (optional)
#   - max_ij(k): the maximum number of unitcells in all direction that is beeing considered to build geometric neighbors
#   - epsilon: The precision up to which distances are considered to be equal
#   - min_NN: minimal number of nearest neighbors that ALL sites have to have (-1: Auto)
#   - max_NN: maximal number of nearest neighbors that ALL sites have to have (-1: Auto)
#   - strength_NN: interaction strength of the new connections
#   - name: name of the lattice (if "AUTO" is selected, a name will be build automatically)
#   - save: if to save the newly constructed unitcell
#
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellFromSites2D(
        sites::Array{Array, 1},
        lattice_vectors::Array{Array, 1};
        max_ij::Int64=3,
        epsilon::Float64=1e-8,
        min_NN::Int64=-1,
        max_NN::Int64=-1,
        strength_NN=1.0,
        name::String="AUTO",
        save::Bool=false
    )
    # make a new list with connections
    connections = Array[]
    # define what "close to" means
    function distance(p1, p2)
        return sqrt(sum((p1.-p2).*(p1.-p2)))
    end
    function closeto(p1, p2)
        return distance(p1, p2) < epsilon
    end
    # function to sort the checklist
    function sortfunction(checklistentry)
        return checklistentry[3]
    end
    # iterate over all sites
    for s in 1:size(sites, 1)
        # build a new checklist for sites
        checklist = Array[]
        # check all other sites in a lattice around this one
        for i in -max_ij:max_ij
        for j in -max_ij:max_ij
        for alpha in 1:size(sites, 1)
            # site at i*a1 + j*a2 + basis[alpha]
            position = i.*lattice_vectors[1] .+ j.*lattice_vectors[2] .+ sites[alpha]
            # check if the position is close to the original position
            if closeto(position, sites[s])
                # just the original site
                continue
            end
            # put the site into the list
            # format: [site_index, (wrap), distance]
            push!(checklist, [alpha, (i,j), distance(position, sites[s])])
        end
        end
        end
        # sort the checklist so that now the first element is the closest
        sort!(checklist, by=sortfunction)
        # evaluate the checklist
        # first: find out the number of nearest neighbors
        number_NN = 0
        # if unrestricted: simply go by distance and search for all sites
        # that are as close as the first (closest) neighbor
        if min_NN <= 0 && max_NN <= 0
            for c in checklist
                if closeto(c[3], checklist[1][3])
                    number_NN = number_NN + 1
                else
                    break
                end
            end
        # if restricted, search only between the min_NN^th and max_NN^th (closest) neighbors
        elseif min_NN > 0 && max_NN <=0
            number_NN = min_NN-1
            for c in checklist[min_NN:end]
                if closeto(c[3], checklist[min_NN][3])
                    number_NN = number_NN + 1
                else
                    break
                end
            end
        else
            if min_NN <= 0
                min_NN = 1
            end
            number_NN = min_NN-1
            for c in checklist[min_NN:max_NN]
                if closeto(c[3], checklist[min_NN][3])
                    number_NN = number_NN + 1
                else
                    break
                end
            end
        end
        # print the number of nearest neighbors
        println("number of nearest neighbors of site $(s): $(number_NN)")
        # check if maybe cut at the wrong element, i.e. if there are more nearest neighbors
        if closeto(checklist[number_NN][3],checklist[number_NN+1][3])
            println("cut at the wrong number of nearest neighbors, there are more that are equally close!")
        end
        # restrict the checklist to that length
        checklist = checklist[1:number_NN]
        # parse into connections
        for c in checklist
            # define a new connection of interaction strength type float or int
            if typeof(strength_NN) == Float64 || typeof(strength_NN) == Int64
                push!(connections, [s, c[1], strength_NN, c[2]])
            elseif typeof(strength_NN) == String
                # automatic connection geometry based strength
                if typeof(strength_NN) == String && strength_NN == "AUTO"
                    push!(connections, [s, c[1], 1/c[3], c[2]])
                end
            else
                println("unknown connection strength type \"$(typeof(strength_NN))\" of connection strength $(strength_NN)")
                return
            end
        end
    end
    # correct all connections - i.e. find out if the returning connection is present, if not: add it
    connections_additional = Array[]
    for c in connections
        # check if "-c" is present
        found = false
        for mc in connections
            if c[1] == mc[2] && c[2] == mc[1] && c[3] == mc[3] && -c[4][1] == mc[4][1] && -c[4][2] == mc[4][2]
                found = true
            end
        end
        if !found
            push!(connections_additional, [c[2], c[1], c[3], (-c[4][1], -c[4][2])])
        end
    end
    # add all additional connections
    for c in connections_additional
        push!(connections, c)
    end
    # build a filename
    if name == "AUTO"
        name = "from_$(size(sites, 1))_sites"
    end
    filename = "$(FOLDER_UNITCELLS)2d_$(name)_unitcell.jld"
    # build the finished unitcell object
    unitcell = Unitcell(
        lattice_vectors,
        sites,
        connections,
        filename
    )
    # maybe save
    if save
        saveUnitcell(unitcell)
    end
    # return the unitcell
    return unitcell
end
export getUnitcellFromSites2D

function getUnitcellFromSites3D(
        sites::Array{Array, 1},
        lattice_vectors::Array{Array, 1};
        max_ij::Int64=3,
        epsilon::Float64=1e-8,
        min_NN::Int64=-1,
        max_NN::Int64=-1,
        strength_NN=1.0,
        name::String="AUTO",
        save::Bool=false
    )
    # make a new list with connections
    connections = Array[]
    # define what "close to" means
    function distance(p1, p2)
        return sqrt(sum((p1.-p2).*(p1.-p2)))
    end
    function closeto(p1, p2)
        return distance(p1, p2) < epsilon
    end
    # function to sort the checklist
    function sortfunction(checklistentry)
        return checklistentry[3]
    end
    # iterate over all sites
    for s in 1:size(sites, 1)
        # build a new checklist for sites
        checklist = Array[]
        # check all other sites in a lattice around this one
        for i in -max_ijk:max_ijk
        for j in -max_ijk:max_ijk
        for k in -max_ijk:max_ijk
        for alpha in 1:size(sites, 1)
            # site at i*a1 + j*a2 + k*a3 + basis[alpha]
            position = i.*lattice_vectors[1] .+ j.*lattice_vectors[2] .+ k.*lattice_vectors[3] .+ sites[alpha]
            # check if the position is close to the original position
            if closeto(position, sites[s])
                # just the original site
                continue
            end
            # put the site into the list
            # format: [site_index, (wrap), distance]
            push!(checklist, [alpha, (i,j,k), distance(position, sites[s])])
        end
        end
        end
        end
        # sort the checklist
        sort!(checklist, by=sortfunction)
        # evaluate the checklist
        # first: find out the number of nearest neighbors
        number_NN = 0
        # if unrestricted: simply go by distance and search for all sites
        # that are as close as the first (closest) neighbor
        if min_NN <= 0 && max_NN <= 0
            for c in checklist
                if closeto(c[3], checklist[1][3])
                    number_NN = number_NN + 1
                else
                    break
                end
            end
        # if restricted, search only between the min_NN^th and max_NN^th (closest) neighbors
        elseif min_NN > 0 && max_NN <=0
            number_NN = min_NN-1
            for c in checklist[min_NN:end]
                if closeto(c[3], checklist[min_NN][3])
                    number_NN = number_NN + 1
                else
                    break
                end
            end
        else
            if min_NN <= 0
                min_NN = 1
            end
            number_NN = min_NN-1
            for c in checklist[min_NN:max_NN]
                if closeto(c[3], checklist[min_NN][3])
                    number_NN = number_NN + 1
                else
                    break
                end
            end
        end
        # print the number of nearest neighbors
        println("number of nearest neighbors of site $(s): $(number_NN)")
        # check if maybe cut at the wrong element, i.e. if there are more nearest neighbors
        if closeto(checklist[number_NN][3],checklist[number_NN+1][3])
            println("cut at the wrong number of nearest neighbors, there are more that are equally close!")
        end
        # restrict the checklist to that length
        checklist = checklist[1:number_NN]
        # parse into connections
        for c in checklist
            # define a new connection of interaction strength type float or int
            if typeof(strength_NN) == Float64 || typeof(strength_NN) == Int64
                push!(connections, [s, c[1], strength_NN, c[2]])
            elseif typeof(strength_NN) == String
                # automatic connection geometry based strength
                if typeof(strength_NN) == String && strength_NN == "AUTO"
                    push!(connections, [s, c[1], 1/c[3], c[2]])
                end
            else
                println("unknown connection strength type \"$(typeof(strength_NN))\" of connection strength $(strength_NN)")
                return
            end
        end
    end
    # correct all connections - i.e. find out if the returning connection is present, if not: add it
    connections_additional = Array[]
    for c in connections
        # check if "-c" is present
        found = false
        for mc in connections
            if c[1] == mc[2] && c[2] == mc[1] && c[3] == mc[3] && -c[4][1] == mc[4][1] && -c[4][2] == mc[4][2] && -c[4][3] == mc[4][3]
                found = true
            end
        end
        if !found
            push!(connections_additional, [c[2], c[1], c[3], (-c[4][1], -c[4][2], -c[4][3])])
        end
    end
    # add all additional connections
    for c in connections_additional
        push!(connections, c)
    end
    # build a filename
    if name == "AUTO"
        name = "from_$(size(sites, 1))_sites"
    end
    filename = "$(FOLDER_UNITCELLS)3d_$(name)_unitcell.jld"
    # build the finished unitcell object
    unitcell = Unitcell(
        lattice_vectors,
        sites,
        connections,
        filename
    )
    # maybe save
    if save
        saveUnitcell(unitcell)
    end
    # return the unitcell
    return unitcell
end
export getUnitcellFromSites3D

function getUnitcellFromSites(
        sites::Array{Array, 1},
        lattice_vectors::Array{Array, 1};
        max_ij::Int64=3,
        epsilon::Float64=1e-8,
        min_NN::Int64=-1,
        max_NN::Int64=-1,
        strength_NN=1.0,
        name::String="AUTO",
        save::Bool=false
    )
    # check the dimension and call the respective function
    if size(lattice_vectors, 1) == 1
        println("building from sites with one lattice vector not supported")
    elseif size(lattice_vectors, 1) == 2
        return getUnitcellFromSites2D(
            sites,
            lattice_vectors;
            max_ij=max_ijk,
            epsilon=epsilon,
            min_NN=min_NN,
            max_NN=max_NN,
            strength_NN=strength_NN,
            name=name,
            save=save
        )
    elseif size(lattice_vectors, 1) == 3
        return getUnitcellFromSites3D(
            sites,
            lattice_vectors;
            max_ijk=max_ijk,
            epsilon=epsilon,
            min_NN=min_NN,
            max_NN=max_NN,
            strength_NN=strength_NN,
            name=name,
            save=save
        )
    elseif size(lattice_vectors, 1) > 3
        println("building from sites with more than 3 lattice vector not supported")
    end
end
export getUnitcellFromSites









"""
    printUnitcellGeneratingCode(unitcell::Unitcell [; version::Int64=1])

Function that prints julia code into the console which can be copied and used from there on.
The printed code is a function definiton for construction of the `Unitcell` object passed.

If the additional `version` parameter is passed as well, the code will be versionized in
the convention that is also used through the library.


# Examples

```julia-repl
julia> printUnitcellGeneratingCode(unitcell)
...

julia> printUnitcellGeneratingCode(unitcell, version=4)
...
```
"""
function printUnitcellGeneratingCode(unitcell::Unitcell; version::Int64=1)
    # distinguish if versioning yes or no
    if version != -1
        # header printing
        println("function getUnitcellCustom(version::Int64=$(version); save::Bool=false)")
        # version distinguish
        println("   # distinguish by version")
        println("   if version == $(version)")
        # lattice vectors
        println("       # the lattice vectors")
        index = 0
        for l in unitcell.lattice_vectors
            index = index+1
            println("       a$(index) = $(l)")
        end
        println("       lattice_vectors = Array[]")
        for index in 1:size(unitcell.lattice_vectors, 1)
            println("       push!(lattice_vectors, a$(index))")
        end
        # Basis
        println("       # Basis Definition")
        println("       basis = Array[")
        for b in unitcell.basis[1:end-1]
            println("           $(b),")
        end
        println("           $(unitcell.basis[end])")
        println("       ]")
        # Connections
        println("       # Connection Definition")
        println("       # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]")
        println("       connections = Array[")
        for c in unitcell.connections[1:end-1]
            if typeof(c[3]) == String
                println("           [$(c[1]); $(c[2]); \"$(c[3])\"; $(c[4])],")
            else
                println("           [$(c[1]); $(c[2]); $(c[3]); $(c[4])],")
            end
        end
        c = unitcell.connections[end]
        if typeof(c[3]) == String
            println("           [$(c[1]); $(c[2]); \"$(c[3])\"; $(c[4])],")
        else
            println("           [$(c[1]); $(c[2]); $(c[3]); $(c[4])],")
        end
        println("       ]")
        # Filename
        println("       # filename")
        if size(unitcell.lattice_vectors,1) == 2
            println("       filename = \"\$(FOLDER_UNITCELLS)2d_custom_unitcell.jld\"")
        elseif size(unitcell.lattice_vectors,1) == 3
            println("       filename = \"\$(FOLDER_UNITCELLS)3d_custom_unitcell.jld\"")
        else
            println("       filename = \"\$(FOLDER_UNITCELLS)Xd_custom_unitcell.jld\"")
        end
        # end
        println("   end")
    else
        # header printing
        println("function getUnitcellCustom(version::Int64=1; save::Bool=false)")
        # lattice vectors
        println("   # the lattice vectors")
        index = 0
        for l in unitcell.lattice_vectors
            index = index+1
            println("   a$(index) = $(l)")
        end
        println("   lattice_vectors = Array[]")
        for index in 1:size(unitcell.lattice_vectors, 1)
            println("   push!(lattice_vectors, a$(index))")
        end
        # Basis
        println("   # Basis Definition")
        println("   basis = Array[")
        for b in unitcell.basis[1:end-1]
            println("       $(b),")
        end
        println("       $(unitcell.basis[end])")
        println("   ]")
        # Connections
        println("   # Connection Definition")
        println("   # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]")
        println("   connections = Array[")
        for c in unitcell.connections[1:end-1]
            if typeof(c[3]) == String
                println("       [$(c[1]); $(c[2]); \"$(c[3])\"; $(c[4])],")
            else
                println("       [$(c[1]); $(c[2]); $(c[3]); $(c[4])],")
            end
        end
        c = unitcell.connections[end]
        if typeof(c[3]) == String
            println("       [$(c[1]); $(c[2]); \"$(c[3])\"; $(c[4])],")
        else
            println("       [$(c[1]); $(c[2]); $(c[3]); $(c[4])],")
        end
        println("   ]")
        # Filename
        println("   # filename")
        if size(unitcell.lattice_vectors,1) == 2
            println("   filename = \"\$(FOLDER_UNITCELLS)2d_custom_unitcell.jld\"")
        elseif size(unitcell.lattice_vectors,1) == 3
            println("   filename = \"\$(FOLDER_UNITCELLS)3d_custom_unitcell.jld\"")
        else
            println("   filename = \"\$(FOLDER_UNITCELLS)Xd_custom_unitcell.jld\"")
        end
    end
    # part which is not version dependent
    println("   # generate unitcell")
    println("   uc = Unitcell(lattice_vectors, basis, connections, filename)")
    println("   if save")
    println("       saveUnitcell(uc)")
    println("   end")
    println("   # return the unitcell")
    println("   return uc")
    println("end")
end
export printUnitcellGeneratingCode
