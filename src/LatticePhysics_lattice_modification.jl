################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT LATTICE MODIFICATION FUNCTIONS
#   (MOST OF THEM ALSO WORK FOR UNITCELLS)
#
#   STRUCTURE OF THE FILE
#
#   1) CONNECTIONS ADDING REMOVING
#       - Add a new connection
#       - remove connections (based on indices)
#       - remove connections (based on strength)
#
#   2) SITES ADDING REMOVING
#       - add a new site
#       - remove a site (based on index)
#
#   3) CONNECTION STRENGTH MODIFICATION
#       - TODO all connections
#       - TODO mapping of strengths
#       - TODO optimize connections
#
################################################################################
















################################################################################
#
#   1) CONNECTIONS ADDING REMOVING
#       - Add a new connection
#       - remove connections (based on indices)
#       - remove connections (based on strength)
#
################################################################################



#-------------------------------------------------------------------------------
#
#  Add connection to given unitcell / lattice objects
#
#-------------------------------------------------------------------------------
"""
    addConnection!(unitcell::Unitcell, index_from::Int64, index_to::Int64, strength [, wrap=-1; overwrite::Bool=false])
    addConnection!(lattice::Lattice,   index_from::Int64, index_to::Int64, strength [, wrap=-1; overwrite::Bool=false])

Function to add a connection to either a `Unitcell` object or a `Lattice` object.
Necessary parameters include the indices of sites between which the connection takes place as well as the strength of the connection.
Additionally (in case of a unitcell or periodic boundaries), a wrap can be specified in terms of lattice vectors.
If no wrap is specified, the wrap will be chosen automatically as without wrapping around periodic boundaries.

The additional option `overwrite` specifies if one wants to simply add the connection
or if one wants to first search if there is a connection that already satisfies the new connection (in which case nothing would be added).

NOTE 1: This function always tries to add two connections, one going forward, one going backward.

NOTE 2: This function modifies the given object and does not create a copy.


# Examples

```julia-repl
julia> addConnection!(unitcell, 1, 3, "tx", (0,1))


julia> addConnection!(lattice, 31, 27, "tx")

```
"""
function addConnection!(unitcell::Unitcell, index_from::Int64, index_to::Int64, strength, wrap=-1; overwrite::Bool=false)

    # maybe default wrap
    if wrap == -1
        if length(unitcell.lattice_vectors) == 3
            wrap = (0,0,0)
        elseif length(unitcell.lattice_vectors) == 2
            wrap = (0,0)
        else
            wrap = (0)
        end
    end

    # construct the returning wrap
    if length(wrap) == 3
        wrap_return = (-wrap[1],-wrap[2],-wrap[3])
    elseif length(wrap) == 2
        wrap_return = (-wrap[1],-wrap[2])
    else
        wrap_return = (-wrap)
    end
    # construct the returning strength
    if typeof(strength) == Complex
        strength_return = copy(conj(strength))
    else
        strength_return = copy(strength)
    end

    # construct the connection as an array
    connection_1 = Any[index_from; index_to; strength; wrap]
    # construct the returning connection as an array
    connection_2 = Any[index_to; index_from; strength_return; wrap_return]

    # if not overwrite, check if the connection exists
    if !overwrite
        # report if found
        found_c1 = false
        found_c2 = false
        # check in all connections
        for c in unitcell.connections
            if  Int(c[1]) == Int(connection_1[1]) && Int(c[2]) == Int(connection_1[2]) && c[3] == connection_1[3] && c[4] == connection_1[4]
                found_c1 = true
            end
            if  Int(c[1]) == Int(connection_2[1]) && Int(c[2]) == Int(connection_2[2]) && c[3] == connection_2[3] && c[4] == connection_2[4]
                found_c2 = true
            end
        end
        # only add if not added already
        if !found_c1
            push!(unitcell.connections, connection_1);
        end
        if !found_c2
            push!(unitcell.connections, connection_2);
        end
    # otherwise, just add the connections
    else
        push!(unitcell.connections, connection_1);
        push!(unitcell.connections, connection_2);
    end

    # return nothing
    return nothing
end
function addConnection!(lattice::Lattice, index_from::Int64, index_to::Int64, strength, wrap=-1; overwrite::Bool=false)

    # maybe default wrap
    if wrap == -1
        if length(lattice.lattice_vectors) == 3
            wrap = (0,0,0)
        elseif length(lattice.lattice_vectors) == 2
            wrap = (0,0)
        else
            wrap = (0)
        end
    end

    # construct the returning wrap
    if length(wrap) == 3
        wrap_return = (-wrap[1],-wrap[2],-wrap[3])
    elseif length(wrap) == 2
        wrap_return = (-wrap[1],-wrap[2])
    else
        wrap_return = (-wrap)
    end
    # construct the returning strength
    if typeof(strength) == Complex
        strength_return = copy(conj(strength))
    else
        strength_return = copy(strength)
    end

    # construct the connection as an array
    connection_1 = Any[index_from; index_to; strength; wrap]
    # construct the returning connection as an array
    connection_2 = Any[index_to; index_from; strength_return; wrap_return]

    # if not overwrite, check if the connection exists
    if !overwrite
        # report if found
        found_c1 = false
        found_c2 = false
        # check in all connections
        for c in lattice.connections
            if  Int(c[1]) == Int(connection_1[1]) && Int(c[2]) == Int(connection_1[2]) && c[3] == connection_1[3] && c[4] == connection_1[4]
                found_c1 = true
            end
            if  Int(c[1]) == Int(connection_2[1]) && Int(c[2]) == Int(connection_2[2]) && c[3] == connection_2[3] && c[4] == connection_2[4]
                found_c2 = true
            end
        end
        # only add if not added already
        if !found_c1
            push!(lattice.connections, connection_1);
        end
        if !found_c2
            push!(lattice.connections, connection_2);
        end
    # otherwise, just add the connections
    else
        push!(lattice.connections, connection_1);
        push!(lattice.connections, connection_2);
    end

    # return nothing
    return nothing
end
export addConnection!






#-------------------------------------------------------------------------------
#
#   Remove interactions with given strengths of the unitcell / lattice
#
#-------------------------------------------------------------------------------
"""
    removeConnections!(unitcell::Unitcell, index_from::Int64, index_to::Int64 [, strength ])
    removeConnections!(lattice::Lattice,   index_from::Int64, index_to::Int64 [, strength ])

Function to remove connections from either a `Unitcell` object or a `Lattice` object based on the sites that they connect.
Optionally, one can pass a connection strength to also filter based on this strength. If no strength is passed, all connnections
between the two sites will be removed.

    removeConnections!(unitcell::Unitcell  [, strength])
    removeConnections!(lattice::Lattice    [, strength])

Function to remove connections from either a `Unitcell` object or a `Lattice` object only based on the strength of connections.
If no strength is passed, all connnections with Float64 strength `0.0` will be removed.

NOTE 1: This function always tries to delte both connections, one going forward, one going backward.

NOTE 2: This function modifies the given object and does not create a copy.


# Examples

```julia-repl
julia> removeConnections!(unitcell, 1, 3, "tx")

julia> removeConnections!(unitcell, 1, 3)

julia> removeConnections!(lattice, 23, 73)

julia> removeConnections!(lattice, "tx")

```
"""
function removeConnections!(lattice::Lattice,   index_from::Int64, index_to::Int64, strength)
    # use the julia filter function to filter out all elements which have the strength passed and the fitting indices
    filter!(x->!(x[3]==strength && Int(x[1])==index_from && Int(x[2])==index_to) && !(x[3]==strength && Int(x[2])==index_from && Int(x[1])==index_to), lattice.connections)
    # return nothing
    return nothing
end
function removeConnections!(unitcell::Unitcell, index_from::Int64, index_to::Int64, strength)
    # use the julia filter function to filter out all elements which have the strength passed and the fitting indices
    filter!(x->!(x[3]==strength && Int(x[1])==index_from && Int(x[2])==index_to) && !(x[3]==strength && Int(x[2])==index_from && Int(x[1])==index_to), unitcell.connections)
    # return nothing
    return nothing
end
function removeConnections!(lattice::Lattice,   index_from::Int64, index_to::Int64)
    # use the julia filter function to filter out all elements which have the fitting indices
    filter!(x->!(Int(x[1])==index_from && Int(x[2])==index_to) && !(Int(x[2])==index_from && Int(x[1])==index_to), lattice.connections)
    # return nothing
    return nothing
end
function removeConnections!(unitcell::Unitcell, index_from::Int64, index_to::Int64)
    # use the julia filter function to filter out all elements which have the fitting indices
    filter!(x->!(Int(x[1])==index_from && Int(x[2])==index_to) && !(Int(x[2])==index_from && Int(x[1])==index_to), unitcell.connections)
    # return nothing
    return nothing
end
function removeConnections!(lattice::Lattice,   strength=0.0)
    # use the julia filter function to filter out all elements which have the strength passed
    filter!(x->x[3]!=strength, lattice.connections)
    # return nothing
    return nothing
end
function removeConnections!(unitcell::Unitcell, strength=0.0)
    # use the julia filter function to filter out all elements which have the strength passed
    filter!(x->x[3]!=strength, unitcell.connections)
    # return nothing
    return nothing
end

export removeConnections!











################################################################################
#
#   2) SITES ADDING REMOVING
#       - add a new site
#       - remove a site (based on index)
#
################################################################################

"""
    addSite!(unitcell::Unitcell, position::Array{Float64,1})
    addSite!(lattice::Lattice,   position::Array{Float64,1})

Function to add a site to either a `Unitcell` object or a `Lattice` object which is located at
the real space position `position`. The function returns the index of the newly added site for further use.


NOTE 1: The newly added site is not connected to any other sites yet.

NOTE 2: This function modifies the given object and does not create a copy.


# Examples

```julia-repl
julia> addSite!(unitcell, Float64[0.0, 0.5])
5

julia> addSite!(lattice, Float64[0.0, 0.5, 1.0])
187
```
"""
function addSite!(lattice::Lattice, position::Array{Float64,1})
    # push the site into the positions list
    push!(lattice.positions, position)
    # return the site index
    return length(lattice.positions)
end
function addSite!(unitcell::Unitcell, position::Array{Float64,1})
    # push the site into the positions list
    push!(unitcell.basis, position)
    # return the site index
    return length(unitcell.basis)
end


"""
    removeSite!(unitcell::Unitcell, index::Int64)
    removeSite!(lattice::Lattice,   index::Int64)

Function to remove a site from either a `Unitcell` object or a `Lattice` object with index `index`.
The function furhtermore removes all connections that are connected to the site and relabels all
other connections as well to account for the index shift introduced by the removal.

NOTE: This function modifies the given object and does not create a copy.


# Examples

```julia-repl
julia> removeSite!(unitcell, 2)

julia> removeSite!(lattice, 128)

```
"""
function removeSite!(lattice::Lattice, index::Int64)
    # delete the element from the positions list
    deleteat!(lattice.positions, index)
    # delete all connections that touch the removed site
    deleteat!(lattice.connections, [c[1]==index || c[2]==index for c in lattice.connections])
    # now shift all indices in the connections array by 1
    for (i,c) in enumerate(lattice.connections)
        if Int(c[1]) > index
            c[1] = Int(c[1]) - 1
        end
        if Int(c[2]) > index
            c[2] = Int(c[2]) - 1
        end
    end
    # return nothing
    return nothing
end
function removeSite!(unitcell::Unitcell, index::Int64)
    # delete the element from the positions list
    deleteat!(unitcell.basis, index)
    # delete all connections that touch the removed site
    deleteat!(unitcell.connections, [c[1]==index || c[2]==index for c in unitcell.connections])
    # now shift all indices in the connections array by 1
    for (i,c) in enumerate(unitcell.connections)
        if Int(c[1]) > index
            c[1] = Int(c[1]) - 1
        end
        if Int(c[2]) > index
            c[2] = Int(c[2]) - 1
        end
    end
    # return nothing
    return nothing
end















################################################################################
#
#   3) CONNECTION STRENGTH MODIFICATION
#       - all connections
#       - mapping of strengths
#       - TODO optimize connections
#
################################################################################





#-----------------------------------------------------------------------------------------------------------------------------
#
#   Set / Map interaction strengths of the unitcell / lattice
#
#-----------------------------------------------------------------------------------------------------------------------------
"""
    setAllInteractionStrengths!(unitcell::Unitcell, strength)
    setAllInteractionStrengths!(lattice::Lattice,   strength)

Function to overwrite ALL connection strength values of either a `Unitcell` object or a `Lattice` object
with the value `strength`.

NOTE: This function modifies the given object and does not create a copy.


# Examples

```julia-repl
julia> setAllInteractionStrengths!(unitcell, "tx")

julia> setAllInteractionStrengths!(lattice, 1.0)

```
"""
function setAllInteractionStrengths!(unitcell::Unitcell, strength)
    # go through all connections and modify their strength
    for c in unitcell.connections
        c[3] = strength
    end
    # return nothing
    return nothing
end
function setAllInteractionStrengths!(lattice::Lattice,   strength)
    # go through all connections and modify their strength
    for c in lattice.connections
        c[3] = strength
    end
    # return nothing
    return nothing
end
export setAllInteractionStrengths!





"""
    mapInteractionStrengths!(unitcell::Unitcell, mapping::Dict [; replace_in_strings::Bool=false, evaluate::Bool=false])
    mapInteractionStrengths!(lattice::Lattice,   mapping::Dict [; replace_in_strings::Bool=false, evaluate::Bool=false])

Function to map connection strength values of either a `Unitcell` object or a `Lattice` object
by replacing them along keys of a dictonary `mapping`.

Furthermore, to be able to resolve strengths like `"tx*ty"` (`String` valued) to a `Float64` value one can use the following two steps:
1) chosing `replace_in_strings=true` to also search for mapping keys in `String` connection strengths.
2) chosing `evaluate=true` to evaluate every `String` connection strength by `eval(parse(.))` after all replacement has taken place.

NOTE: This function modifies the given object and does not create a copy.


# Examples

Example mapping could read

    mapping = Dict("tx"=>1.0, "ty"=>1.0, "tz"=>2.0)

which modifies all `String` valued connection strengths `"tx"`, `"ty"` and `"tz"` to be `Float64` with values `1.0` and `2.0`

```julia-repl
julia> mapInteractionStrengths!(unitcell, mapping)

julia> mapInteractionStrengths!(lattice, mapping, replace_in_strings=true)

```
"""
function mapInteractionStrengths!(unitcell::Unitcell, mapping::Dict; replace_in_strings::Bool=false, evaluate::Bool=false)
    # get the old strengths
    old_strengths = keys(mapping)
    # iterate for replacement
    for old_strength in old_strengths
        # new strength is given by mapping
        new_strength = mapping[old_strength]
        # check all connections that are identical
        for c in unitcell.connections
            if c[3] == old_strength
                c[3] = new_strength
            end
        end
        # check for string replacement
        if replace_in_strings
            for c in unitcell.connections
                if typeof(c[3]) == String && contains(c[3], old_strength)
                    c[3] = replace(c[3], old_strength, new_strength)
                end
            end
        end
    end
    # maybe even evaluate
    if evaluate
        for c in unitcell.connections
            if typeof(c[3]) == String
                c[3] = eval(parse(c[3]))
            end
        end
    end
    # return nothing
    return nothing
end
function mapInteractionStrengths!(lattice::Lattice,   mapping::Dict; replace_in_strings::Bool=false, evaluate::Bool=false)
    # get the old strengths
    old_strengths = keys(mapping)
    # iterate for replacement
    for old_strength in old_strengths
        # new strength is given by mapping
        new_strength = mapping[old_strength]
        # check all connections that are identical
        for c in lattice.connections
            if c[3] == old_strength
                c[3] = new_strength
            end
        end
        # check for string replacement
        if replace_in_strings
            for c in lattice.connections
                if typeof(c[3]) == String && contains(c[3], old_strength)
                    c[3] = replace(c[3], old_strength, new_strength)
                end
            end
        end
    end
    # maybe even evaluate
    if evaluate
        for c in lattice.connections
            if typeof(c[3]) == String
                c[3] = eval(parse(c[3]))
            end
        end
    end
    # return nothing
    return nothing
end
export mapInteractionStrengths!





#-----------------------------------------------------------------------------------------------------------------------------
#
#   Optimize the lattice / unitcell by collapsing redundant connections into fewer connections and modify given object
#
#-----------------------------------------------------------------------------------------------------------------------------
function optimizeConnections!(lattice::Lattice)
    # build up new connections
    connections_new = Array[]
    # go through all old connections
    for c in lattice.connections
        # check if the connection already is present in the list
        found = false
        for c_new in connections_new
            if c[1] == c_new[1] && c[2] == c_new[2] && c[4] == c_new[4]
                # found a redundent connection
                found = true
                # add connection strength
                if typeof(c[3]) == String || typeof(c_new[3]) == String
                    c_new[3] = "$(c_new[3])+$(c[3])"
                else
                    c_new[3] = c_new[3] + c[3]
                end
                # break the inner loop
                break
            end
        end
        # if not, add the current connection
        if !found
            push!(connections_new, c)
        end
    end
    # check if the connection strength is 0
    connections_new_2 = Array[]
    for c in connections_new
        # check if close to 0 connection strength
        if typeof(c[3]) == Float64 || typeof(c[3]) == Int64
            if abs(c[3]) > 1e-18
                push!(connections_new_2, c)
            else
                # do not push! c[3] too small!
            end
        elseif isnull(tryparse(Float64,c[3]))
            # not a float and cannot be parsed to float, just push in list as it is
            push!(connections_new_2, c)
        else
            # check if the strength is non-zero
            if abs(parse(Float64,c[3])) > 1e-18
                push!(connections_new_2, c)
            else
                # do not push! c[3] too small!
            end
        end
    end
    # overwrite the connections in the given lattice
    lattice.connections = connections_new_2
end
function optimizeConnections!(unitcell::Unitcell)
    # build up new connections
    connections_new = Array[]
    # go through all old connections
    for c in unitcell.connections
        # check if the connection already is present in the list
        found = false
        for c_new in connections_new
            if c[1] == c_new[1] && c[2] == c_new[2] && c[4] == c_new[4]
                # found a redundent connection
                found = true
                # add connection strength
                if typeof(c[3]) == String || typeof(c_new[3]) == String
                    c_new[3] = "$(c_new[3])+$(c[3])"
                else
                    c_new[3] = c_new[3] + c[3]
                end
                # break the inner loop
                break
            end
        end
        # if not, add the current connection
        if !found
            push!(connections_new, c)
        end
    end
    # check if the connection strength is 0
    connections_new_2 = Array[]
    for c in connections_new
        # check if close to 0 connection strength
        if typeof(c[3]) == Float64 || typeof(c[3]) == Int64
            if abs(c[3]) > 1e-18
                push!(connections_new_2, c)
            else
                # do not push! c[3] too small!
            end
        elseif isnull(tryparse(Float64,c[3]))
            # not a float and cannot be parsed to float, just push in list as it is
            push!(connections_new_2, c)
        else
            # check if the strength is non-zero
            if abs(parse(Float64,c[3])) > 1e-18
                push!(connections_new_2, c)
            else
                # do not push! c[3] too small!
            end
        end
    end
    # overwrite the connections in the given lattice
    unitcell.connections = connections_new_2
end
export optimizeConnections!






#-----------------------------------------------------------------------------------------------------------------------------
#
#   Add interactions along next-nearest neighbors
#
#-----------------------------------------------------------------------------------------------------------------------------
function addNextNearestNeighborsToConnections!(lattice::Lattice; strengthNNN="AUTO", connectionsNN=["ALL"])
    # get the connectivity matrix of the lattice
    connectivity = getConnectionList(lattice)
    # define a list of new connections
    connections_new = Array[]
    # iterate over all sites and check if they host a NNN connection
    for i in 1:size(lattice.positions,1)
        # check the neighbors of site i
        for (i1,c1) in enumerate(connectivity[i])
        for (i2,c2) in enumerate(connectivity[i])
            # check if the connections are valid
            if !(in("ALL",connectionsNN) || (in(c1[3],connectionsNN) && in(c2[3],connectionsNN)))
                continue
            end
            # propose a connection c1+c2 if they are different
            if i1 < i2
                # build the offsets
                if length(c1[4]) == 1
                    off_1 = (c2[4] - c1[4])
                    off_2 = (c1[4] - c2[4])
                elseif length(c1[4]) == 2
                    off_1 = (c2[4][1] - c1[4][1], c2[4][2] - c1[4][2])
                    off_2 = (c1[4][1] - c2[4][1], c1[4][2] - c2[4][2])
                elseif length(c1[4]) == 3
                    off_1 = (c2[4][1] - c1[4][1], c2[4][2] - c1[4][2], c2[4][3] - c1[4][3])
                    off_2 = (c1[4][1] - c2[4][1], c1[4][2] - c2[4][2], c1[4][3] - c2[4][3])
                end
                # build two new connections
                if strengthNNN == "AUTO"
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); "NNN$(c1[2])to$(c2[2])"; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); "NNN$(c1[2])to$(c2[2])"; off_2]
                elseif strengthNNN == "MULTIPLY"
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); "$(c1[3])*$(c2[3])"; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); "$(c1[3])*$(c2[3])"; off_2]
                else
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); strengthNNN; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); strengthNNN; off_2]
                end
                # push them to the add list
                push!(connections_new, connection_new_1)
                push!(connections_new, connection_new_2)
            end
        end
        end
    end
    # push all connections to the lattice
    for c in connections_new
        push!(lattice.connections, c)
    end
end
function addNextNearestNeighborsToConnections!(unitcell::Unitcell; strengthNNN="AUTO", connectionsNN=["ALL"])
    # get the connectivity matrix of the unitcell
    connectivity = getConnectionList(unitcell)
    # define a list of new connections
    connections_new = Array[]
    # iterate over all sites and check if they host a NNN connection
    for i in 1:size(unitcell.basis,1)
        # check the neighbors of site i
        for (i1,c1) in enumerate(connectivity[i])
        for (i2,c2) in enumerate(connectivity[i])
            # check if the connections are valid
            if !(in("ALL",connectionsNN) || (in(c1[3],connectionsNN) && in(c2[3],connectionsNN)))
                continue
            end
            # propose a connection c1+c2 if they are different
            if i1 < i2
                # build the offsets
                if length(c1[4]) == 1
                    off_1 = (c2[4] - c1[4])
                    off_2 = (c1[4] - c2[4])
                elseif length(c1[4]) == 2
                    off_1 = (c2[4][1] - c1[4][1], c2[4][2] - c1[4][2])
                    off_2 = (c1[4][1] - c2[4][1], c1[4][2] - c2[4][2])
                elseif length(c1[4]) == 3
                    off_1 = (c2[4][1] - c1[4][1], c2[4][2] - c1[4][2], c2[4][3] - c1[4][3])
                    off_2 = (c1[4][1] - c2[4][1], c1[4][2] - c2[4][2], c1[4][3] - c2[4][3])
                end
                # build two new connections
                if strengthNNN == "AUTO"
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); "NNN$(c1[2])to$(c2[2])"; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); "NNN$(c1[2])to$(c2[2])"; off_2]
                elseif strengthNNN == "MULTIPLY"
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); "$(c1[3])*$(c2[3])"; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); "$(c1[3])*$(c2[3])"; off_2]
                else
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); strengthNNN; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); strengthNNN; off_2]
                end
                # push them to the add list
                push!(connections_new, connection_new_1)
                push!(connections_new, connection_new_2)
            end
        end
        end
    end
    # push all connections to the unitcell
    for c in connections_new
        push!(unitcell.connections, c)
    end
end

export addNextNearestNeighborsToConnections!
