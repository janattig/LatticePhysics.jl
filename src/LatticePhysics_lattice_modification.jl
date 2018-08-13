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
#       - all connections
#       - mapping of strengths
#       - evaluate strengths
#       - optimize connections
#
#   4) REAL SPACE ROTATION / SCALING / SHIFTING
#       - Rotation around some axis
#       - Scaling along some axis
#       - Global (isotropic) scaling
#       - Global (isotropic) scaling to mean / minimum bond length
#       - Shifting along some axis
#       - TODO Shifting along some lattice vector
#       - TODO Shifting center of lattice to origin
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
        strength_return = deepcopy(conj(strength))
    else
        strength_return = deepcopy(strength)
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
        strength_return = deepcopy(conj(strength))
    else
        strength_return = deepcopy(strength)
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
export addSite!


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
export removeSite!














################################################################################
#
#   3) CONNECTION STRENGTH MODIFICATION
#       - all connections
#       - mapping of strengths
#       - evaluate strengths
#       - optimize connections
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






"""
    evaluateInteractionStrengths!(unitcell::Unitcell [, parameters::Dict])
    evaluateInteractionStrengths!(lattice::Lattice   [, parameters::Dict])

Function to evaluate connection strength values of either a `Unitcell` object or a `Lattice` object
by first optionally replacing them along keys of a dictonary `parameters` and then calling `eval(parse(.))`.

NOTE: This function modifies the given object and does not create a copy.


# Examples

Example parameters could read

    parameters = Dict("tx"=>1.0, "ty"=>1.0, "tz"=>2.0)

which modifies all `String` valued connection strengths `"tx"`, `"ty"` and `"tz"` to be `Float64` with values `1.0` and `2.0`

```julia-repl
julia> evaluateInteractionStrengths!(unitcell)

julia> evaluateInteractionStrengths!(lattice, parameters)

```
"""
function evaluateInteractionStrengths!(unitcell::Unitcell, parameters::Dict=Dict())
    # just call the mapping functions
    mapInteractionStrengths!(unitcell, parameters, replace_in_strings=true, evaluate=true)
end
function evaluateInteractionStrengths!(lattice::Lattice,   parameters::Dict=Dict())
    # just call the mapping functions
    mapInteractionStrengths!(lattice, parameters, replace_in_strings=true, evaluate=true)
end
export evaluateInteractionStrengths!






#-----------------------------------------------------------------------------------------------------------------------------
#
#   Optimize the lattice / unitcell by collapsing redundant connections into fewer connections and modify given object
#
#-----------------------------------------------------------------------------------------------------------------------------

"""
    optimizeConnections!(unitcell::Unitcell)
    optimizeConnections!(lattice::Lattice)

Function to optimize connections of either a `Unitcell` object or a `Lattice` object
by checking which ones can be compressed and put together.

NOTE: This function modifies the given object and does not create a copy.


# Examples

```julia-repl
julia> optimizeConnections!(unitcell)

```
"""
function optimizeConnections!(lattice::Lattice)
    # build up new connections
    connections_new = Array{Any,1}[]
    # go through all old connections to check for redundant entries
    for c in lattice.connections
        # check if the connection already is present in the new list
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
    connections_new_2 = Array{Any,1}[]
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
    # return nothing
    return nothing
end
function optimizeConnections!(unitcell::Unitcell)
    # build up new connections
    connections_new = Array{Any,1}[]
    # go through all old connections to check for redundant entries
    for c in unitcell.connections
        # check if the connection already is present in the new list
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
    connections_new_2 = Array{Any,1}[]
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
    # return nothing
    return nothing
end
export optimizeConnections!







################################################################################
#
#   4) REAL SPACE ROTATION / SCALING / SHIFTING
#       - Rotation around some axis
#       - Scaling along some axis
#       - Global (isotropic) scaling
#       - Global (isotropic) scaling to mean / minimum bond length
#       - Shifting along some axis
#       - TODO Shifting along some lattice vector
#       - TODO Shifting center of lattice to origin
#
################################################################################




###################
#   ROTATION
###################


# ROTATION AROUND X AXIS
function rotateAroundXAxis!(unitcell::Unitcell, angle::Float64)
    # rotate the lattice vectors
    for lv in unitcell.lattice_vectors
        # rotate and overwrite simulateously all components
        lv[2], lv[3] = cos(angle)*lv[2] + sin(angle)*lv[3]  ,   -sin(angle)*lv[2] + cos(angle)*lv[3]
    end
    # rotate the basis sites
    for bs in unitcell.basis
        # rotate and overwrite simulateously all components
        bs[2], bs[3] = cos(angle)*bs[2] + sin(angle)*bs[3]  ,   -sin(angle)*bs[2] + cos(angle)*bs[3]
    end
end
function rotateAroundXAxis!(lattice::Lattice, angle::Float64)
    # rotate the lattice vectors
    for lv in lattice.lattice_vectors
        # rotate and overwrite simulateously all components
        lv[2], lv[3] = cos(angle)*lv[2] + sin(angle)*lv[3]  ,   -sin(angle)*lv[2] + cos(angle)*lv[3]
    end
    # rotate the basis sites
    for p in lattice.positions
        # rotate and overwrite simulateously all components
        p[2], p[3] = cos(angle)*p[2] + sin(angle)*p[3]  ,   -sin(angle)*p[2] + cos(angle)*p[3]
    end
end
export rotateAroundXAxis!

function rotateAroundXAxisDeg!(unitcell::Unitcell, angle::Float64)
    rotateAroundXAxis!(unitcell, angle*pi/180.0)
end
function rotateAroundXAxisDeg!(lattice::Lattice, angle::Float64)
    rotateAroundXAxis!(lattice, angle*pi/180.0)
end
export rotateAroundXAxisDeg!


# ROTATION AROUND Y AXIS
function rotateAroundYAxis!(unitcell::Unitcell, angle::Float64)
    # rotate the lattice vectors
    for lv in unitcell.lattice_vectors
        # rotate and overwrite simulateously all components
        lv[3], lv[1] = cos(angle)*lv[3] + sin(angle)*lv[1]  ,   -sin(angle)*lv[3] + cos(angle)*lv[1]
    end
    # rotate the basis sites
    for bs in unitcell.basis
        # rotate and overwrite simulateously all components
        bs[3], bs[1] = cos(angle)*bs[3] + sin(angle)*bs[1]  ,   -sin(angle)*bs[3] + cos(angle)*bs[1]
    end
end
function rotateAroundYAxis!(lattice::Lattice, angle::Float64)
    # rotate the lattice vectors
    for lv in lattice.lattice_vectors
        # rotate and overwrite simulateously all components
        lv[3], lv[1] = cos(angle)*lv[3] + sin(angle)*lv[1]  ,   -sin(angle)*lv[3] + cos(angle)*lv[1]
    end
    # rotate the basis sites
    for p in lattice.positions
        # rotate and overwrite simulateously all components
        p[3], p[1] = cos(angle)*p[3] + sin(angle)*p[1]  ,   -sin(angle)*p[3] + cos(angle)*p[1]
    end
end
export rotateAroundYAxis!

function rotateAroundYAxisDeg!(unitcell::Unitcell, angle::Float64)
    rotateAroundYAxis!(unitcell, angle*pi/180.0)
end
function rotateAroundYAxisDeg!(lattice::Lattice, angle::Float64)
    rotateAroundYAxis!(lattice, angle*pi/180.0)
end
export rotateAroundYAxisDeg!


# ROTATION AROUND Z AXIS
function rotateAroundZAxis!(unitcell::Unitcell, angle::Float64)
    # rotate the lattice vectors
    for lv in unitcell.lattice_vectors
        # rotate and overwrite simulateously all components
        lv[1], lv[2] = cos(angle)*lv[1] + sin(angle)*lv[2]  ,   -sin(angle)*lv[1] + cos(angle)*lv[2]
    end
    # rotate the basis sites
    for bs in unitcell.basis
        # rotate and overwrite simulateously all components
        bs[1], bs[2] = cos(angle)*bs[1] + sin(angle)*bs[2]  ,   -sin(angle)*bs[1] + cos(angle)*bs[2]
    end
end
function rotateAroundZAxis!(lattice::Lattice, anglTODOe::Float64)
    # rotate the lattice vectors
    for lv in lattice.lattice_vectors
        # rotate and overwrite simulateously all components
        lv[1], lv[2] = cos(angle)*lv[1] + sin(angle)*lv[2]  ,   -sin(angle)*lv[1] + cos(angle)*lv[2]
    end
    # rotate the basis sites
    for p in lattice.positions
        # rotate and overwrite simulateously all components
        p[1], p[2] = cos(angle)*p[1] + sin(angle)*p[2]  ,   -sin(angle)*p[1] + cos(angle)*p[2]
    end
end
export rotateAroundZAxis!

function rotateAroundZAxisDeg!(unitcell::Unitcell, angle::Float64)
    rotateAroundZAxis!(unitcell, angle*pi/180.0)
end
function rotateAroundZAxisDeg!(lattice::Lattice, angle::Float64)
    rotateAroundZAxis!(lattice, angle*pi/180.0)
end
export rotateAroundZAxisDeg!





###################
#   SCALING
###################

# SCALING ALONG X AXIS
function scaleAlongXAxis!(unitcell::Unitcell, factor::Float64)
    # scale all lattice vectors
    for l in unitcell.lattice_vectors
        l[1] *= factor
    end
    # scale all positions
    for p in unitcell.basis
        p[1] *= factor
    end
end
function scaleAlongXAxis!(lattice::Lattice, factor::Float64)
    # scale all lattice vectors
    for l in lattice.lattice_vectors
        l[1] *= factor
    end
    # scale all positions
    for p in lattice.positions
        p[1] *= factor
    end
end
export scaleAlongXAxis!

# SCALING ALONG Y AXIS
function scaleAlongYAxis!(unitcell::Unitcell, factor::Float64)
    # scale all lattice vectors
    for l in unitcell.lattice_vectors
        l[2] *= factor
    end
    # scale all positions
    for p in unitcell.basis
        p[2] *= factor
    end
end
function scaleAlongYAxis!(lattice::Lattice, factor::Float64)
    # scale all lattice vectors
    for l in lattice.lattice_vectors
        l[2] *= factor
    end
    # scale all positions
    for p in lattice.positions
        p[2] *= factor
    end
end
export scaleAlongYAxis!

# SCALING ALONG Z AXIS
function scaleAlongZAxis!(unitcell::Unitcell, factor::Float64)
    # scale all lattice vectors
    for l in unitcell.lattice_vectors
        l[3] *= factor
    end
    # scale all positions
    for p in unitcell.basis
        p[3] *= factor
    end
end
function scaleAlongZAxis!(lattice::Lattice, factor::Float64)
    # scale all lattice vectors
    for l in lattice.lattice_vectors
        l[3] *= factor
    end
    # scale all positions
    for p in lattice.positions
        p[3] *= factor
    end
end
export scaleAlongZAxis!


# SCALING ISOTROPIC
function scaleIsotropic!(unitcell::Unitcell, factor::Float64)
    # scale all lattice vectors
    for l in unitcell.lattice_vectors
        l .*= factor
    end
    # scale all positions
    for p in unitcell.basis
        p .*= factor
    end
end
function scaleIsotropic!(lattice::Lattice, factor::Float64)
    # scale all lattice vectors
    for l in lattice.lattice_vectors
        l .*= factor
    end
    # scale all positions
    for p in lattice.positions
        p .*= factor
    end
end
export scaleIsotropic!


# SCALING TO MEAN BOND LENGTH
function scaleToMeanBondLength!(unitcell::Unitcell, bond_length::Float64=1.0)
    # get all bond lengths
    bond_lengths = zeros(length(unitcell.connections))
    for i in 1:length(bond_lengths)
        # get the connection
        c = unitcell.connections[i]
        # get the connecting vector
        c_vector = unitcell.basis[Int(c[2])] .- unitcell.basis[Int(c[1])]
        for (j,l) in enumerate(unitcell.lattice_vectors)
            c_vector .+= l.*c[4][j]
        end
        # get how long it is
        bond_lengths[i] = sqrt(sum(c_vector.*c_vector))
    end
    # calculate the mean bond length
    mean_bond_length = mean(bond_lengths)
    # calculate how the scale factor looks like
    factor = bond_length / mean_bond_length
    # scale the object
    scaleIsotropic!(unitcell, factor)
end
function scaleToMeanBondLength!(lattice::Lattice, bond_length::Float64=1.0)
    # get all bond lengths
    bond_lengths = zeros(length(lattice.connections))
    for i in 1:length(bond_lengths)
        # get the connection
        c = lattice.connections[i]
        # get the connecting vector
        c_vector = lattice.positions[Int(c[2])] .- lattice.positions[Int(c[1])]
        for (j,l) in enumerate(lattice.lattice_vectors)
            c_vector .+= l.*c[4][j]
        end
        # get how long it is
        bond_lengths[i] = sqrt(sum(c_vector.*c_vector))
    end
    # calculate the mean bond length
    mean_bond_length = mean(bond_lengths)
    # calculate how the scale factor looks like
    factor = bond_length / mean_bond_length
    # scale the object
    scaleIsotropic!(lattice, factor)
end
export scaleToMeanBondLength!

# SCALING TO MINIMUM BOND LENGTH
function scaleToMinimumBondLength!(unitcell::Unitcell, bond_length::Float64=1.0)
    # get all bond lengths
    bond_lengths = zeros(length(unitcell.connections))
    for i in 1:length(bond_lengths)
        # get the connection
        c = unitcell.connections[i]
        # get the connecting vector
        c_vector = unitcell.basis[Int(c[2])] .- unitcell.basis[Int(c[1])]
        for (j,l) in enumerate(unitcell.lattice_vectors)
            c_vector .+= l.*c[4][j]
        end
        # get how long it is
        bond_lengths[i] = sqrt(sum(c_vector.*c_vector))
    end
    # calculate the minimum bond length
    minimum_bond_length = minimum(bond_lengths)
    # calculate how the scale factor looks like
    factor = bond_length / minimum_bond_length
    # scale the object
    scaleIsotropic!(unitcell, factor)
end
function scaleToMinimumBondLength!(lattice::Lattice, bond_length::Float64=1.0)
    # get all bond lengths
    bond_lengths = zeros(length(lattice.connections))
    for i in 1:length(bond_lengths)
        # get the connection
        c = lattice.connections[i]
        # get the connecting vector
        c_vector = lattice.positions[Int(c[2])] .- lattice.positions[Int(c[1])]
        for (j,l) in enumerate(lattice.lattice_vectors)
            c_vector .+= l.*c[4][j]
        end
        # get how long it is
        bond_lengths[i] = sqrt(sum(c_vector.*c_vector))
    end
    # calculate the minimum bond length
    minimum_bond_length = minimum(bond_lengths)
    # calculate how the scale factor looks like
    factor = bond_length / minimum_bond_length
    # scale the object
    scaleIsotropic!(lattice, factor)
end
export scaleToMinimumBondLength!





###################
#   SHIFTING
###################


# SHIFT ALONG X AXIS
function shiftAlongXAxis!(unitcell::Unitcell, offset::Float64)
    # shift all positions
    for p in unitcell.basis
        p[1] += offset
    end
end
function shiftAlongXAxis!(lattice::Lattice, offset::Float64)
    # shift all positions
    for p in lattice.positions
        p[1] += offset
    end
end
export shiftAlongXAxis!

# SHIFT ALONG Y AXIS
function shiftAlongYAxis!(unitcell::Unitcell, offset::Float64)
    # shift all positions
    for p in unitcell.basis
        p[2] += offset
    end
end
function shiftAlongYAxis!(lattice::Lattice, offset::Float64)
    # shift all positions
    for p in lattice.positions
        p[2] += offset
    end
end
export shiftAlongYAxis!

# SHIFT ALONG Z AXIS
function shiftAlongZAxis!(unitcell::Unitcell, offset::Float64)
    # shift all positions
    for p in unitcell.basis
        p[3] += offset
    end
end
function shiftAlongZAxis!(lattice::Lattice, offset::Float64)
    # shift all positions
    for p in lattice.positions
        p[3] += offset
    end
end
export shiftAlongZAxis!


# SHIFT BY VECTOR
function shiftByVector!(unitcell::Unitcell, vector::Array{Float64,1})
    # shift all positions
    for p in unitcell.basis
        p .+= vector
    end
end
function shiftByVector!(lattice::Lattice, vector::Array{Float64,1})
    # shift all positions
    for p in lattice.positions
        p .+= vector
    end
end
export shiftByVector!


# SHIFT ALONG ANY LATTICE VECTOR
function shiftAlongLatticeVector!(unitcell::Unitcell, lattice_vector::Int64, offset::Float64)
    # define a vector for shifting
    vector = unitcell.lattice_vectors[lattice_vector] .* (offset / sqrt(sum(unitcell.lattice_vectors[lattice_vector].*unitcell.lattice_vectors[lattice_vector])))
    # shift by that vector
    shiftByVector!(unitcell, vector)
end
function shiftAlongLatticeVector!(lattice::Lattice, lattice_vector::Int64, offset::Float64)
    # define a vector for shifting
    vector = lattice.lattice_vectors[lattice_vector] .* (offset / sqrt(sum(lattice.lattice_vectors[lattice_vector].*lattice.lattice_vectors[lattice_vector])))
    # shift by that vector
    shiftByVector!(lattice, vector)
end
export shiftAlongLatticeVector!


# SHIFT ALONG A1
function shiftAlongA1!(unitcell::Unitcell, offset::Float64)
    # shift along the first lattice vector
    shiftAlongLatticeVector!(unitcell, 1, offset)
end
function shiftAlongA1!(lattice::Lattice, offset::Float64)
    # shift along the first lattice vector
    shiftAlongLatticeVector!(lattice, 1, offset)
end
export shiftAlongA1!

# SHIFT ALONG A2
function shiftAlongA2!(unitcell::Unitcell, offset::Float64)
    # shift along the second lattice vector
    shiftAlongLatticeVector!(unitcell, 2, offset)
end
function shiftAlongA2!(lattice::Lattice, offset::Float64)
    # shift along the second lattice vector
    shiftAlongLatticeVector!(lattice, 2, offset)
end
export shiftAlongA2!

# SHIFT ALONG A3
function shiftAlongA3!(unitcell::Unitcell, offset::Float64)
    # shift along the third lattice vector
    shiftAlongLatticeVector!(unitcell, 3, offset)
end
function shiftAlongA3!(lattice::Lattice, offset::Float64)
    # shift along the third lattice vector
    shiftAlongLatticeVector!(lattice, 3, offset)
end
export shiftAlongA3!


# SHIFT ALONG ANY LATTICE VECTOR (Relative)
function shiftAlongLatticeVectorRelative!(unitcell::Unitcell, lattice_vector::Int64, offset::Float64)
    # define a vector for shifting
    vector = unitcell.lattice_vectors[lattice_vector] .* offset
    # shift by that vector
    shiftByVector!(unitcell, vector)
end
function shiftAlongLatticeVectorRelative!(lattice::Lattice, lattice_vector::Int64, offset::Float64)
    # define a vector for shifting
    vector = lattice.lattice_vectors[lattice_vector] .* offset
    # shift by that vector
    shiftByVector!(lattice, vector)
end
export shiftAlongLatticeVectorRelative!


# SHIFT ALONG A1 (Relative)
function shiftAlongA1Relative!(unitcell::Unitcell, offset::Float64)
    # shift along the first lattice vector
    shiftAlongLatticeVectorRelative!(unitcell, 1, offset)
end
function shiftAlongA1Relative!(lattice::Lattice, offset::Float64)
    # shift along the first lattice vector
    shiftAlongLatticeVectorRelative!(lattice, 1, offset)
end
export shiftAlongA1Relative!

# SHIFT ALONG A2 (Relative)
function shiftAlongA2Relative!(unitcell::Unitcell, offset::Float64)
    # shift along the second lattice vector
    shiftAlongLatticeVectorRelative!(unitcell, 2, offset)
end
function shiftAlongA2Relative!(lattice::Lattice, offset::Float64)
    # shift along the second lattice vector
    shiftAlongLatticeVectorRelative!(lattice, 2, offset)
end
export shiftAlongA2Relative!

# SHIFT ALONG A3 (Relative)
function shiftAlongA3Relative!(unitcell::Unitcell, offset::Float64)
    # shift along the third lattice vector
    shiftAlongLatticeVectorRelative!(unitcell, 3, offset)
end
function shiftAlongA3Relative!(lattice::Lattice, offset::Float64)
    # shift along the third lattice vector
    shiftAlongLatticeVectorRelative!(lattice, 3, offset)
end
export shiftAlongA3Relative!
