

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
#   Set / Map interaction strengths of the unitcell / lattice
#
#-----------------------------------------------------------------------------------------------------------------------------
function setAllInteractionStrengths!(unitcell::Unitcell, strengthNew)
    # go through all connections and modify their strength
    for c in unitcell.connections
        c[3] = strengthNew
    end
end
function setAllInteractionStrengths!(lattice::Lattice, strengthNew)
    # go through all connections and modify their strength
    for c in lattice.connections
        c[3] = strengthNew
    end
end
export setAllInteractionStrengths!



function mapInteractionStrengths!(unitcell::Unitcell, mapping; replace_in_strings=true, evaluate=false)
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
end
function mapInteractionStrengths!(lattice::Lattice, mapping; replace_in_strings=true, evaluate=false)
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
end
export mapInteractionStrengths!




#-----------------------------------------------------------------------------------------------------------------------------
#
#   Remove interactions with given strengths of the unitcell / lattice
#
#-----------------------------------------------------------------------------------------------------------------------------
function removeConnections!(lattice::Lattice, strengthToRemove=0.0)
    connections_new = Array[]
    for c in lattice.connections
        if c[3] != strengthToRemove
            push!(connections_new, c)
        end
    end
    lattice.connections = connections_new;
end
function removeConnections!(unitcell::Unitcell, strengthToRemove=0.0)
    connections_new = Array[]
    for c in unitcell.connections
        if c[3] != strengthToRemove
            push!(connections_new, c)
        end
    end
    unitcell.connections = connections_new;
end

export removeConnections!




# TODO remove sites
