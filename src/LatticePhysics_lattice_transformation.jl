#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#
#   LATTICE OPERATIONS (MODIFY LATTICE CONNECTIONS AND POSITIONS)
#
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------------------
#
#   Build a lattice-X / unitcell-X version of the lattice / unitcell
#   Procedure replaces all bonds by sites
#
#-----------------------------------------------------------------------------------------------------------------------------
function getTransformedLatticeBondToSite(lattice::Lattice; connection_strength="AUTO")
    # new positions and connections
    positions = copy(lattice.positions)
    connections = Array{Any,1}[]
    connectionsTreated = Array[]
    # generate a new filename
    filename_new = replace(lattice.filename, ".jld", "_X.jld")
    # neutral wrap
    neutral_wrap = (0,0,0)
    if size(lattice.lattice_vectors,1) == 2
   	    neutral_wrap = (0,0)
    elseif size(lattice.lattice_vectors,1) == 1
   	    neutral_wrap = (0)
    end
    # iterate over all connections
    for c in lattice.connections
        # only treat one of the two
        treated = false
        for ct in connectionsTreated
            if ct[1] == c[2] && ct[2] == c[1] && ct[3] == c[3]
                treated = true
                break
            end
        end
        if treated
            continue
        end
        # add a new position
        pnew = positions[Int(c[1])].+positions[Int(c[2])]
        for i in 1:length(c[4])
            pnew .+= c[4][i].*lattice.lattice_vectors[i]
        end
        pnew = pnew .* 0.5
        push!(positions, pnew)
        # add new connections
        if typeof(connection_strength) == String && connection_strength == "AUTO"
            con_strength = "con$(size(connectionsTreated,1))"
        elseif typeof(connection_strength) == String && connection_strength == "SAME"
            con_strength = c[3]
        elseif typeof(connection_strength) == String && connection_strength == "SQRT"
            if typeof(c[3]) == String
                con_strength = "sqrt($(c[3]))"
            else
                con_strength = sqrt(c[3])
            end
        else
            con_strength = connection_strength
        end
        push!(connections, [c[1], size(positions,1), con_strength, neutral_wrap])
        push!(connections, [size(positions,1), c[1], con_strength, neutral_wrap])
        push!(connections, [size(positions,1), c[2], con_strength, c[4]])
        if length(c[4]) == 3
            push!(connections, [c[2], size(positions,1), con_strength, (-c[4][1],-c[4][2],-c[4][3])])
        elseif length(c[4]) == 2
            push!(connections, [c[2], size(positions,1), con_strength, (-c[4][1],-c[4][2])])
        else
            push!(connections, [c[2], size(positions,1), con_strength, (-c[4])])
        end
        # add to treated
        push!(connectionsTreated, c)
    end
    # build a new lattice
    lattice_X = Lattice(lattice.lattice_vectors, positions, connections, filename_new)
    # return the new lattice
    return lattice_X
end
export getTransformedLatticeBondToSite

function getTransformedUnitcellBondToSite(unitcell::Unitcell; connection_strength="AUTO")
    # new positions and connections
    positions = copy(unitcell.basis)
    connections = Array{Any,1}[]
    connectionsTreated = Array[]
    # generate a new filename
    filename_new = replace(unitcell.filename, ".jld", "_X.jld")
    # neutral wrap
    neutral_wrap = (0,0,0)
    if size(unitcell.lattice_vectors,1) == 2
   	    neutral_wrap = (0,0)
    elseif size(unitcell.lattice_vectors,1) == 1
   	    neutral_wrap = (0)
    end
    # iterate over all connections
    for c in unitcell.connections
        # only treat one of the two
        treated = false
        for ct in connectionsTreated
            if ct[1] == c[2] && ct[2] == c[1] && ct[3] == c[3]
                treated = true
                break
            end
        end
        if treated
            continue
        end
        # add a new position
        pnew = positions[Int(c[1])].+positions[Int(c[2])]
        for i in 1:length(c[4])
            pnew .+= c[4][i].*unitcell.lattice_vectors[i]
        end
        pnew = pnew .* 0.5
        push!(positions, pnew)
        # add new connections
        if typeof(connection_strength) == String && connection_strength == "AUTO"
            con_strength = "con$(size(connectionsTreated,1))"
        elseif typeof(connection_strength) == String && connection_strength == "SAME"
            con_strength = c[3]
        elseif typeof(connection_strength) == String && connection_strength == "SQRT"
            if typeof(c[3]) == String
                con_strength = "sqrt($(c[3]))"
            else
                con_strength = sqrt(c[3])
            end
        else
            con_strength = connection_strength
        end
        push!(connections, [c[1], size(positions,1), con_strength, neutral_wrap])
        push!(connections, [size(positions,1), c[1], con_strength, neutral_wrap])
        push!(connections, [size(positions,1), c[2], con_strength, c[4]])
        if length(c[4]) == 3
            push!(connections, [c[2], size(positions,1), con_strength, (-c[4][1],-c[4][2],-c[4][3])])
        elseif length(c[4]) == 2
            push!(connections, [c[2], size(positions,1), con_strength, (-c[4][1],-c[4][2])])
        else
            push!(connections, [c[2], size(positions,1), con_strength, (-c[4])])
        end
        # add to treated
        push!(connectionsTreated, c)
    end
    # build a new unitcell
    uc_X = Unitcell(unitcell.lattice_vectors, positions, connections, filename_new)
    # return the new unitcell
    return uc_X
end
export getTransformedUnitcellBondToSite




#-----------------------------------------------------------------------------------------------------------------------------
#
#   Build a the squareroot version of the unitcell
#   Procedure replaces all fully connected plaquettes by sites and then sets the interaction strengths (TODO)
#
#-----------------------------------------------------------------------------------------------------------------------------
function getTransformedUnitcellFCPToSite(unitcell::Unitcell; printFCP=false, connection_strength="AUTO")
    # look for all fully connected plaquettes
    fcp_list = []
    connectionlist = getConnectionList(unitcell)
    neutral_wrap = (0,0,0)
    if length(connectionlist[1][1][4]) == 2
        neutral_wrap = (0,0)
    elseif length(connectionlist[1][1][4]) == 1
        neutral_wrap = (0)
    end
    # define a recursive looking method
    function getFCPs(fcp_current)
        # a site with index current_index is asked if it is in a fully connected plaquette with more sits than given by fcp
        # returns the fcps that the site is in (minimum: the fcp given)
        fcp_list_new = Array[]
        # iterate over all neighbors of the current site
        for nc in connectionlist[Int(fcp_current[end][1])]
            # current test site and test wrap
            test_site = Int(nc[2])
            test_wrap = collect(nc[4]) .+ fcp_current[end][2]
            # check, if "step back"
            if test_site == fcp_current[end-1][1] && test_wrap == fcp_current[end-1][2]
                continue
            end
            # no step back, check if this site is connected to all fcp sites that are already in the fcp
            completely_connected = true
            for site_in_fcp in fcp_current
                # check if new site is connected to this site
                connected = false
                for c in connectionlist[test_site]
                    # check if the destination and wrap agree
                    if Int(c[2]) == Int(site_in_fcp[1]) && collect(c[4]) .+ test_wrap == site_in_fcp[2]
                        # found a match
                        connected = true
                        break
                    end
                end
                # if not connected: complete connection broken
                if !connected
                    completely_connected = false
                    break
                end
            end
            # if it is connected to everything: make a new fcp out of it and go deper into the recursion
            if completely_connected
                fcp_new = []
                for fcp_current_element in fcp_current
                    push!(fcp_new, fcp_current_element)
                end
                push!(fcp_new, [test_site, test_wrap])
                # recursion
                fcp_list_recursion = getFCPs(fcp_new)
                # add those to the list
                for fcp in fcp_list_recursion
                    push!(fcp_list_new, fcp)
                end
            end
        end
        # if empty, insert the current fcp and return
        if length(fcp_list_new) == 0
            push!(fcp_list_new, fcp_current)
        end
        # return fcp list
        return fcp_list_new
    end
    # define a function to give connections for a certain fcp
    function getFCPConnections(fcp)
        # list
        fcp_connections = []
        # find all connections
        for site1 in fcp
        for site2 in fcp
            # add connection from 1 to 2 if found
            for c in connectionlist[Int(site1[1])]
                if c[1] == Int(site1[1]) && c[2] == Int(site2[1]) && collect(c[4]) == site2[2] .- site1[2]
                    push!(fcp_connections, c)
                end
            end
        end
        end
        # return the list of connections
        return fcp_connections
    end
    # define a function to compare to lists
    function sameConnections(list1, list2)
        # check if every element in list 1 is in list 2
        for element1 in list1
            # not found yet
            found = false
            # check all elements of list2
            for element2 in list2
                # check if equal
                if Int(element1[1]) == Int(element2[1]) && Int(element1[2]) == Int(element2[2]) && element1[4] == element2[4]
                    found = true
                    break
                end
            end
            # if not found, return false
            if !found
                return false
            end
        end
        # check if every element in list 2 is in list 1
        for element2 in list2
            # not found yet
            found = false
            # check all elements of list1
            for element1 in list1
                # check if equal
                if Int(element1[1]) == Int(element2[1]) && Int(element1[2]) == Int(element2[2]) && element1[4] == element2[4]
                    found = true
                    break
                end
            end
            # if not found, return false
            if !found
                return false
            end
        end
        # every test was successfull: 2 lists identical
        return true
    end
    # define a function to get the center of the fcp
    function getFCPCenter(fcp)
        # add all sites
        center = unitcell.basis[Int(fcp[1][1])]
        for i in 2:size(fcp,1)
            center = center .+ unitcell.basis[Int(fcp[i][1])]
        end
        # add all wraps
        for i in 1:size(fcp, 1)
            for l in 1:size(unitcell.lattice_vectors,1)
                center = center .+ (fcp[i][2][l] .* unitcell.lattice_vectors[l])
            end
        end
        # devide by number of sites
        center = center ./ size(fcp,1)
    end
    # iterate for all sites
    for s in 1:size(unitcell.basis,1)
        # looking for fully connected plaquettes of site s
        # check if neighbors are part of bigger plaquettes
        for neighbor_connection in connectionlist[s]
            # ask if the neighbor is part of an FCP
            fcp_new = [
                [s, collect(neutral_wrap)],
                [Int(neighbor_connection[2]), collect(neighbor_connection[4])]
            ]
            fcp_list_of_neighbor_and_site = getFCPs(fcp_new)
            # go through all fcps that have been identified and add new ones if necessary
            for fcp in fcp_list_of_neighbor_and_site
                # just add to global list, filtering later on
                push!(fcp_list, fcp)
            end
        end
    end
    # up to here: found all fcps
    fcp_list_raw = fcp_list
    fcp_list = []
    fcp_list_connections = []
    # now: filter the identical fcps out
    for fcp_test in fcp_list_raw
        # boolean if found in the ready to use list
        found = false
        # get the connections
        fcp_test_connections = getFCPConnections(fcp_test)
        # check all elements of the ready to use list
        for (i,fcp) in enumerate(fcp_list)
            # first test: dimensions identical?
            if size(fcp, 1) != size(fcp_test, 1)
                continue
            end
            # second test: get the connections associated with this certain fcp
            fcp_connections = fcp_list_connections[i]
            # test if the same lists, i.e. if both are in the other
            if sameConnections(fcp_connections, fcp_test_connections)
                found = true
                break
            end
        end
        # if not found, add to the list
        if !found
            push!(fcp_list, fcp_test)
            push!(fcp_list_connections, fcp_test_connections)
        end
    end
    # maybe print the FCP list
    if printFCP
        println("In total: $(size(fcp_list,1)) FCPs found")
        for i in 1:size(fcp_list,1)
            println("FCP $(i) with $(size(fcp_list[i], 1)) sites:")
            println(" - sites: $(fcp_list[i])")
            println(" - connections: $(fcp_list_connections[i])")
        end
    end
    # replace all fcps by inserting new sites
    positions_new = []
    for p in unitcell.basis
        push!(positions_new, p)
    end
    # push new sites
    for fcp in fcp_list
        # get the center of the fcp and add it as a new position
        push!(positions_new, getFCPCenter(fcp))
    end
    # define new connections
    connections_new = []
    # iterate over all sites / all fcps
    for (i, fcp) in enumerate(fcp_list)
        # new position index
        pos_index = size(unitcell.basis, 1) + i
        # define new connections for every site in the FCP (new position wrap is neutral wrap)
        for (s,site) in enumerate(fcp)
            # build the fields of the new 2 connections to insert
            con_field_1 = Int(site[1])
            con_field_2 = pos_index
            if typeof(connection_strength) == String && connection_strength == "AUTO"
                con_field_3 = "fcp$(i)con$(s)"
            else
                con_field_3 = connection_strength
            end
            con_field_4_array = collect(neutral_wrap) .- site[2]
            if length(con_field_4_array) == 3
                con_field_4_1 = (con_field_4_array[1], con_field_4_array[2], con_field_4_array[3])
                con_field_4_2 = (-con_field_4_array[1], -con_field_4_array[2], -con_field_4_array[3])
            elseif length(con_field_4_array) == 2
                con_field_4_1 = (con_field_4_array[1], con_field_4_array[2])
                con_field_4_2 = (-con_field_4_array[1], -con_field_4_array[2])
            elseif length(con_field_4_array) == 1
                con_field_4_1 = (con_field_4_array[1])
                con_field_4_2 = (-con_field_4_array[1])
            else
                println("wrong length for field 4: $(con_field_4_array)")
            end
            # build new connections
            con_1 = [con_field_1; con_field_2; con_field_3; con_field_4_1]
            con_2 = [con_field_2; con_field_1; con_field_3; con_field_4_2]
            push!(connections_new, con_1)
            push!(connections_new, con_2)
        end
    end
    # build new unitcell
    unitcell_fcp_replaced = Unitcell(
        unitcell.lattice_vectors,
        positions_new,
        connections_new,
        replace(unitcell.filename, ".jld", "_fcp_to_site.jld"))
    # return the new unitcell
    return unitcell_fcp_replaced
end
export getTransformedUnitcellFCPToSite

#function getTransformedLatticeFCPToSite(lattice::Lattice)
#
#end
#export getTransformedLatticeFCPToSite



#-----------------------------------------------------------------------------------------------------------------------------
#
#   Square the lattice / unitcell
#   This maps all connections to NNN connections by multiplication
#
#-----------------------------------------------------------------------------------------------------------------------------
function getSquaredLattice(lattice::Lattice)
    # get the connectivity matrix of the lattice
    connectivity = getConnectionList(lattice)
    # define a list of new connections
    connections_new = Array[]
    # iterate over all sites and check if they host a NNN connection
    for i in 1:size(lattice.positions,1)
        # check the neighbors of site i
        for (i1,c1) in enumerate(connectivity[i])
        for (i2,c2) in enumerate(connectivity[i])
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
                if typeof(c1[3]) == String || typeof(c2[3]) == String
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); "$(c1[3])*$(c2[3])"; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); "$(c1[3])*$(c2[3])"; off_2]
                else
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); (c1[3])*(c2[3]); off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); (c1[3])*(c2[3]); off_2]
                end
                # push them to the add list
                push!(connections_new, connection_new_1)
                push!(connections_new, connection_new_2)
            end
        end
        end
    end
    # create new Lattice with only new connections (and leave out old ones)
    return getLatticeWithOptimizedConnections(
        Lattice(
            lattice.unitcell,
            lattice.unitcellRepetitions,
            lattice.lattice_vectors,
            lattice.positions,
            lattice.positions_indices,
            connections_new,
            replace(lattice.filename, ".jld", "_squared.jld"))
    )
end
export getSquaredLattice

function getSquaredUnitcell(unitcell::Unitcell)
    # get the connectivity matrix of the unitcell
    connectivity = getConnectionList(unitcell)
    # define a list of new connections
    connections_new = Array[]
    # iterate over all sites and check if they host a NNN connection
    for i in 1:size(unitcell.basis,1)
        # check the neighbors of site i
        for (i1,c1) in enumerate(connectivity[i])
        for (i2,c2) in enumerate(connectivity[i])
            # propose a connection c1+c2 if the connections c1 and c2 are different
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
                if typeof(c1[3]) == String || typeof(c2[3]) == String
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); "$(c1[3])*$(c2[3])"; off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); "$(c1[3])*$(c2[3])"; off_2]
                else
                    connection_new_1 = [Int(c1[2]); Int(c2[2]); (c1[3])*(c2[3]); off_1]
                    connection_new_2 = [Int(c2[2]); Int(c1[2]); (c1[3])*(c2[3]); off_2]
                end
                # push them to the add list
                push!(connections_new, connection_new_1)
                push!(connections_new, connection_new_2)
            end
        end
        end
    end
    # create new Lattice with only new connections (and leave out old ones)
    return getUnitcellWithOptimizedConnections(
        Unitcell(
            unitcell.lattice_vectors,
            unitcell.basis,
            connections_new,
            replace(unitcell.filename, ".jld", "_squared.jld"))
    )
end
export getSquaredUnitcell






#-----------------------------------------------------------------------------------------------------------------------------
#
#   Get a list of all independent sublattices / sublattices of a lattice / unitcell (i.e. not connected parts of a lattice)
#   Not connected parts can occur when "squaring" a bipartite lattice
#
#-----------------------------------------------------------------------------------------------------------------------------
function getIndependentSublattices(lattice::Lattice)
    # save relevant lattice data
    connections = lattice.connections
    positions = lattice.positions
    positions_indices = lattice.positions_indices
    # generate a list for the sublattices
    sublattices = []
    # define an array of labels
    labels = zeros(Int64, size(lattice.positions,1))
    currentlabel = 1
    # get the connectivity list of the lattice
    connectivity = getConnectionList(lattice)
    # iterate over all positions
    for i in 1:size(lattice.positions, 1)
        # get the label of site i
        if labels[i] == 0
            labels[i] = currentlabel
            currentlabel = currentlabel + 1
        end
        # iterate over all connections
        for c in connectivity[i]
            # get the connected site
            j = Int(c[2])
            # check the label of j
            if labels[j] == 0
                # set as connected to i
                labels[j] = labels[i]
            elseif labels[j] == labels[i]
                # just ignore
            else
                # merging of two clusters here!
                newlabel = min(labels[i], labels[j])
                # scan through the whole lattice
                for n in 1:length(labels)
                    if labels[n] == labels[i] || labels[n] == labels[j]
                        labels[n] = newlabel
                    end
                end
                # get a new currentlabel
                currentlabel = Int(maximum(labels)) + 1
            end
        end
    end
    # get a list of labels
    labellist = []
    for l in labels
        if !(l in labellist)
            push!(labellist, l)
        end
    end
    # print how many sublattices found
    println("$(length(labellist)) independent sublattice(s) found")
    # check how many sublattices
    if length(labellist) == 1
        # just push the original lattice into the sublattices list
        sublattice = Lattice(
            lattice.unitcell,
            lattice.unitcellRepetitions,
            lattice.lattice_vectors,
            lattice.positions,
            lattice.positions_indices,
            lattice.connections,
            lattice.filename_output
        )
        push!(sublattices, sublattice)
        # return
        return sublattices
    end
    # iterate over labels
    for l in labellist
        # find the new positions
        positions_new = Array[]
        positions_indices_new = []
        connections_new = Array[]
        mapping_indices_old = []
        # iterate over all connections to identify relevant positions
        for c in connections
            # check if belonging to the sublattice
            if labels[Int(c[1])] != l
                # not of this sublattice
                continue
            end
            # get the starting and finishing position
            index_from = Int(c[1])
            index_to = Int(c[2])
            # check where they map to
            if index_from in mapping_indices_old
                index_from_new = findfirst(mapping_indices_old, index_from)
            else
                # add the position to the list of positions
                push!(positions_new, positions[index_from])
                push!(positions_indices_new, positions_indices[index_from])
                # add the index mapping to the list
                push!(mapping_indices_old, index_from)
                # get the new index
                index_from_new = size(positions_new,1)
            end
            if index_to in mapping_indices_old
                index_to_new = findfirst(mapping_indices_old, index_to)
            else
                # add the position to the list of positions
                push!(positions_new, positions[index_to])
                push!(positions_indices_new, positions_indices[index_to])
                # add the index mapping to the list
                push!(mapping_indices_old, index_to)
                # get the new index
                index_to_new = size(positions_new,1)
            end
            # push the new connection
            # format [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            connection_new = [index_from_new; index_to_new; c[3]; c[4]]
			# register as connection
			push!(connections_new, connection_new)
        end
        # generate a new sublattice
        sublattice = Lattice(
            lattice.unitcell,
            lattice.unitcellRepetitions,
            lattice.lattice_vectors,
            positions_new,
            positions_indices_new,
            connections_new,
            "$(lattice.filename[1:end-4])_sublattice_$(l).jld"
        )
        # push it the list
        push!(sublattices, sublattice)
    end
    # return the list
    return sublattices
end
export getIndependentSublattices

function getIndependentSubunitcells(unitcell::Unitcell)
    # save relevant lattice data
    connections = unitcell.connections
    positions = unitcell.basis
    # generate a list for the sublattices
    subunitcells = []
    # define an array of labels
    labels = zeros(Int64, size(positions,1))
    currentlabel = 1
    # get the connectivity list of the unitcell
    connectivity = getConnectionList(unitcell)
    # iterate over all positions
    for i in 1:size(positions, 1)
        # get the label of site i
        if labels[i] == 0
            labels[i] = currentlabel
            currentlabel = currentlabel + 1
        end
        # iterate over all connections
        for c in connectivity[i]
            # get the connected site
            j = Int(c[2])
            # check the label of j
            if labels[j] == 0
                # set as connected to i
                labels[j] = labels[i]
            elseif labels[j] == labels[i]
                # just ignore
            else
                # merging of two clusters here!
                newlabel = min(labels[i], labels[j])
                # scan through the whole lattice
                for n in 1:length(labels)
                    if labels[n] == labels[i] || labels[n] == labels[j]
                        labels[n] = newlabel
                    end
                end
                # get a new currentlabel
                currentlabel = Int(maximum(labels)) + 1
            end
        end
    end
    # get a list of labels
    labellist = []
    for l in labels
        if !(l in labellist)
            push!(labellist, l)
        end
    end
    # print how many sublattices found
    println("$(length(labellist)) independent sublattice(s) found")
    # check how many sublattices
    if length(labellist) == 1
        # just push the original lattice into the sublattices list
        subunitcell = Unitcell(
            unitcell.lattice_vectors,
            unitcell.basis,
            unitcell.connections,
            unitcell.filename
        )
        push!(subunitcells, subunitcell)
        # return
        return subunitcells
    end
    # iterate over labels
    for l in labellist
        # find the new positions
        positions_new = Array[]
        connections_new = Array[]
        mapping_indices_old = []
        # iterate over all connections to identify relevant positions
        for c in connections
            # check if belonging to the sublattice
            if labels[Int(c[1])] != l
                # not of this sublattice
                continue
            end
            # get the starting and finishing position
            index_from = Int(c[1])
            index_to = Int(c[2])
            # check where they map to
            if index_from in mapping_indices_old
                index_from_new = findfirst(mapping_indices_old, index_from)
            else
                # add the position to the list of positions
                push!(positions_new, positions[index_from])
                # add the index mapping to the list
                push!(mapping_indices_old, index_from)
                # get the new index
                index_from_new = size(positions_new,1)
            end
            if index_to in mapping_indices_old
                index_to_new = findfirst(mapping_indices_old, index_to)
            else
                # add the position to the list of positions
                push!(positions_new, positions[index_to])
                # add the index mapping to the list
                push!(mapping_indices_old, index_to)
                # get the new index
                index_to_new = size(positions_new,1)
            end
            # push the new connection
            # format [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            connection_new = [index_from_new; index_to_new; c[3]; c[4]]
			# register as connection
			push!(connections_new, connection_new)
        end
        # generate a new sublattice
        subunitcell = Unitcell(
            unitcell.lattice_vectors,
            positions_new,
            connections_new,
            replace(unitcell.filename, ".jld", "_subunitcell_$(l).jld")
        )
        # push it the list
        push!(subunitcells, subunitcell)
    end
    # return the list
    return subunitcells
end
export getIndependentSubunitcells






#-----------------------------------------------------------------------------------------------------------------------------
#
#   Optimize the lattice / unitcell by collapsing redundant connections into fewer connections and return a new object
#
#-----------------------------------------------------------------------------------------------------------------------------
function getLatticeWithOptimizedConnections(lattice::Lattice)
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
    # return a new lattice
    return Lattice(
        lattice.unitcell,
        lattice.unitcellRepetitions,
        lattice.lattice_vectors,
        lattice.positions,
        lattice.positions_indices,
        connections_new_2,
        lattice.filename)
end
export getLatticeWithOptimizedConnections

function getUnitcellWithOptimizedConnections(unitcell::Unitcell)
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
    # return a new unitcell
    return Unitcell(
        unitcell.lattice_vectors,
        unitcell.basis,
        connections_new_2,
        unitcell.filename)
end
export getUnitcellWithOptimizedConnections
