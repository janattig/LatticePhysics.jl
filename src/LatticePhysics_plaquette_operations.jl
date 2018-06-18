################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT PLAQUETTE RELATED FUNCTIONS
#
#   STRUCTURE OF THE FILE
#
#   1) FINDING PLAQUETTES IN LATTICE OBJECTS
#   - of a site
#   - of complete lattice
#
#   2) PRINTING PLAQUETTE STATISTICS
#   - of a site
#   - of complete lattice
#
#
#   BUG: Plaquettes dont implement periodic boundaries so far
#
################################################################################





#-----------------------------------------------------------------------------------------------------------------------------
#
#   Find Plaquettes of some site in the lattice with desired length len
#
#-----------------------------------------------------------------------------------------------------------------------------


"""
    getPlaquettesOfSite(lattice::Lattice, site::Int64, len::Int64)

Function to find the plaquettes that loop through a given site with index `site`
of a passed `Lattice` object and have length `len`.

NOTE: Plaquettes cannot wrap around periodic boundaries correctly so far (in the sense that
they cannot visit sites with the same label in different periodic copies of the lattice)

The plaquettes are returned as a list of Arrays, where each element in the array is a site-index.



# Examples

```julia-repl
julia> plaquettes = getPlaquettesOfSite(lattice, 1, 6)
Array{Array{Int64, 1}}[
    [1, 5, 7, 22, 10, 17, 1],
    ...
]
```
"""
function getPlaquettesOfSite(lattice::Lattice, site::Int64, len::Int64)
    # get the connectivity of the lattice
    cl = getConnectionList(lattice)
    # history includes all sites that are already visited, including the current site
    function givePlaquettePaths(history::Array{Int64,1}, steps::Int64, patharray::Array{Array{Int64,1}, 1})
        # check if the maximum depth of the search has been reached
        if length(history) == steps+1
            # check if a closed loop exists for the maximum depth
            if history[end] == history[1]
                push!(patharray, history)
            end
            return
        elseif history[end] in history[1:end-1]
            # building a closed loop of not proper length
            return
        end
        # current site
        currentSite = history[end]
        # get all options
        #options = []
        #for c in cl[currentSite]
        #    push!(options, Int(c[2]))
        #end
        options = Int64[Int(c[2]) for c in cl[currentSite]]
        # distinguish between origin and further apart sites
        if length(history) == 1
            # go for all options
            for o in options
                # get the new history
                history_new = copy(history)
                push!(history_new, o)
                # get all subpaths
                givePlaquettePaths(history_new, steps, patharray)
            end
        else
            # go for all options that dont go back
            for o in options
                # check if going back
                if o == history[end-1]
                    continue
                end
                # get the new history
                history_new = copy(history)
                push!(history_new, o)
                # get all subpaths
                givePlaquettePaths(history_new, steps, patharray)
            end
        end
        # return
        return
    end
    # call the function to generate array of plaquettes
    paths = Array{Int64,1}[]
    givePlaquettePaths([site], len, paths)
    # check all paths if they are plaquettes
    plaquettes = Array{Int64,1}[]
    for p in paths
        if p[1] == p[end]
            # check if the inverse path is already in the plaquettes list
            if !(p[end:-1:1] in plaquettes)
                push!(plaquettes, p)
            end
        end
    end
    # return the array
    return plaquettes
end
export getPlaquettesOfSite



"""
    getPlaquettesOfLattice(lattice::Lattice, len::Int64)

Function to find all plaquettes which have length `len` of a given `Lattice` object.

NOTE: Plaquettes cannot wrap around periodic boundaries correctly so far (in the sense that
they cannot visit sites with the same label in different periodic copies of the lattice)

The plaquettes are returned as a list of Arrays, where each element in the array is a site-index.



# Examples

```julia-repl
julia> plaquettes = getPlaquettesOfLattice(lattice, 6)
Array{Array{Int64, 1}}[
    [1, 5, 7, 22, 10, 17, 1],
    ...
]
```
"""
function getPlaquettesOfLattice(lattice::Lattice, len::Int64)
    # get all plaquettes of all sites
    plaquettes = Array{Int64,1}[]
    plaquettes_sorted = Array{Int64,1}[]
    for i in 1:size(lattice.positions,1)
        # obtain plaquettes
        plaquettes_site = getPlaquettesOfSite(lattice, i, len)
        # push plaquettes into array
        for p in plaquettes_site
            # check if present
            p_sorted = copy(p[1:end-1])
            sort!(p_sorted)
            # search if present
            found = false
            for ps in plaquettes_sorted
                # check if lists equal
                if findfirst([ps[i] == p_sorted[i] for i in 1:len], false) == 0
                    found = true
                    break
                end
            end
            if !(found)
                push!(plaquettes, p)
                push!(plaquettes_sorted, p_sorted)
            end
        end
    end
    # return the list
    return plaquettes
end
export getPlaquettesOfLattice









#-----------------------------------------------------------------------------------------------------------------------------
#
#   Obtain and print information on plaquettes
#
#-----------------------------------------------------------------------------------------------------------------------------

"""
    printPlaquetteStatisticsOfSite(lattice::Lattice, site::Int64 [ ; detailed::Bool=false, max_length::Int64=12])

Function to find and print information about all plaquettes going through a given site of a passed `Lattice` object.
If `detailed=true`, the individual plaquettes are printed as well.
The option `max_length` specifies what maximal length of plaquettes should be detected.

NOTE: Plaquettes cannot wrap around periodic boundaries correctly so far (in the sense that
they cannot visit sites with the same label in different periodic copies of the lattice)


# Examples

```julia-repl
julia> plaquettes = printPlaquetteStatisticsOfSite(lattice, 1)
...

julia> plaquettes = printPlaquetteStatisticsOfSite(lattice, 1, detailed=true, max_length=10)
...
```
"""
function printPlaquetteStatisticsOfSite(lattice::Lattice, site::Int64; detailed::Bool=false, max_length::Int64=12)
    # print header information
    println("printing plaquette information for site $(site) of lattice:")
    # number of total plaquettes
    plaquettes_total = 0
    # iterate over all lengths
    for l in 3:max_length
        # get plaquettes
        plaquettes = getPlaquettesOfSite(lattice, site, l)
        # print statistics
        if detailed
            # definately print the number of plaquettes
            println("- length $(l) bonds: $(size(plaquettes,1)) plaquettes")
            # print more information
            for p in plaquettes
                println("    - $(p)")
            end
        else
            # maybe print the number of plaquettes
            if size(plaquettes,1) > 0
                println("- length $(l) bonds: $(size(plaquettes,1)) plaquettes")
            end
        end
        # add up how many plaquettes were found
        plaquettes_total += size(plaquettes,1)
    end
    # print total number
    println("In total: $(plaquettes_total) plaquettes found")
end
export printPlaquetteStatisticsOfSite


"""
    printPlaquetteStatisticsOfLattice(lattice::Lattice [ ; detailed::Bool=false, max_length::Int64=12])

Function to find and print information about all plaquettes of a passed `Lattice` object.
If `detailed=true`, the individual plaquettes are printed as well.
The option `max_length` specifies what maximal length of plaquettes should be detected.

NOTE: Plaquettes cannot wrap around periodic boundaries correctly so far (in the sense that
they cannot visit sites with the same label in different periodic copies of the lattice)


# Examples

```julia-repl
julia> plaquettes = printPlaquetteStatisticsOfLattice(lattice)
...

julia> plaquettes = printPlaquetteStatisticsOfLattice(lattice, detailed=true, max_length=10)
...
```
"""
function printPlaquetteStatisticsOfLattice(lattice::Lattice; detailed::Bool=false, max_length::Int64=12)
    # print header information
    println("printing plaquette information for lattice:")
    # number of total plaquettes
    plaquettes_total = 0
    # iterate over all lengths
    for l in 3:max_length
        # get plaquettes
        plaquettes = getPlaquettesOfLattice(lattice, l)
        # print statistics
        if detailed
            # definately print the number of plaquettes
            println("- length $(l) bonds: $(size(plaquettes,1)) plaquettes")
            # print more information
            for p in plaquettes
                println("    - $(p)")
            end
        else
            # maybe print the number of plaquettes
            if size(plaquettes,1) > 0
                println("- length $(l) bonds: $(size(plaquettes,1)) plaquettes")
            end
        end
        # add up how many plaquettes were found
        plaquettes_total += size(plaquettes,1)
    end
    # print total number
    println("In total: $(plaquettes_total) plaquettes found")
end
export printPlaquetteStatisticsOfLattice
