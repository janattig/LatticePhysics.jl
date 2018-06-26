#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#
#   CALCULATIONS OF BAND STRUCTURE
#
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------







#-----------------------------------------------------------------------------------------------------------------------------
#
#   METHODS FOR CALCULATING THE BAND STRUCTURE OF A MATRIX IN K SPACE
#
#-----------------------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------------------
#
#   BAND STRUCTURE ALONG A PATH IN K SPACE
#
#   path has to be of following format
#   Array[
#       ["name1", [coordinates1]],
#       ["name2", [coordinates2]],
#       ...
#   ]
#   path does not close! to close, insert first point again
#
#   Parameters (necessary):
#   - lattice: The lattice object of which the band structure should be calculated
#   - path: The path along which the band structure should be calculated (in format given above)
#
#   Parameters (optional):
#   - reduceLattice: If the lattice should be reduced to a 1x1(x1) lattice of the original unitcell (i.e. for purposes of Luttinger Tisza, etc.)
#   - percentages: Either array of floats giving the individual percentages of the segments or the String "EQUAL" for equal length in the plot
#   - resolution: How many k-points to calculate in total
#   - enforce_hermitian: If the matrix should be enforced to be hermitian
#   - limits_energy: The y-axis (energy axis) limits of the plot
#   - plot_title: A title for the whole plot ("AUTO" for automatic title)
#   - plot_color: The color of bands
#   - figsize: the figure size as given to PyPlot
#
#-----------------------------------------------------------------------------------------------------------------------------
# BAND STRUCTURE ALONG PATH
function calculateBandStructureAlongPath(
        interaction_matrix::Array{Complex,2},
        path;
        percentages="EQUAL",
        resolution=1000,
        enforce_hermitian=false,
        limits_energy="AUTO",
        plot_title="",
        plot_color="b",
        figsize=(6,4),
        showPlot=true
            )

    # normalize percentages
    if percentages == "EQUAL"
        percentages = ones(size(path,1)-1)
    end
    percentages = percentages ./ sum(percentages)
    # build up the segment list
    segments = Array[]
    for i in 1:size(path,1)-1
        segment = [i, i+1, percentages[i]]
        push!(segments, segment)
    end
    # if LT is checked
    #if check_LT
    #    LT_k = Array[]
    #end
    # segment data, i.e. the bandstructure over the segments
    segments_data = Array[]
    resolution_actual = 0
    hlines = []
    # iterate over all segments
    for segment in segments
        # get the grid in between two points
        resolution_local = Int(floor(segment[3]*resolution))
        multipliers = linspace(0, 1, resolution_local)
        resolution_actual += resolution_local
        push!(hlines, resolution_actual+0)
        #println(segment)
        k1 = convert(Array{Float64,1}, path[Int(segment[1])][2:end])
        k2 = convert(Array{Float64,1}, path[Int(segment[2])][2:end])
        #println(k1)
        #println(k2)
        # insert bands
        bands = Array[]
        for b in 1:size(interaction_matrix,1)
            push!(bands, zeros(resolution_local))
        end
        # calculate all energies
        for i in 1:resolution_local
            # get the current k
            k = k2 .* multipliers[i] .+ k1 .* (1-multipliers[i])
            # if LT is checked, push current k
            #if check_LT
            #    push!(LT_k, k)
            #end
            # get the interaction matrix for this k
            matrix = copy(interaction_matrix)
            # diagonalize the matrix
            eigenvalues = eigvals(matrix)
            # save all the eigenvalues
            for b in 1:size(bands, 1)
                if imag(eigenvalues[b]) > 0
                    if imag(eigenvalues[b]) > 1e-15
                        println(imag(eigenvalues[b]))
                        println(matrix)
                        bands[b][i] = eigenvalues[b]
                    else
                        bands[b][i] = real(eigenvalues[b])
                    end
                else
                    bands[b][i] = eigenvalues[b]
                end
            end
        end
        # push the obtained back structure into the data array
        push!(segments_data, bands)
    end
    # generate the complete band structure
    bandstructure = Array[zeros(resolution_actual) for b in segments_data[1]]
    index = 1
    for i in 1:size(segments_data,1)
        segment = segments[i]
        data = segments_data[i]
        for b in 1:size(bandstructure,1)
            bandstructure[b][index:hlines[i]] = data[b]
        end
        index = hlines[i]+1
    end
    # if LT is checked, give the results
    #if check_LT
    #    LT_v = checkLuttingerTisza(lattice, LT_k, only_GS=false)
    #    println("$(100.0*sum(LT_v)/length(LT_v)) % of all eigenvalues are valid in LT")
    #end
    # plot the eigenvalues
    rc("font", family="serif")
    fig = figure(figsize=figsize)
    if plot_title == "AUTO"
        title("energy spectrum along path of interaction matrix")
    elseif plot_title == ""
        # do nothing title related
    else
        title(plot_title)
    end
    for l in hlines[1:end-1]
        axvline(l,color=[0.6, 0.6, 0.6], linestyle="--")
    end
    xlabel("momentum")
    ylabel("energy")
    for b in bandstructure
        plot(collect(1:resolution_actual), b, "-$(plot_color)")
    end
    ax = gca()
    axx = ax[:get_xaxis]()
    xtpos = []
    push!(xtpos, 0)
    for h in hlines
        push!(xtpos, h)
    end
    xtlabs = [p[1] for p in path]
    xticks(xtpos, xtlabs)
    #axx[:set_ticks]([])
    axx[:set_tick_params](which="both", direction="out")
    axx[:set_tick_params](which="top", color="none")
    axy = ax[:get_yaxis]()
    axy[:set_tick_params](which="both", direction="out")
    # check if specific boundaries are desired
    if !(limits_energy == "AUTO")
        ylim(limits_energy[1], limits_energy[2])
    end
    # tighten the layout
    tight_layout()
    # save the plot
    figurename = "interaction_matrix"
    figurename1 = "$(FOLDER_SPECTRA)bandstructure_path_$(figurename[1:end-4]).pdf"
    figurename2 = "$(FOLDER_SPECTRA)bandstructure_path_$(figurename[1:end-4]).png"
    buildFolderSpectra()
    savefig(figurename1)
    savefig(figurename2)
    if showPlot
        show()
        #print("Continue? ")
        #readline()
    end
    return fig
end
function calculateBandStructureAlongPath(
        lattice::Lattice,
        path;
        reduceLattice=true,
        percentages="EQUAL",
        resolution=1000,
        enforce_hermitian=false,
        limits_energy="AUTO",
        plot_title="",
        plot_color="b",
        figsize=(6,4),
        showPlot=true,
        majorana=false
            )

    # check if to reduce the lattice
    if reduceLattice && lattice.unitcell.filename != UNITCELL_DUMMY_FILENAME
        lattice = getLatticePeriodic(lattice.unitcell, ones(Int64, size(lattice.unitcell.lattice_vectors,1)), save=false)
    end
    # normalize percentages
    if percentages == "EQUAL"
        percentages = ones(size(path,1)-1)
    end
    percentages = percentages ./ sum(percentages)
    # build up the segment list
    segments = Array[]
    for i in 1:size(path,1)-1
        segment = [i, i+1, percentages[i]]
        push!(segments, segment)
    end
    # if LT is checked
    #if check_LT
    #    LT_k = Array[]
    #end
    # segment data, i.e. the bandstructure over the segments
    segments_data = Array[]
    resolution_actual = 0
    hlines = []
    # iterate over all segments
    for segment in segments
        # get the grid in between two points
        resolution_local = Int(floor(segment[3]*resolution))
        multipliers = linspace(0, 1, resolution_local)
        resolution_actual += resolution_local
        push!(hlines, resolution_actual+0)
        #println(segment)
        k1 = convert(Array{Float64,1}, path[Int(segment[1])][2:end])
        k2 = convert(Array{Float64,1}, path[Int(segment[2])][2:end])
        #println(k1)
        #println(k2)
        # insert bands
        bands = Array[]
        for b in 1:size(lattice.positions,1)
            push!(bands, zeros(resolution_local))
        end
        # calculate all energies
        for i in 1:resolution_local
            # get the current k
            k = k2 .* multipliers[i] .+ k1 .* (1-multipliers[i])
            # if LT is checked, push current k
            #if check_LT
            #    push!(LT_k, k)
            #end
            # get the interaction matrix for this k
            matrix = getInteractionMatrixKSpace(lattice, k, enforce_hermitian=enforce_hermitian, majorana=majorana)
            # diagonalize the matrix
            eigenvalues = eigvals(matrix)
            # save all the eigenvalues
            for b in 1:size(bands, 1)
                if imag(eigenvalues[b]) > 0
                    if imag(eigenvalues[b]) > 1e-15
                        println(imag(eigenvalues[b]))
                        println(matrix)
                        bands[b][i] = eigenvalues[b]
                    else
                        bands[b][i] = real(eigenvalues[b])
                    end
                else
                    bands[b][i] = eigenvalues[b]
                end
            end
        end
        # push the obtained back structure into the data array
        push!(segments_data, bands)
    end
    # generate the complete band structure
    bandstructure = Array[zeros(resolution_actual) for b in segments_data[1]]
    index = 1
    for i in 1:size(segments_data,1)
        segment = segments[i]
        data = segments_data[i]
        for b in 1:size(bandstructure,1)
            bandstructure[b][index:hlines[i]] = data[b]
        end
        index = hlines[i]+1
    end
    # if LT is checked, give the results
    #if check_LT
    #    LT_v = checkLuttingerTisza(lattice, LT_k, only_GS=false)
    #    println("$(100.0*sum(LT_v)/length(LT_v)) % of all eigenvalues are valid in LT")
    #end
    # plot the eigenvalues
    rc("font", family="serif")
    fig = figure(figsize=figsize)
    if plot_title == "AUTO"
        if majorana
            title("majorana energy spectrum along path of lattice \"$(lattice.filename)\"")
        else
            title("energy spectrum along path of lattice \"$(lattice.filename)\"")
        end
    elseif plot_title == ""
        # do nothing title related
    else
        title(plot_title)
    end
    for l in hlines[1:end-1]
        axvline(l,color=[0.6, 0.6, 0.6], linestyle="--")
    end
    xlabel("momentum")
    ylabel("energy")
    for b in bandstructure
        plot(collect(1:resolution_actual), b, "-$(plot_color)")
    end
    ax = gca()
    axx = ax[:get_xaxis]()
    xtpos = []
    push!(xtpos, 0)
    for h in hlines
        push!(xtpos, h)
    end
    xtlabs = [p[1] for p in path]
    xticks(xtpos, xtlabs)
    #axx[:set_ticks]([])
    axx[:set_tick_params](which="both", direction="out")
    axx[:set_tick_params](which="top", color="none")
    axy = ax[:get_yaxis]()
    axy[:set_tick_params](which="both", direction="out")
    # check if specific boundaries are desired
    if !(limits_energy == "AUTO")
        ylim(limits_energy[1], limits_energy[2])
    end
    # tighten the layout
    tight_layout()
    # save the plot
    figurename = split(lattice.filename, FOLDER_SPECTRA[end])[end]
    if majorana
        figurename1 = "$(FOLDER_SPECTRA)majorana_bandstructure_path_$(figurename[1:end-4]).pdf"
        figurename2 = "$(FOLDER_SPECTRA)majorana_bandstructure_path_$(figurename[1:end-4]).png"
    else
        figurename1 = "$(FOLDER_SPECTRA)bandstructure_path_$(figurename[1:end-4]).pdf"
        figurename2 = "$(FOLDER_SPECTRA)bandstructure_path_$(figurename[1:end-4]).png"
    end
    buildFolderSpectra()
    savefig(figurename1)
    savefig(figurename2)
    if showPlot
        show()
        #print("Continue? ")
        #readline()
    end
    return fig
end
function calculateBandStructureAlongPath(
        unitcell::Unitcell,
        path;
        percentages="EQUAL",
        resolution=1000,
        enforce_hermitian=false,
        limits_energy="AUTO",
        plot_title="",
        plot_color="b",
        figsize=(6,4),
        showPlot=true,
        majorana=false
            )

    # make a lattice from the unitcell
    lattice = getLatticePeriodic(unitcell, ones(Int64, size(unitcell.lattice_vectors,1)), save=false)
    lattice.filename = replace(identity(unitcell.filename), FOLDER_UNITCELLS, FOLDER_LATTICES)
    # push to lattice based method and return the result
    return calculateBandStructureAlongPath(
        lattice,
        path;
        reduceLattice=false,
        percentages=percentages,
        resolution=resolution,
        enforce_hermitian=enforce_hermitian,
        limits_energy=limits_energy,
        plot_title=plot_title,
        plot_color=plot_color,
        figsize=figsize,
        showPlot=showPlot,
        majorana=majorana
            )
end
export calculateBandStructureAlongPath

# BAND STRUCTURE OF STRIP (1D periodic)
function calculateBandStructureOfStrip(
        unitcell::Unitcell,
        periodic_direction::Int64,
        finite_direction_N::Int64;
        percentages="EQUAL",
        resolution=1000,
        enforce_hermitian=false,
        limits_energy="AUTO",
        plot_title="",
        plot_color="b",
        figsize=(6,4),
        showPlot=true,
        majorana=false
    )

    # obtain the lattice
    dimensions = ones(Int64, size(unitcell.lattice_vectors, 1)) .* finite_direction_N
    dimensions[periodic_direction] = -1
    lattice = getLattice(unitcell, dimensions)

    # obtain the path
    last_vector = lattice.lattice_vectors[1]
    last_vector = last_vector
    a = sqrt(sum(last_vector.*last_vector))
    gamma_zero  = last_vector .* 0.0
    gamma_minus = last_vector .* (-2*pi) / (a*a)
    gamma_plus  = last_vector .* ( 2*pi) / (a*a)
    K_minus     = last_vector .* (-1*pi) / (a*a)
    K_plus      = last_vector .* ( 1*pi) / (a*a)

    path = Array[
        ["Gamma (-1)";  gamma_minus],
        ["-pi";         K_minus],
        ["Gamma (0)";   gamma_zero],
        ["pi";          K_plus],
        ["Gamma (+1)";  gamma_plus]
    ]

    # calculate the band structure
    calculateBandStructureAlongPath(
        lattice,
        path;
        reduceLattice=false,
        percentages=percentages,
        resolution=resolution,
        enforce_hermitian=enforce_hermitian,
        limits_energy=limits_energy,
        plot_title=plot_title,
        plot_color=plot_color,
        figsize=figsize,
        showPlot=showPlot,
        majorana=majorana
    )
end
export calculateBandStructureOfStrip


# FULL 2D Bandstructure
function calculateBandStructure2D(
        unitcell::Unitcell,
        kx, ky;
        enforce_hermitian=false,
        limits_energy="AUTO",
        plot_title="",
        plot_color="b",
        figsize=(6,4),
        showPlot=true,
        majorana=false
            )


    # insert bands
    bands = Array[]
    kx_vals = zeros(length(kx), length(ky))
    ky_vals = zeros(length(kx), length(ky))
    for b in 1:size(unitcell.basis,1)
        push!(bands, zeros(length(kx), length(ky)))
    end
    # calculate all energies
    for i in 1:length(kx)
    for j in 1:length(ky)
        # get the current k
        k = [kx[i], ky[j]]
        kx_vals[i,j] = k[1]
        ky_vals[i,j] = k[2]
        # get the interaction matrix for this k
        matrix = getInteractionMatrixKSpace(unitcell, k, enforce_hermitian=enforce_hermitian, majorana=majorana)
        # diagonalize the matrix
        eigenvalues = eigvals(matrix)
        # save all the eigenvalues
        for b in 1:size(bands, 1)
            if imag(eigenvalues[b]) > 0
                if imag(eigenvalues[b]) > 1e-15
                    println(imag(eigenvalues[b]))
                    println(matrix)
                    bands[b][i,j] = eigenvalues[b]
                else
                    bands[b][i,j] = real(eigenvalues[b])
                end
            else
                bands[b][i,j] = eigenvalues[b]
            end
        end
    end
    end
    # generate the complete band structure
    bandstructure = bands
    # if LT is checked, give the results
    #if check_LT
    #    LT_v = checkLuttingerTisza(lattice, LT_k, only_GS=false)
    #    println("$(100.0*sum(LT_v)/length(LT_v)) % of all eigenvalues are valid in LT")
    #end
    # plot the eigenvalues
    rc("font", family="serif")
    fig = figure(figsize=figsize)
    ax = fig[:add_subplot](111, projection="3d")
    if plot_title == "AUTO"
        if majorana
            title("majorana energy spectrum along path of unitcell \"$(unitcell.filename)\"")
        else
            title("energy spectrum along path of unitcell \"$(unitcell.filename)\"")
        end
    elseif plot_title == ""
        # do nothing title related
    else
        title(plot_title)
    end
    xlabel("momentum kx")
    ylabel("momentum ky")
    zlabel("energy")
    for b in bandstructure
        ax[:plot_surface](kx_vals,ky_vals,b, rstride=1, cstride=1, cmap="coolwarm", linewidth=0)
    end
    axx = ax[:get_xaxis]()
    # check if specific boundaries are desired
    if !(limits_energy == "AUTO")
        ylim(limits_energy[1], limits_energy[2])
    end
    # tighten the layout
    tight_layout()
    # save the plot
    figurename = split(unitcell.filename, FOLDER_SPECTRA[end])[end]
    if majorana
        figurename1 = "$(FOLDER_SPECTRA)majorana_bandstructure_path_$(figurename[1:end-4]).pdf"
        figurename2 = "$(FOLDER_SPECTRA)majorana_bandstructure_path_$(figurename[1:end-4]).png"
    else
        figurename1 = "$(FOLDER_SPECTRA)bandstructure_path_$(figurename[1:end-4]).pdf"
        figurename2 = "$(FOLDER_SPECTRA)bandstructure_path_$(figurename[1:end-4]).png"
    end
    buildFolderSpectra()
    savefig(figurename1)
    savefig(figurename2)
    if showPlot
        show()
        #print("Continue? ")
        #readline()
    end
    return fig
end
export calculateBandStructure2D






# SOME DEFAULT PATHS
DEFAULT_PATH_FCC = Array[
    ["gamma"; [0,0,0]],
    ["X"; [2*pi, 0, 0]],
    ["W"; [2*pi, pi, 0]],
    ["L"; [pi, pi, pi]],
    ["gamma"; [0,0,0]],
    ["K"; [3*pi/2, 3*pi/2, 0]],
    ["X"; [2*pi, 0, 0]]
]
export DEFAULT_PATH_FCC

DEFAULT_PATH_TRIANGULAR = Array[
    ["gamma"; [0,0]],
    ["K"; [2*pi/sqrt(3.0), 2*pi/3]],
    ["M"; [2*pi/sqrt(3.0), 0]],
    ["gamma"; [0,0]]
]
export DEFAULT_PATH_TRIANGULAR

DEFAULT_PATH_SQUAREOCTAGON_2 = Array[
    ["gamma"; [0,0]],
    ["K"; [2,  0].*(pi / (1.0 + 1.0/sqrt(2.0)))],
    ["M"; [1, -1].*(pi / (1.0 + 1.0/sqrt(2.0)))],
    ["gamma"; [0,0]]
]
export DEFAULT_PATH_SQUAREOCTAGON_2


DEFAULT_PATH_SQUARE_LONG = Array[
    ["M";     [pi,   0.0]],
    ["Gamma"; [0.0,  0.0]],
    ["K'";    [pi,   -pi]],
    ["M";     [pi,   0.0]],
    ["M'";    [0.0,   pi]],
    ["K";     [pi,    pi]],
    ["Gamma"; [0.0,  0.0]]
]
export DEFAULT_PATH_SQUARE_LONG

DEFAULT_PATH_SQUARE_SHORT = Array[
    ["Gamma"; [0.0,  0.0]],
    ["M";     [pi,   0.0]],
    ["K";     [pi,    pi]],
    ["Gamma"; [0.0,  0.0]]
]
export DEFAULT_PATH_SQUARE_SHORT

DEFAULT_PATH_SQUARE = DEFAULT_PATH_SQUARE_SHORT
export DEFAULT_PATH_SQUARE
