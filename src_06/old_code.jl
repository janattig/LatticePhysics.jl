
#-----------------------------------------------------------------------------------------------------------------------------
#
#   Interaction matrix INFORMATION
#   (returns two matrices that indicate the number of periodic and non-periodic connections between respective sites)
#
#   Parameters are
#   - lattice: The complete lattice of which the interaction matrix should be constructed
#   or
#   - unitcell: The unitcell of which the interaction matrix should be constructed
#
#-----------------------------------------------------------------------------------------------------------------------------
function getInteractionMatrixInformation(lattice::Lattice)
    # matrices indicating the number of connections
    con_periodic =  zeros(Int64, size(lattice.positions,1),size(lattice.positions,1))
    con_inside   =  zeros(Int64, size(lattice.positions,1),size(lattice.positions,1))
    # go through all connections of the lattice
    for c in lattice.connections
        # get the indices
        index_from  = Int(c[1])
        index_to    = Int(c[2])
        strength    = c[3]
        wrap        = c[4]
        # check if inside or periodic and add to the counter
        if sum([abs(el) for el in wrap]) == 0
            con_inside[index_from, index_to] = con_inside[index_from, index_to] + 1
        else
            con_periodic[index_from, index_to] = con_periodic[index_from, index_to] + 1
        end
    end
    # return the matrices
    return con_inside, con_periodic
end
function getInteractionMatrixInformation(unitcell::Unitcell)
    # matrices indicating the number of connections
    con_periodic =  zeros(Int64, size(unitcell.basis,1),size(unitcell.basis,1))
    con_inside   =  zeros(Int64, size(unitcell.basis,1),size(unitcell.basis,1))
    # go through all connections of the lattice
    for c in unitcell.connections
        # get the indices
        index_from  = Int(c[1])
        index_to    = Int(c[2])
        strength    = c[3]
        wrap        = c[4]
        # check if inside or periodic and add to the counter
        if sum([abs(el) for el in wrap]) == 0
            con_inside[index_from, index_to] = con_inside[index_from, index_to] + 1
        else
            con_periodic[index_from, index_to] = con_periodic[index_from, index_to] + 1
        end
    end
    # return the matrices
    return con_inside, con_periodic
end


# print the interaction matrix information in compressed LaTeX format
function printInteractionMatrixInformation(connections_inside, connections_periodic)
    # get the dimension of the matrices
    len = size(connections_periodic,1)
    # starting lines
    println("\\begin{equation*}")
    println("   \\mathbf{H} = \\bordermatrix{")
    print("       ")
    for i in 1:len
        print("& $(i) ")
    end
    println("\\cr")
    # iterate over all lines of the matrix
    for j in 1:len
        # start line by printing the current index
        print("       $(j) ")
        # go through all respective numbers
        for i in 1:len
            # check if periodic or internal connections
            if connections_inside[i,j] == 0 && connections_periodic[i,j] == 0
                print("& 0 ")
            elseif connections_inside[i,j] > 0 && connections_periodic[i,j] == 0
                #print("& {\\color{blue} 0},{\\color{red} $(connections_inside[i,j])} ")
                print("& {\\color{blue} $(connections_inside[i,j])} ")
            elseif connections_inside[i,j] == 0 && connections_periodic[i,j] > 0
                #print("& {\\color{blue} $(connections_periodic[i,j])},{\\color{red} 0} ")
                print("& {\\color{red} $(connections_periodic[i,j])} ")
            else
                print("& {\\color{red} $(connections_periodic[i,j])},{\\color{blue} $(connections_inside[i,j])} ")
            end
        end
        # end the line by printing a linebreak
        println("\\cr")
    end
    # end lines
    println("   }")
    println("\\end{equation*}")
end
function printInteractionMatrixInformation(lattice::Lattice)
    println("Interaction matrix for lattice \"$(lattice.filename)\"")
    connections_inside, connections_periodic = getInteractionMatrixInformation(lattice)
    printInteractionMatrixInformation(connections_inside, connections_periodic)
end
function printInteractionMatrixInformation(unitcell::Unitcell)
    println("Interaction matrix for unitcell \"$(unitcell.filename)\"")
    connections_inside, connections_periodic = getInteractionMatrixInformation(unitcell)
    printInteractionMatrixInformation(connections_inside, connections_periodic)
end
export printInteractionMatrixInformation




# function to show the lattice
function showLattice(
		lattice::Lattice;
		conversion = 160,
		site_radius=25,
        site_labels="OFF",
		bond_thickness::Int64=8,
		visualize_periodic=false,
		colorcode_sites = Dict(0 => [255,255,255], 1 => [255,255,255]),
		colorcode_bonds = Dict("0" => [0,0,0], "1.0" => [0,0,0]),
        colorcode_bonds_automation::String = "OFF",
        background_color = [200,200,200],
        forcePyPlot=false
    )
    # if forcing PyPlot, do this and exit
    if forcePyPlot
        showLatticePyPlot(
            lattice,
            conversion = conversion,
            site_radius=site_radius,
            site_labels=site_labels,
            bond_thickness=bond_thickness,
            visualize_periodic=visualize_periodic,
            colorcode_sites = colorcode_sites,
            colorcode_bonds = colorcode_bonds,
            colorcode_bonds_automation=colorcode_bonds_automation,
            background_color = background_color
        )
        return
    end
    # CHECK IF MAYAVI IS PRESENT
    try
        @pyimport mayavi.mlab as mlab
        MAYAVI_AVAILABLE = true
        try
        showLatticeMayavi(
            lattice,
            conversion = conversion,
            site_radius=site_radius,
            site_labels=site_labels,
            bond_thickness=bond_thickness,
            visualize_periodic=visualize_periodic,
            colorcode_sites = colorcode_sites,
            colorcode_bonds = colorcode_bonds,
            colorcode_bonds_automation=colorcode_bonds_automation,
            background_color = background_color
        )
        catch x
        end
        return
    catch error
        if isa(error, PyCall.PyError)
            #println("PyError occured when importing:")
            #println(error)
            MAYAVI_AVAILABLE = false
            #println("Using PyPlot now")
            try
                showLatticePyPlot(
                    lattice,
                    conversion = conversion,
                    site_radius=site_radius,
                    site_labels=site_labels,
                    bond_thickness=bond_thickness,
                    visualize_periodic=visualize_periodic,
                    colorcode_sites = colorcode_sites,
                    colorcode_bonds = colorcode_bonds,
                    colorcode_bonds_automation=colorcode_bonds_automation,
                    background_color = background_color
                )
            catch error2
                println("Error occured when plotting:")
                println(error2)
            end
        else
            println("Error occured when loading:")
            println(error)
        end
    end

end
export showLattice


# show Lattice function to plot in pyplot
function showLatticePyPlot(
		lattice::Lattice;
		conversion = 160,
		site_radius=25,
        site_labels="OFF",
		bond_thickness::Int64=8,
		visualize_periodic=false,
		colorcode_sites = Dict(0 => [255,255,255], 1 => [255,255,255]),
		colorcode_bonds = Dict("0" => [0,0,0], "1.0" => [0,0,0]),
        colorcode_bonds_automation::String = "OFF",
        background_color = [200,200,200]
    )

    # maybe overwrite dictonary
    if colorcode_bonds_automation == "GREY"
        # construct the list of interaction strengths
        cs_list = getConnectionStrengthList(lattice)
        # get the color code list
        cc_list = getGreySequence(size(cs_list, 1))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[string(cs_list[i])] = cc_list[i]
        end
    elseif colorcode_bonds_automation == "COLOR"
        # construct the list of interaction strengths
        cs_list = getConnectionStrengthList(lattice)
        # get the color code list
        cc_list = getColorSequence(size(cs_list, 1))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[string(cs_list[i])] = cc_list[i]
        end
    end

	# repair color dictonary
	colorcode_bonds["0"] = get(colorcode_bonds, "0", [0,0,0])
	colorcode_sites[0] = get(colorcode_sites, 0, [255,255,255])

    # PLOTTING WITH PYPLOT
    # create a figure
    fig = figure(figsize=(10,10), facecolor=background_color ./ 255.0)
    # add a 3d projection subplot
    ax = fig[:add_subplot](111, projection="3d", axisbg=background_color ./ 255.0)
    # define a function to plot a sphere for every site
    function plotSite(site,radius,color; detail=10)
        x = site[1]
        y = site[2]
        if length(site) == 3
            z = site[3]
        else
            z = 0.0
        end
        # Make data
        n = detail
        u = linspace(0,2*pi,n);
        v = linspace(0,pi,n);
        px = cos(u) * sin(v)';
        py = sin(u) * sin(v)';
        pz = ones(n) * cos(v)';
        px = (px.*radius) .+ x
        py = (py.*radius) .+ y
        pz = (pz.*radius) .+ z
        # Plot the surface
        surf(px,py,pz, rstride=1, cstride=1, linewidth=0, antialiased=true, color=color)
    end
    # define a function to plot a bond as a tube
    function plotBond(from, to, radius, color; detail_ring=8, detail_length=3)
        # check to be 3d
        if length(from) == 2
            from = [from[1], from[2], 0.0]
        end
        if length(to) == 2
            to = [to[1], to[2], 0.0]
        end
        # https://de.mathworks.com/matlabcentral/answers/151235-tube-plot-with-x-y-coordinates-and-radius?requestedDomain=www.mathworks.com
        # plot bond
        vec_u = to .- from;
        t = nullspace(vec_u')';
        vec_v = t[1,:]
        vec_w = t[2,:]
        # normalize
        vec_v = radius.* vec_v ./ norm(vec_v)
        vec_w = radius.* vec_w ./ norm(vec_w)
        #println("v=$(vec_v)    w=$(vec_w)")
        #[S,T] = linspace(0,1,detail)*linspace(0,2*pi,detail)';
        S = linspace(0,1,detail_length)
        T = linspace(0,2*pi,detail_ring)
        px = zeros(detail_length, detail_ring)
        py = zeros(detail_length, detail_ring)
        pz = zeros(detail_length, detail_ring)
        for i in 1:detail_length
        for j in 1:detail_ring
            #P = repmat(from',m*n,1) .+ S.*vec_u .+ cos(T).*vec_v .+ sin(T).*vec_w;
            px[i,j] = from[1] + S[i]*vec_u[1] + cos(T[j])*vec_v[1] + sin(T[j])*vec_w[1];
            py[i,j] = from[2] + S[i]*vec_u[2] + cos(T[j])*vec_v[2] + sin(T[j])*vec_w[2];
            pz[i,j] = from[3] + S[i]*vec_u[3] + cos(T[j])*vec_v[3] + sin(T[j])*vec_w[3];
        end
        end
        surf(px,py,pz, rstride=1, cstride=1, linewidth=0, antialiased=true, color=color)
    end

    # plot all bonds
    for c in lattice.connections
        # skip if wrong directioin
        if c[1] < c[2] || (!visualize_periodic && sum([abs(el) for el in c[4]]) != 0)
            continue
        end
        # get data
        from = lattice.positions[Int(c[1])] .* conversion
        to = lattice.positions[Int(c[2])] .* conversion
        color = get(colorcode_bonds, string(c[3]), colorcode_bonds["0"]) ./ 255.0
        # plot
        plotBond(from, to, bond_thickness, color)
    end
    # plot all sites
    for index in 1:size(lattice.positions,1)
        p = lattice.positions[index]
        if length(p) == 2
            p = [p[1], p[2], 0.0]
        end
        p = p.*conversion
        plotSite(p, site_radius, get(colorcode_sites, lattice.positions_indices[index], colorcode_sites[0])./255.0)
    end
    # plot all site labels
    for index in 1:size(lattice.positions,1)
        p = lattice.positions[index]
        if length(p) == 2
            p = [p[1], p[2], 0.0]
        end
        p = p.*conversion
        # add label
        if site_labels == "LATTICE INDEX"
            ax[:text](p[1]+site_radius, p[2]+site_radius, p[3]+site_radius,  "$(index)", zorder=1, size=site_radius, color="k")
        elseif site_labels == "POSITION INDEX"
            ax[:text](p[1]+site_radius, p[2]+site_radius, p[3]+site_radius,  "$(lattice.positions_indices[index])", zorder=1, size=site_radius, color="k")
        end
    end


    # turn off everthing but the plot
    ax[:set_axis_off]()
    # get the current limits
    x_limits = ax[:get_xlim3d]()
    y_limits = ax[:get_ylim3d]()
    z_limits = ax[:get_zlim3d]()
    # build the current ranges
    x_range = abs(x_limits[2] - x_limits[1])
    x_middle = (x_limits[2] + x_limits[1]) / 2.0
    y_range = abs(y_limits[2] - y_limits[1])
    y_middle = (y_limits[2] + y_limits[1]) / 2.0
    z_range = abs(z_limits[2] - z_limits[1])
    z_middle = (z_limits[2] + z_limits[1]) / 2.0
    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*maximum([x_range, y_range, z_range])
    ax[:set_xlim3d](x_middle - plot_radius, x_middle + plot_radius)
    ax[:set_ylim3d](y_middle - plot_radius, y_middle + plot_radius)
    ax[:set_zlim3d](z_middle - plot_radius, z_middle + plot_radius)
    # tighten layout
    tight_layout()
    # show
    show()
end
# show Lattice function to plot in mayavi
function showLatticeMayavi(
		lattice::Lattice;
		conversion = 160,
		site_radius=25,
        site_labels="OFF",
		bond_thickness::Int64=8,
		visualize_periodic=false,
		colorcode_sites = Dict(0 => [255,255,255], 1 => [255,255,255]),
		colorcode_bonds = Dict("0" => [0,0,0], "1.0" => [0,0,0]),
        colorcode_bonds_automation::String = "OFF",
        background_color = [200,200,200]
    )

    # CHECK IF MAYAVI IS PRESENT
    @pyimport mayavi.mlab as mlab

    # maybe overwrite dictonary
    if colorcode_bonds_automation == "GREY"
        # construct the list of interaction strengths
        cs_list = getConnectionStrengthList(lattice)
        # get the color code list
        cc_list = getGreySequence(size(cs_list, 1))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[string(cs_list[i])] = cc_list[i]
        end
    elseif colorcode_bonds_automation == "COLOR"
        # construct the list of interaction strengths
        cs_list = getConnectionStrengthList(lattice)
        # get the color code list
        cc_list = getColorSequence(size(cs_list, 1))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[string(cs_list[i])] = cc_list[i]
        end
    end

	# repair color dictonary
	colorcode_bonds["0"] = get(colorcode_bonds, "0", [0,0,0])
	colorcode_sites[0] = get(colorcode_sites, 0, [255,255,255])


    # PLOTTING WITH MAYAVI
    # create a figure
    mlab.figure(bgcolor=(background_color[1], background_color[2], background_color[3]))
    # define a function to plot a sphere for every site
    function plotSite(site,radius,color; detail=10)
        x = site[1]
        y = site[2]
        if length(site) == 3
            z = site[3]
        else
            z = 0.0
        end
        mlab.points3d([x],[y],[z],scale_factor=radius, color=(color[1], color[2], color[3]))
    end
    # define a function to plot a bond as a tube
    function plotBond(from, to, radius, color; detail_ring=8, detail_length=3)
        # check to be 3d
        if length(from) == 2
            from = [from[1], from[2], 0.0]
        end
        if length(to) == 2
            to = [to[1], to[2], 0.0]
        end
        mlab.plot3d([from[1], to[1]], [from[2], to[2]], [from[3], to[3]], color=(color[1], color[2], color[3]), tube_radius=radius)
    end

    # plot all bonds
    for c in lattice.connections
        # skip if wrong directioin
        if c[1] < c[2] || (!visualize_periodic && sum([abs(el) for el in c[4]]) != 0)
            continue
        end
        # get data
        from = lattice.positions[Int(c[1])] .* conversion
        to = lattice.positions[Int(c[2])] .* conversion
        color = get(colorcode_bonds, string(c[3]), colorcode_bonds["0"]) ./ 255.0
        # plot
        plotBond(from, to, bond_thickness, color)
    end
    # plot all sites
    for index in 1:size(lattice.positions,1)
        p = lattice.positions[index]
        if length(p) == 2
            p = [p[1], p[2], 0.0]
        end
        p = p.*conversion
        plotSite(p, site_radius, get(colorcode_sites, lattice.positions_indices[index], colorcode_sites[0])./255.0)
    end
    # plot all site labels
    for index in 1:size(lattice.positions,1)
        p = lattice.positions[index]
        if length(p) == 2
            p = [p[1], p[2], 0.0]
        end
        p = p.*conversion
        # add label
        if site_labels == "LATTICE INDEX"
            mlab.text3d(p[1]+site_radius, p[2]+site_radius, p[3]+site_radius, "$(index)", scale=site_radius, color=(0,0,0))
        elseif site_labels == "POSITION INDEX"
            mlab.text3d(p[1]+site_radius, p[2]+site_radius, p[3]+site_radius, "$(lattice.positions_indices[index])", scale=site_radius, color=(0,0,0))
        end
    end


    # show plot
    mlab.show()


    # return the figure
    return fig
end

































#-----------------------------------------------------------------------------------------------------------------------------
#
#   METHODS FOR JOB SCHEDULING
#
#-----------------------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------------------
#
#   FUNCTION TO CREATE A JOB DIRECTORY FOR SCHEDULING
#
#   USAGE: Give a jobdirectory name and a main program to run as well as a default input file and the function creates
#   a job directory containing all input files as well as two scripts to run the program in parallel or single core to
#   calculate all input files.
#   Usage should be in the following steps:
#   1) provide an input file template of the form: <parameter name>\t<default value>
#   2) provide a main program which can be called by a command similar to "PROGRAM_NAME INPUTFILE.txt"
#   3) call this function to create job directory and all input files
#   4) Enter job directory and run one of the two julia execution scripts
#       - Single core execution one after another: "julia script_single_execution.jl"
#       - MPI based multi-core execution: "mpirun -x PATH --host <host1, host2, host3,...> --map-by node -np <number of processes> nice -n19 julia script_MPI_execution.jl"
#   5) collect output files
#
#
#   Parameters (necessary):
#   - directory_path:       The name of the directory to be created
#   - input_file_template:  The input file template from which parameters are extracted
#   - main_program_file:    The file which contains the main program to be run
#
#   Parameters (optional):
#   - parameters:           A dictonary of parameter values to be used
#   - paramters_use_ARGS:   Bool to determine if the ARGS string is searched for parameters
#   - input_file_seperator: The seperator used for input file seperation
#   - main_program_execution_command:   The command to launch the program PROGRAM with a given input file INPUTFILE. Default: julia PROGRAM INPUTFILE
#   - files_to_copy:        A list of file that are copied to the destination folder as well
#
#-----------------------------------------------------------------------------------------------------------------------------
function createJobDirectory(
    directory_path::String,
    input_file_template_name::String,
    main_program_file::String
    ;
    parameters::Dict=Dict(),
    paramters_use_ARGS::Bool=true,
    input_file_seperator::String=":\t",
    main_program_execution_command::String="julia PROGRAM INPUTFILE",
    files_to_copy::Array{String,1}=String[]
        )


    # -----------------
    # START INFORMATION
    # -----------------

    # correct the directory input
    if !endswith(directory_path, "/")
        directory_path = "$(directory_path)/"
    end
    # Start: Print the directory path
    println("Job directory \"$(directory_path)\" will be created")
    # Print the input template filename
    println("Input template is located at \"$(input_file_template_name)\"\n")




    # ---------------------------------------------------
    # OBTAIN THE INPUT FILE TEMPLATE AND EXTRACT KEYWORDS
    # ---------------------------------------------------

    # Open the template file and read all lines
    input_file_template         = open(input_file_template_name)
    input_file_template_lines   = readlines(input_file_template)
    close(input_file_template)

    # Create a new list for all parameters
    parameters_RAW              = Array[]
    parameters_passed_values    = Array[]

    # push in all relevant lines
    for line in input_file_template_lines
        line = strip(line)
        if startswith(line, "#") || line==""
            continue
        end
        push!(parameters_RAW, split(line, "\t"))
    end

    # print all parameters that are found
    println("Found $(size(parameters_RAW, 1)) parameters in template")
    for p in parameters_RAW
        if length(p) == 1
            println("- \"$(p[1])\" without default value")
        elseif length(p) == 2
            println("- \"$(p[1])\" with default value $(p[2])")
        end
    end
    println("")




    # ---------------------------------------------------
    # OBTAIN THE CORRECT VALUES TO BE USED
    # ---------------------------------------------------

    # Iterate over all accessible parameters
    for p in parameters_RAW

        # list of values to use
        values_to_use   = []

        # OPTION 1, check in Command line arguments
        if paramters_use_ARGS
            # look through all command line arguments
            for arg in ARGS
                # check if it fulfils the structure
                if startswith(arg, "$(p[1])=")
                    # get the string representation
                    arg_value = replace(arg, "$(p[1])=", "")
                    # push to the list of values to use
                    push!(values_to_use, arg_value)
                    # break the loop
                    break
                end
            end
        end

        # OPTION 2, check in dictonary
        if length(values_to_use)==0 && p[1] in keys(parameters)
            # get the string representation
            dict_value = parameters[p[1]]
            # push to the list of values to use
            push!(values_to_use, dict_value)
        end

        # OPTION 3, check default values of template
        if length(values_to_use)==0 && length(p) >= 2
            # get the string representation
            default_value = p[2]
            # push to the list of values to use
            push!(values_to_use, default_value)
        end


        # PARSING OF VALUES

        # create a list with parsed entries
        values_to_use_parsed = []

        # iterate over all values to be used
        for v in values_to_use
            # determine which type of value
            if typeof(v) != String && typeof(v) != SubString{String}
                # already parsed, push into array
                push!(values_to_use_parsed, v)
            elseif startswith(v, "[") && endswith(v, "]")
                # list of format [el1, el2, el3, ...]
                v = split(strip(v[2:end-1]), ",")
                for i in 1:length(v)
                    v[i] = strip(v[i])
                    push!(values_to_use_parsed, parse(Float64,v[i]))
                end
            elseif startswith(v, "linspace(") && endswith(v, ")")
                # list via linspace(start, end, steps)
                v = split(v[10:end-1], ",")
                for i in 1:length(v)
                    v[i] = strip(v[i])
                end
                v = linspace(parse(Float64,v[1]),parse(Float64,v[2]),parse(Float64,v[3]))
                for i in 1:length(v)
                    push!(values_to_use_parsed, v[i])
                end
            else
                # it is probably either a number or a string, try parsing to float
                if isnull(tryparse(Float64,v))
                    # not a float, just push in list as it is
                    push!(values_to_use_parsed, v)
                else
                    # push parsed into list
                    push!(values_to_use_parsed, parse(Float64,v))
                end
            end
        end


        # append the list of paramters to the global list
        push!(parameters_passed_values, values_to_use_parsed)

    end

    # print all values that are found
    println("Found values for parameters:")
    for (p_index,p) in enumerate(parameters_RAW)
        println("- \"$(p[1])\", $(length(parameters_passed_values[p_index])) values: $(parameters_passed_values[p_index])")
    end
    println("")






    # -----------
    # CREATE JOBS
    # -----------

    # create a new list of dictonaries which are the jobs
    joblist = []
    # add a dummy job without information
    push!(joblist, Dict())

    # iterate over all parameters
    for p in 1:length(parameters_RAW)

        # find out if single value or multiple values
        if length(parameters_passed_values[p]) == 1
            # only one value, all jobs get this value
            for j in joblist
                j["$(parameters_RAW[p][1])"] = parameters_passed_values[p][1]
            end
        elseif length(parameters_passed_values[p]) >= 1
            # more than one value, every job gets copied for every value
            joblist_backup = joblist
            joblist = []
            for j in joblist_backup
                for i in 1:length(parameters_passed_values[p])
                    j_tmp = copy(j)
                    j_tmp["$(parameters_RAW[p][1])"] = parameters_passed_values[p][i]
                    push!(joblist, j_tmp)
                end
            end
        else
            # PROBLEM: NO VALUE FOR PARAMETER, exlude from job
            println("DANGER: Excluding parameter $(parameters_RAW[p][1]) from jobs because it has no value")
        end

    end

    # print the success
    println("Job building is done, $(length(joblist)) jobs are constructed\n")

    # check if only one job without parameters
    if length(joblist) == 1 && length(keys(joblist[1])) == 0
        println("Some mistake occured, only job has no parameters. Aborting now...\n")
        return
    end





    # --------------------
    # CREATE THE DIRECTORY
    # --------------------

    # check if the directory exists
    if isdir(directory_path)
        # Directory exists
        println("Directory \"$(directory_path)\" already exists")
        print("Delete and reconstruct? [y/n]: ")
        answer = strip(readline()[end-1:end])
        while !in(answer, ["y", "n"])
            print("Delete and reconstruct? [y/n]: ")
            answer = strip(readline()[end-1:end])
        end
        # depending on answer
        if answer == "y"
            # ask for recursive deletion
            print("Delete recursively? [y/n]: ")
            answer = strip(readline()[end-1:end])
            while !in(answer, ["y", "n"])
                print("Delete and reconstruct? [y/n]: ")
                answer = strip(readline()[end-1:end])
            end
            if answer == "y"
                # delete the complete directory
                rm(directory_path, recursive=true)
            else
                # delete the complete directory
                rm(directory_path, recursive=false)
            end
            # print
            println("Deleted already existing directory $(directory_path)")
            # create the directory
            mkdir(directory_path)
            println("Recreated directory $(directory_path)\n")

        else
            # abort the program
            println("Not deleting directory, aborting program now...")
            return
        end
    else
        # Directory does not exist, create it
        mkdir(directory_path)
        println("Created directory $(directory_path)\n")
    end




    # ----------------------
    # CREATE THE INPUT FILES
    # ----------------------

    # iterate over all jobs
    for i in 1:length(joblist)
        # dump job i to file input_i.txt
        jobfilename = "input_$(i).txt"
        f = open("$(directory_path)$(jobfilename)", "w")
        # write entries of the dictonary into the file
        for k in keys(joblist[i])
            write(f, "$(k)$(input_file_seperator)$(joblist[i][k])\n")
        end
        # close the file
        close(f)
        # set the correct input filename in the job
        joblist[i]["input_filename"] = jobfilename
    end
    # print the success
    println("Created all $(length(joblist)) jobfiles\n")




    # ------------------------------
    # COPY / CREATE THE SCRIPT FILES
    # ------------------------------

    # Copy main program file
    cp(main_program_file, "$(directory_path)$(main_program_file)")
    println("Copied main program file \"$(main_program_file)\"")

    # Create single execution julia script
    script_single_exe = open("$(directory_path)script_single_execution.jl", "w")
    # write all lines
    for i in 1:length(joblist)
        # compose a line
        line = "run(`$(main_program_execution_command)`)\n"
        # change the program part
        line = replace(line, "PROGRAM", main_program_file)
        # change the input part
        line = replace(line, "INPUTFILE", joblist[i]["input_filename"])
        # set the correct input filename in the job
        write(script_single_exe, line)
    end
    # close the script
    close(script_single_exe)
    println("Created single execution script \"script_single_execution.jl\"")


    # Create single execution julia script
    script_MPI_exe = open("$(directory_path)script_MPI_execution.jl", "w")


    # write all lines, starting with imports
    write(script_MPI_exe, "# use with:\n")
    write(script_MPI_exe, "# mpirun -x PATH --host <host1,host2,...> --map-by node -np <number of processes> nice -n19 julia script_MPI_execution.jl\n")
    write(script_MPI_exe, "using MPI\n")
    write(script_MPI_exe, "MPI.Init()\n")
    # get the number of the MPI process
    write(script_MPI_exe, "commworld = MPI.COMM_WORLD\n")
    write(script_MPI_exe, "commrank  = MPI.Comm_rank(commworld)\n")
    write(script_MPI_exe, "commsize  = MPI.Comm_size(commworld)\n")
    write(script_MPI_exe, "max_i     = $(length(joblist))\n")
    write(script_MPI_exe, "init_i    = commrank+1\n")
    write(script_MPI_exe, "for i in init_i:commsize:max_i\n")
    # compose main line
    line = "\trun(pipeline(`$(main_program_execution_command)`, stdout=\"bash_stdout_INPUTFILE\", stderr=\"bash_stderr_INPUTFILE\"))\n\tprintln(\"Done input file \\\"INPUTFILE\\\"\")\n"
    # change the program part
    line = replace(line, "PROGRAM", main_program_file)
    # change the input part
    line = replace(line, "INPUTFILE", "input_\$(i).txt")
    # set the correct input filename in the job
    write(script_MPI_exe, line)
    write(script_MPI_exe, "end\n")
    # Finalize MPI script
    write(script_MPI_exe, "MPI.Finalize()\n")
    # close the script
    close(script_MPI_exe)
    println("Created MPI execution script \"script_MPI_execution.jl\"")

    # Copy additional files
    for fn in files_to_copy
        cp(fn, "$(directory_path)$(fn)")
        println("Copied additional file \"$(fn)\"")
    end

end
export createJobDirectory
