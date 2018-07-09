# function to dump lattice data to a file for blender to read
function dumpBlenderFile(
            lattice::Lattice,
            filename_output::String="AUTO";
            visualize_periodic::Bool=false,
            site_radius::Float64=-0.2,
            bond_thickness::Float64=-0.1,
            print_used_options::Bool=false,
            site_labels::String="OFF",
    		colorcode_sites::Dict = Dict(0 => [255,255,255]),
    		colorcode_bonds::Dict = Dict("0" => [0,0,0]),
            colorcode_bonds_automation::String = "OFF"
        )

    ##########---------------------------------
    # STEP 1 #  INITIALIZATION OF PARAMTERS
    ##########---------------------------------

    # Maybe print what is being done
    if print_used_options
        println("Dumping lattice object with filename \"$(lattice.filename)\" to blender file")
    end

    # define the FILENAME of the output file if it is set to "AUTO"
    if filename_output=="AUTO"
        filename_output = "$(lattice.filename[1:end-4])_blender_dump.txt"
        filename_output = replace(filename_output, FOLDER_LATTICES, "")
        # maybe print
        if print_used_options
            println("- Changing Blender filename to \"$(filename_output)\"")
        end
    elseif print_used_options
        println("- Using given SVG filename \"$(filename_output)\"")
    end






    ############---------------------------------
    # STEP 2.1 #  LOADING OF CONNECTION DATA
    ############---------------------------------

	# connections to plot (all)
	connections     = deepcopy(lattice.connections)
	positions	    = deepcopy(lattice.positions)
    indices_to_plot = deepcopy(lattice.positions_indices)

    # Function to determine if a connection should be plotted
    # 0 yes, no periodic
    # 1 yes, periodic
    # -1 no
    function shallBePlotted(c::Array{Any,1})
        # check if the indices are okay (exclude double bonds)'
        if c[1] > c[2]
            return -1
        end
        # check the sum of all wraps
        if sum(abs.(c[4])) > 0
            if visualize_periodic
                return 1
            else
                return -1
            end
        else
            return 0
        end
    end

    # print plotting of connections
    if print_used_options
        if visualize_periodic
            println("- Periodic connections will be plotted, based on visualize_periodic=$(visualize_periodic)")
        else
            println("- Periodic connections will not be plotted, based on visualize_periodic=$(visualize_periodic)")
        end
    end




    ############-------------------------------------
    # STEP 2.2 #  DETERMINATION OF SITE & BOND DATA
    ############-------------------------------------


    ##############--------------
    # STEP 2.2.1 #  SITE DATA
    ##############--------------

    # check if site radius has to be chosen
    if site_radius < 0
        # site radius is given in relative units, first maybe print
        if print_used_options
            println("- Site radius is given in relative units to be site_radius=$(site_radius)")
        end
        # then, calculate the minimum length of connections
        min_length = 999999999
        for c in connections
            # check if the connection is NOT PERIODIC and NOT FROM SITE TO ITSELF
            if sum(abs.(c[4]))==0 && c[1]!=c[2]
                # determine the length of the position on the canvas
                c_vector = positions[Int(c[1])] .- positions[Int(c[2])]
                c_length = sqrt(sum(c_vector .* c_vector))
                # check if smaller then minimum length
                min_length = min(min_length, c_length)
            end
        end
        # set the radius of a site to be -site_radius*min_length
        site_radius = -site_radius*min_length
        # maybe print
        if print_used_options
            println("- Site radius changed to absolute units to be site_radius=$(site_radius) as minimal bond length is $(min_length)")
        end
    else
        # site radius is given in absolute units, maybe print
        if print_used_options
            println("- Site radius is given in absolute units to be site_radius=$(site_radius)")
        end
    end


    # repair color dictonary by making sure default elements exist
    colorcode_sites[0]   = get(colorcode_sites, 0, [255,255,255])
    # maybe print
    if print_used_options
        println("- Site colors passed with dictonary:")
        for k in keys(colorcode_sites)
            println("    - $(k) -> $(colorcode_sites[k])")
        end
    end




    ##############--------------
    # STEP 2.2.2 #  BOND DATA
    ##############--------------

    # check if bond thickness has to be chosen
    if bond_thickness < 0
        # bond thickness is given in relative units, first maybe print
        if print_used_options
            println("- Bond thickness is given in relative units to be bond_thickness=$(bond_thickness)")
        end
        # then, calculate the minimum length of connections
        min_length = 999999999
        for c in connections
            # check if the connection is NOT PERIODIC and NOT FROM SITE TO ITSELF
            if sum(abs.(c[4]))==0 && c[1]!=c[2]
                # determine the length of the position on the canvas
                c_vector = positions[Int(c[1])] .- positions[Int(c[2])]
                c_length = sqrt(sum(c_vector .* c_vector))
                # check if smaller then minimum length
                min_length = min(min_length, c_length)
            end
        end
        # set the thickness of a bond to be -bond_thickness*min_length
        bond_thickness = -bond_thickness*min_length
        # maybe print
        if print_used_options
            println("- Bond thickness changed to absolute units to be bond_thickness=$(bond_thickness) as minimal bond length is $(min_length)")
        end
    else
        # bond thickness is given in absolute units, maybe print
        if print_used_options
            println("- Bond thickness is given in absolute units to be bond_thickness=$(bond_thickness)")
        end
    end


    # HANDLING BOND COLORS

    # get list of all connections strengths
    cs_list = getConnectionStrengthList(lattice)
    # parse all elements to strings
    cs_list = String["$(element)" for element in cs_list]
    # sort the list
    sort!(cs_list)

    # BOND COLOR DICTONARY
    # maybe overwrite dictonary
    if colorcode_bonds_automation == "GREY"
        # get the color code list
        cc_list = LatticePhysics.getRandomGreys(length(cs_list))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[cs_list[i]] = cc_list[i]
        end
        # maybe print
        if print_used_options
            println("- Using random GREY for bonds, $(length(cc_list)) greys used")
        end
    elseif colorcode_bonds_automation == "COLOR"
        # get the color code list
        cc_list = LatticePhysics.getRandomColors(length(cs_list))
        # put in a new dictonary
        colorcode_bonds = Dict()
        # insert all pairs
        for i in 1:size(cc_list, 1)
            colorcode_bonds[cs_list[i]] = cc_list[i]
        end
        # maybe print
        if print_used_options
            println("- Using random COLOR for bonds, $(length(cc_list)) colors used")
        end
    elseif colorcode_bonds_automation == "OFF"
        # maybe print
        if print_used_options
            println("- No redefinition of bond colors, using provided dictonary:")
            for k in keys(colorcode_bonds)
                println("    - \"$(k)\" -> $(colorcode_bonds[k])")
            end
        end
    else
        # maybe print
        if print_used_options
            println("- Could not identify redefinition of bond colors based on colorcode_bonds_automation=\"$(colorcode_bonds_automation)\", using provided dictonary:")
            for k in keys(colorcode_bonds)
                println("    - \"$(k)\" -> $(colorcode_bonds[k])")
            end
        end
    end

	# repair color dictonary by making sure default elements exist
	colorcode_bonds["0"] = get(colorcode_bonds, "0", [0,0,0])







    ##########---------------------
    # STEP 3 #  BLENDER FILE CREATION
    ##########---------------------


    # open the file in write mode
    f = open(filename_output, "w")

    # write all materials (bond / site color directories)
    for k in keys(colorcode_sites)
        # write start of the line
        write(f, "material:\t")
        # write the name
        write(f, "site_material_$(k)\t")
        # write the color
        site_color = colorcode_sites[k]
        write(f, "$(site_color[1]), $(site_color[2]), $(site_color[3])")
        # write line break
        write(f, "\n")
    end
    for k in keys(colorcode_bonds)
        # write start of the line
        write(f, "material:\t")
        # write the name
        write(f, "bond_material_$(k)\t")
        # write the color
        connection_color = colorcode_bonds[k]
        write(f, "$(connection_color[1]), $(connection_color[2]), $(connection_color[3])")
        # write line break
        write(f, "\n")
    end

    # write all sites
    for (i,p) in enumerate(positions)
        # write start of the line
        write(f, "site:\t$(i), ")
        # write the position
        if length(p) == 2
            write(f, "$(p[1]), $(p[2]), 0.0")
        elseif length(p) == 3
            write(f, "$(p[1]), $(p[2]), $(p[3])")
        else
            println("site has not correct amount of dimensions!")
        end
        # write the radius
        write(f, ", $(site_radius)")
        # write the color
        if indices_to_plot[i] in keys(colorcode_sites)
            write(f, "\tsite_material_$(indices_to_plot[i])")
        else
            write(f, "\tsite_material_$(0)")
        end
        # write line break
        write(f, "\n")
    end

    # write all bonds
	for (i,c) in enumerate(connections)
        # non periodic or periodic and wanted to be plotted
		if shallBePlotted(c) >= 0
            # write start of the line
            write(f, "bond:\t$(i), ")
            # write the position from where
            p = positions[Int(c[1])]
            if length(p) == 2
                write(f, "$(p[1]), $(p[2]), 0.0")
            elseif length(p) == 3
                write(f, "$(p[1]), $(p[2]), $(p[3])")
            else
                println("from site has not correct amount of dimensions!")
            end
            # write the position to where
            p = positions[Int(c[2])]
            if length(p) == 2
                write(f, ", $(p[1]), $(p[2]), 0.0")
            elseif length(p) == 3
                write(f, ", $(p[1]), $(p[2]), $(p[3])")
            else
                println("from site has not correct amount of dimensions!")
            end
            # write the thickness
            write(f, ", $(bond_thickness)")
            # write the color
            if string(c[3]) in keys(colorcode_bonds)
                write(f, "\tbond_material_$(string(c[3]))")
            else
                write(f, "\tbond_material_0")
            end
            # write line break
            write(f, "\n")
        end
	end

    # close the file
    close(f)


    ##########---------------------
    # STEP 4 #  FILE FINISH
    ##########---------------------


    # return the filename
    return filename_output

end
