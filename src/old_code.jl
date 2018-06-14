
#-----------------------------------------------------------------------------------------------------------------------------
#
#   Find Plaquettes of some site in the lattice with desired length len
#
#-----------------------------------------------------------------------------------------------------------------------------
function getPlaquettesOfSite(lattice::Lattice, site, len)
    # get the connectivity of the lattice
    cl = getConnectionList(lattice)
    # history includes all sites that are already visited, including the current site
    function givePlaquettePaths(history, steps, patharray)
        # check if steps = 0, if yes, return history only
        if length(history) == steps+1
            if history[end] == history[1]
                push!(patharray, history)
            end
            return
        elseif history[end] in history[1:end-1]
            # building a closed loop
            return
        end
        # current site
        currentSite = history[end]
        # get all options
        options = []
        for c in cl[currentSite]
            push!(options, Int(c[2]))
        end
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
    paths = Array[]
    givePlaquettePaths([site], len, paths)
    # check all paths if they are plaquettes
    plaquettes = Array[]
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

function getPlaquettesOfLattice(lattice::Lattice, len)
    # get all plaquettes of all sites
    plaquettes = Array[]
    plaquettes_sorted = Array[]
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



function printPlaquetteStatisticsOfLattice(lattice::Lattice; detailed=false, max_length=12)
    # print header information
    println("printing plaquette information for lattice:")
    # iterate over all lengths
    for l in 3:max_length
        # get plaquettes
        plaquettes = getPlaquettesOfLattice(lattice, l)
        # print statistics
        println("- length $(l) bonds: $(size(plaquettes,1)) plaquettes")
        # print more information
        if detailed
            for p in plaquettes
                println("    - $(p)")
            end
        end
    end
end
export printPlaquetteStatisticsOfLattice

function printPlaquetteStatisticsOfSite(lattice::Lattice, site; detailed=false, max_length=12)
    # print header information
    println("printing plaquette information for site $(site) of lattice:")
    # iterate over all lengths
    for l in 3:max_length
        # get plaquettes
        plaquettes = getPlaquettesOfSite(lattice, site, l)
        # print statistics
        println("- length $(l) bonds: $(size(plaquettes,1)) plaquettes")
        # print more information
        if detailed
            for p in plaquettes
                println("    - $(p)")
            end
        end
    end
end
export printPlaquetteStatisticsOfSite

















#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#
#   PLOTTING TO SVG IMAGES
#
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------------------------------
#
#   HELPER METHODS FOR SVG CREATION
#   Get the correct SVG code lines for SVG beginning and end of document
#   Get the correct SVG code lines for certain objects
#
#-----------------------------------------------------------------------------------------------------------------------------

# HEADER STRING (must be first in every SVG file)
# width and height denote the dimensions of the image in px
function getSVGHeaderString(width::Int64, height::Int64)
	headerstring = """<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>
<svg
	xmlns:dc=\"http://purl.org/dc/elements/1.1/\"
	xmlns:cc=\"http://creativecommons.org/ns#\"
	xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"
	xmlns:svg=\"http://www.w3.org/2000/svg\"
	xmlns=\"http://www.w3.org/2000/svg\"
	xmlns:xlink=\"http://www.w3.org/1999/xlink\"
	version=\"1.1\"
	id=\"svg2\"
	viewBox=\"0 0 $(width) $(height)\"
	height=\"$(height)\"
	width=\"$(width)\">
	<defs id=\"defs4\">
	</defs>
	<metadata
		id=\"metadata7\">
		<rdf:RDF>
			<cc:Work rdf:about=\"\">
				<dc:format>image/svg+xml</dc:format>
				<dc:type
					rdf:resource=\"http://purl.org/dc/dcmitype/StillImage\" />
				<dc:title></dc:title>
			</cc:Work>
		</rdf:RDF>
	</metadata>

"""
	return headerstring
end

# FOOTER STRING (must be end of every SVG file)
function getSVGFooterString()
    return """
</svg>
"""
end

# STRING FOR AN ELLIPSE
# Parameters are
# - id: The id which the object has later on
# - centerX / centerY: The coordinates of the ellipse center in the canvas
# - radiusX / radiusY: The radii of the ellipse
# - color: hex string of the ellipse color
# formerly: getEllipseString
function getSVGStringEllipse(id, centerX, centerY, radiusX, radiusY, color; label=987654321, labelcolor="#000000")
	es = """
	<ellipse
		style=\"color:$(color);fill:$(color);fill-opacity:1;fill-rule:nonzero\"
		id=\"$(id)\"
		cx=\"$(centerX)\"
		cy=\"$(centerY)\"
		rx=\"$(radiusX)px\"
		ry=\"$(radiusY)px\" />

"""
    # check if to add a label
    if label != 987654321
        es = """
    $(es)

    <text
    x=\"$(centerX)\"
    y=\"$(centerY+radiusY/2)\"
    style=\"text-anchor: middle; font-size: $(radiusX*1.2)px; fill:$(labelcolor)\"
    >
    $(label)
    </text>
"""
    end
	return es
end

# STRING FOR A STROKED ELLIPSE
# Parameters are
# - id: The id which the object has later on
# - centerX / centerY: The coordinates of the ellipse center in the canvas
# - radiusX / radiusY: The radii of the ellipse
# - colorFill: hex string of the ellipse fill color
# - colorStroke: hex string of the ellipse stroke color
# - strokewidth: Float or Int giving the width of the surrounding stroke
# formerly: getStrokedEllipseString
function getSVGStringEllipseStroked(id, centerX, centerY, radiusX, radiusY, colorFill, colorStroke, strokewidth; label=987654321, labelcolor="#000000")
	es = """
	<ellipse
		style=\"color:$(colorFill);fill:$(colorFill);fill-opacity:1;fill-rule:nonzero;stroke:$(colorStroke);stroke-width:$(strokewidth);stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-dasharray:none;stroke-dashoffset:0;stroke-opacity:1\"
		id=\"$(id)\"
		cx=\"$(centerX)\"
		cy=\"$(centerY)\"
		rx=\"$(radiusX)px\"
		ry=\"$(radiusY)px\" />

"""
    # check if to add a label
    if label != 987654321
        es = """
    $(es)

    <text
    x=\"$(centerX)\"
    y=\"$(centerY+radiusY/2)\"
    style=\"text-anchor: middle; font-size: $(radiusX*1.2)px; fill:$(labelcolor)\"
    >
    $(label)
    </text>
"""
    end
	return es
end

# STRING FOR A LINE
# Parameters are
# - id: The id which the object has later on
# - from: coordinate array [x,y] of the starting point
# - to:   coordinate array [x,y] of the end point
# - colorStroke: hex string of the color that the line has
# - strokewidth: Float or Int giving the width of the line
# - dashed: is the line dashed or not
# formerly: getLineString
function getSVGStringLine(id, from, to, colorStroke, strokewidth; dashed=false, opacity::Float64=1.0)
	if dashed
        ls = """
	<path
		style=\"fill:none;fill-rule:evenodd;stroke:$(colorStroke);stroke-width:$(strokewidth);stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-linecap:butt;stroke-dasharray:$(strokewidth),$(strokewidth*2);stroke-dashoffset:0;stroke-opacity:$(opacity)\"
		id=\"$(id)\"
		d=\"m $(from[1]),$(from[2]) $(to[1]-from[1]),$(to[2]-from[2])\"/>
"""
    else
        ls = """
	<path
		style=\"fill:none;fill-rule:evenodd;stroke:$(colorStroke);stroke-width:$(strokewidth);stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-dasharray:none;stroke-dashoffset:0;stroke-opacity:$(opacity)\"
		id=\"$(id)\"
		d=\"m $(from[1]),$(from[2]) $(to[1]-from[1]),$(to[2]-from[2])\"/>
"""
    end
	return ls
end


# STRING FOR A PLAQUETTE
# Parameters are
# - id: The id which the object has later on
# - plaquette_points: list of point coordinates of points that form the plaquette
# - colorFill: hex string of the color that the plaquette has
# - opacity: Opacity of the plaquette
function getSVGStringPlaquette(id, plaquette_points, colorFill; opacity=0.5)
    plaquette_string = "M"
    for p in plaquette_points
        plaquette_string = "$(plaquette_string) $(p[1]),$(p[2])"
    end
    plaquette_string = "$(plaquette_string) z"
    ps = """
	<path
		style=\"fill:$(colorFill);fill-opacity:$(opacity);fill-rule:evenodd;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"
		d=\"$(plaquette_string)\"
		id=\"$(id)\" />
"""
end

# CONVERSION OF RGB COLORS TO HEX STRINGS
function color_hex(r::Int64,g::Int64,b::Int64)
    return "#$(hex(r,2))$(hex(g,2))$(hex(b,2))"
end
function color_hex(rgb::Array{Int64})
    return color_hex(rgb[1], rgb[2], rgb[3])
end



# GET AUTOMATICAL SEQUENCE OF COLORS
function getColorSequence(len::Int64)
    # define a new list for the sequence
    colors = Array[]
    # fill the sequence depending on the number of requested elements
    if len <= 1
        # only one color --> black
        push!(colors, [0,0,0])
    elseif len == 2
        # only two colors --> black and red
        push!(colors, [0,0,0])
        push!(colors, [255,0,0])
    elseif len == 3
        # only three colors --> red, green, blue (Kitaev)
        push!(colors, [255,0,0])
        push!(colors, [0,255,0])
        push!(colors, [0,0,255])
    else
        # continuous sequence of random colors
        for i in 1:len
            push!(colors, [rand(1:255),rand(1:255),rand(1:255)])
        end
    end
    # return the sequence
    return colors
end
function getGreySequence(len::Int64)
    # define a new list for the sequence
    colors = Array[]
    # fill the sequence depending on the number of requested elements
    if len <= 1
        # only one color --> black
        push!(colors, [0,0,0])
    elseif len == 2
        # only two colors --> black and grey
        push!(colors, [0,0,0])
        push!(colors, [125,125,125])
    elseif len == 3
        # only three colors --> black and two grey
        push!(colors, [0,0,0])
        push!(colors, [90,90,90])
        push!(colors, [180,180,180])
    else
        # continuous sequence
        for i in 1:len
            c = round(Int64, i/len * 200)
            push!(colors, [c,c,c])
        end
    end
    # return the sequence
    return colors
end

















#-----------------------------------------------------------------------------------------------------------------------------
#
#   SVG PLOTTING OF 2D AND 3D LATTICES
#
#   Generates a SVG file that shows the lattice
#   SVG file can optionally be opened and converted to pdf
#
#   Parameters (necessary 2D and 3D):
#   -   lattice: The lattice object to plot
#
#   Parameters (optional 2D and 3D):
#   -   conversion: The factor of how many pixels represent a length 1 of real space distance (Default: 160)
#   -   border_percentage: percentage of the width which is allocated as border (Default: 0.1)
#   -   filename_output: The output filename, "AUTO" meaning it is generated automatically (Default: "AUTO")
#   -   site_radius: Radius of sites (Default: 25)
#   -   site_border_width_percentage: Border width of sites in percent of the diameter (Default: 0.2)
#   -   bond_thickness: The thickness of bonds (Default: 8)
#   -   visualize_periodic: Decide if periodic connections are also shown (Default: false)
#   -   colorcode_sites: Dictonary containing a colorcode for all sites matching the site index to a RGB color array
#   -   colorcode_bonds: Dictonary containing a colorcode for all bonds matching the interaction strength to a RGB color array
#   -   colorcode_bonds_automation: String specifying if bonds are colored automatically. Select either "COLOR" or "GREY" or "OFF"
#   -   openfile: Decide if the newly created SVG should already be opened by external viewer (Default: false)
#   -   export_pdf: Decide if svg should be converted to pdf (Default: true)
#
#   Parameters (optional but only for 3D):
#   -   bond_color_BG: The color of bonds that are in the farthest back of the plot (i.e. to which color to converge towards the back)
#   -   site_color_BG: The color of sites that are in the farthest back of the plot (i.e. to which color to converge towards the back)
#   -   lattice_rotation: Array of Floats indicating the rotation around the XYZ axis
#   -   camera_position_relative: Array of Floats indicating the position of the camera
#   -   DOF: Float indicating the strength of the Depth of Field in color gradient
#
#-----------------------------------------------------------------------------------------------------------------------------
function plotLattice2D(
		lattice::Lattice;
		conversion = 160,
		border_percentage=0.1,
		filename_output::String="AUTO",
		site_radius=25,
		site_border_width_percentage::Float64=0.2,
        site_labels="OFF",
		bond_thickness::Int64=8,
		visualize_periodic=false,
		colorcode_sites = Dict(0 => [255,255,255], 1 => [255,255,255]),
		colorcode_bonds = Dict("0" => [0,0,0], "1.0" => [0,0,0]),
        colorcode_bonds_automation::String = "OFF",
		openfile=false,
        export_pdf=true
		)

    # define the filename_output if it is set to AUTO
    if filename_output=="AUTO"
        filename_output = "$(lattice.filename[1:end-4])_plot.svg"
        filename_output = replace(filename_output, FOLDER_LATTICES, "")
    end

	# load positions and connections
	positions	= lattice.positions
	connections = lattice.connections

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


	# define styles for the different sites, i.e. strings that are saved into the svg strings
	site_r 				= site_radius
	site_border			= "#000000"
	site_border_width	= "$(site_border_width_percentage*site_radius)px"

	# sites to plot
	pos_x = zeros(size(positions,1))
	pos_y = zeros(size(positions,1))
	for i in 1:size(positions,1)
		pos_x[i] = positions[i][1]
		pos_y[i] = positions[i][2]
	end
	sites_to_plot = positions
	indices_to_plot = lattice.positions_indices
	xvals = pos_x
	yvals = pos_y
	border 			= border_percentage * (maximum(pos_x) + maximum(pos_y) - minimum(pos_x) - minimum(pos_y))/4

	# connections to plot (all)
	connections_to_plot = copy(connections)
    # the neutral connection wrap, i.e. which wrap identifies a non-periodic connection
	neutral_connection_wrap = (0,0)
	if size(lattice.lattice_vectors,1) == 1
		neutral_connection_wrap = (0)
	end


	# define the width and height of the canvas
	width_UC 	= (maximum(xvals) + site_radius/conversion + border) - (minimum(xvals) - site_radius/conversion - border)
	height_UC	= (maximum(yvals) + site_radius/conversion + border) - (minimum(yvals) - site_radius/conversion - border)
	min_x	= minimum(xvals) - site_radius/conversion - border
	min_y	= minimum(yvals) - site_radius/conversion - border
	width	= conversion*(width_UC)
	height	= conversion*(height_UC)

	# define the conversion functions for coordinates
	function X(x)
		return conversion*(x - min_x)
	end
	function Y(y)
		return + conversion*(y - min_y)
	end

	# open the SVG file
	file = open(filename_output, "w")

	# write the headerstring
	write(file, getSVGHeaderString(round(Int64, width), round(Int64, height)))

	# write all connections
	for (i,c) in enumerate(connections_to_plot)
		if (c[4] == neutral_connection_wrap)
			connection_color = color_hex(get(colorcode_bonds, string(c[3]), colorcode_bonds["0"]))
			write(file, getSVGStringLine("path$(i)", [X(positions[Int(c[1])][1]), Y(positions[Int(c[1])][2])], [X(positions[Int(c[2])][1]), Y(positions[Int(c[2])][2])], connection_color, bond_thickness))
		elseif visualize_periodic
			connection_color = color_hex(get(colorcode_bonds, string(c[3]), colorcode_bonds["0"]))
			write(file, getSVGStringLine("path$(i)", [X(positions[Int(c[1])][1]), Y(positions[Int(c[1])][2])], [X(positions[Int(c[2])][1]), Y(positions[Int(c[2])][2])], connection_color, bond_thickness, dashed=true))
        end
	end

	# write all sites
	for (i,s) in enumerate(sites_to_plot)
        site_color_basic = get(colorcode_sites, indices_to_plot[i], colorcode_sites[0])
		site_color = color_hex(site_color_basic)
        label_color = color_hex([colorelment < 100 ? 255 : 0 for colorelment in site_color_basic])
        if site_labels == "POSITION INDEX"
		    write(file, getSVGStringEllipseStroked("el$(i)", X(s[1]), Y(s[2]), site_r, site_r, site_color, site_border, site_border_width, label=indices_to_plot[i], labelcolor=label_color))
        elseif site_labels == "LATTICE INDEX"
		    write(file, getSVGStringEllipseStroked("el$(i)", X(s[1]), Y(s[2]), site_r, site_r, site_color, site_border, site_border_width, label=i, labelcolor=label_color))
        else
		    write(file, getSVGStringEllipseStroked("el$(i)", X(s[1]), Y(s[2]), site_r, site_r, site_color, site_border, site_border_width))
        end
	end

	# write the footerstring to close the svg file
	write(file, getSVGFooterString())
	# close the file
	close(file)

	# convert to pdf
    if export_pdf
	    run(`inkscape $(filename_output) --export-pdf $(filename_output[1:end-4]).pdf`)
    end
	# if file shall be opened
	if openfile
        if is_linux()
		    run(`ristretto $(filename_output)`)
        elseif is_windows()
            try
                if export_pdf
                    run(`explorer $(filename_output[1:end-4]).pdf`)
                else
                    run(`explorer $(filename_output)`)
                end
            end
        elseif is_apple()
            run(`open $(filename_output)`)
        else
            println("wanted to show file but operating system not supported for showing file")
        end
	end

	# return the output filename
	return filename_output

end
export plotLattice2D


function plotLattice3D(
		lattice::Lattice;
		border_percentage=0.1,
		filename_output::String="AUTO",
		site_radius=25,
		site_border_width_percentage::Float64=0.2,
        site_labels="OFF",
		bond_thickness::Int64=8,
		visualize_periodic=false,
		colorcode_sites = Dict(0 => [255,255,255], 1 => [255,255,255]),
		colorcode_bonds = Dict("0" => [0,0,0], "1.0" => [0,0,0]),
        colorcode_bonds_automation::String = "OFF",
		bond_color_BG = [100,100,100],
		site_color_BG = [125,125,125],
		lattice_rotation::Array{Float64}=[0.0,0.0,0.0],
		camera_position_relative::Array{Float64}=[0.0,0.0,5.0],
		DOF::Float64 = 0.8,
		openfile=false,
        export_pdf=true
		)

    # define the filename_output if it is set to AUTO
    if filename_output=="AUTO"
        filename_output = "$(lattice.filename[1:end-4])_plot.svg"
        filename_output = replace(filename_output, FOLDER_LATTICES, "")
    end

	# load positions and connections
	positions = copy(lattice.positions)
	connections = copy(lattice.connections)

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

	# shift the positions so they are centered in (0,0,0)
	position_center = [0,0,0]
	for pos in positions
		position_center = position_center .+ pos
	end
	position_center = position_center ./ size(positions,1)
	for index in 1:size(positions,1)
		positions[index] = positions[index] .- position_center
	end

	# rotate around all axis
	for (index,angle_degrees) in enumerate(lattice_rotation)
		# get the angle in radians
		angle = angle_degrees *pi/180.0
		# define the rotation matrix
		R = eye(3) * cos(angle)
		R[index, index] = 1
		R[(index+1)%3+1,(index)%3+1] = sin(angle)
		R[(index)%3+1,(index+1)%3+1] = -sin(angle)
		#println(R)
		# rotate the positions
		for index in 1:size(positions,1)
			positions[index] = R*positions[index]
		end
	end


	# vector absolute values
	function vecabs(r)
		return sqrt(r[1]*r[1] + r[2]*r[2] + r[3]*r[3])
	end

	# determine the position on screen
	# with reference to https://en.wikipedia.org/wiki/3D_projection

	# get the arrays of screen positions
	b_x = zeros(size(positions, 1))
	b_y = zeros(size(positions, 1))
	# the array of distance to camera
	distances_to_camera = zeros(size(positions,1))

	# determine the constants to be used

	# camera position in xyz coordinates
	max_x = maximum([abs(pos[1]) for pos in positions])
	max_y = maximum([abs(pos[2]) for pos in positions])
	max_z = maximum([abs(pos[3]) for pos in positions])
	c_x = max_x * camera_position_relative[1]
	c_y = max_y * camera_position_relative[2]
	c_z = max_z * camera_position_relative[3]
	c_vec = [c_x, c_y, c_z]

	# camera angles
	theta_x = 2*pi * (0)
	theta_y = 2*pi * (0)
	theta_z = 2*pi * (0)

	# viewer position
	e_x = 0
	e_y = 0
	e_z = 1

	# calculate the projection of each points
	for index in 1:length(b_x)
		# determine the vector x
		x = positions[index][1] - c_x
		y = positions[index][2] - c_y
		z = positions[index][3] - c_z
		# calculate the distance to the camera of all positions
		distances_to_camera[index] = vecabs(positions[index] .- c_vec)
		# determine the displacements
		d_x = cos(theta_y)*(sin(theta_z)*y + cos(theta_z)*x)  -  sin(theta_y)*(z)
		d_y = sin(theta_x)*(cos(theta_y)*z + sin(theta_y)*(sin(theta_z)*y + cos(theta_z)*x))  +  cos(theta_x)*(cos(theta_z)*y - sin(theta_z)*x)
		d_z = cos(theta_x)*(cos(theta_y)*z + sin(theta_y)*(sin(theta_z)*y + cos(theta_z)*x))  -  sin(theta_x)*(cos(theta_z)*y - sin(theta_z)*x)
		# determine the positions on screen
		b_x[index] = (e_z/d_z) * d_x - e_x
		b_y[index] = (e_z/d_z) * d_y - e_y
	end


	min_i, max_i 	= -1,1
	min_j, max_j 	= -1,1
	border 			= border_percentage * (maximum(b_x) + maximum(b_y) - minimum(b_x) - minimum(b_y))/4


	# die farbe der sites
	function getProperSiteColor(distance, site_index)
		# get the color of that particular site
		site_color_bare = get(colorcode_sites, site_index, colorcode_sites[0])
		site_color_back = (site_color_bare .* (1-0.3*DOF)) .+ (site_color_BG .* 0.3*DOF)
		# the desired color of the site
		site_color_desired = [255, 255, 255]
		for i in 1:3
			# get the color max and min values
			site_color_max = site_color_back[i]
			site_color_min = site_color_bare[i]
			# parameters of the fit
			a = -(site_color_max-site_color_min)/((1/min_distance) - (1/max_distance))
			b = site_color_max - a * 1/max_distance
			# get the best grey value
			grey = (a / distance) + b
			grey = Int(floor(grey))
			# set the value
			site_color_desired[i] = grey
		end
		# return the hex version of the proper site color
		return site_color_desired
	end

	# proper radius
	radius_min = site_radius
	radius_max = Int(floor(radius_min/2.0))
	function getProperSiteRadius(distance)
		# overwrite the radius
		radius_max = radius_min*min_distance/max_distance
		# parameters of the fit
		a = -(radius_max-radius_min)/((1/min_distance) - (1/max_distance))
		b = radius_max - a * 1/max_distance
		# get the best radius
		radius = Int(floor((a / distance) + b))
		return radius
	end
	function getProperSiteBorderRadius(distance)
		# overwrite the radius
		radius_max = radius_min*min_distance/max_distance
		# parameters of the fit
		a = -(radius_max-radius_min)/((1/min_distance) - (1/max_distance))
		b = radius_max - a * 1/max_distance
		# get the best radius
		radius = ((a / distance) + b)*site_border_width_percentage
		return "$(radius)px"
	end

	# proper connection color
	function getProperConnectionColor(distance, strength)
		# parse to string
		strength = string(strength)
		# get the color of that particular bond
		bond_color_bare = get(colorcode_bonds, strength, colorcode_bonds["0"])
		bond_color_back = (bond_color_bare .* (1-DOF)) .+ (bond_color_BG .* DOF)
		# the desired color of the site
		bond_color_desired = [255, 255, 255]
		for i in 1:3
			# get the color max and min values
			bond_color_max = bond_color_back[i]
			bond_color_min = bond_color_bare[i]
			# parameters of the fit
			a = -(bond_color_max-bond_color_min)/((1/min_distance) - (1/max_distance))
			b = bond_color_max - a * 1/max_distance
			# get the best grey value
			grey = (a / distance) + b
			grey = Int(floor(grey))
			# set the value
			bond_color_desired[i] = grey
		end
		# return the hex version of the proper site color
		return color_hex(bond_color_desired)
	end
	function getProperSiteBorderColor(distance)
		# get the color of that 0 bond
		bond_color_bare = get(colorcode_bonds, "0", colorcode_bonds["0"])
		bond_color_back = (bond_color_bare .* (1-0.1*DOF)) .+ (bond_color_BG .* 0.1*DOF)
		# the desired color of the site
		bond_color_desired = [255, 255, 255]
		for i in 1:3
			# get the color max and min values
			bond_color_max = bond_color_back[i]
			bond_color_min = bond_color_bare[i]
			# parameters of the fit
			a = -(bond_color_max-bond_color_min)/((1/min_distance) - (1/max_distance))
			b = bond_color_max - a * 1/max_distance
			# get the best grey value
			grey = (a / distance) + b
			grey = Int(floor(grey))
			# set the value
			bond_color_desired[i] = grey
		end
		# return the hex version of the proper site color
		return color_hex(bond_color_desired)
	end

	# proper width
	width_min = bond_thickness
	width_max = width_min * 0.6
	function getProperConnectionWidth(distance)
		# parameters of the fit
		a = -(width_max-width_min)/((1/min_distance) - (1/max_distance))
		b = width_max - a * 1/max_distance
		# get the best radius
		width = (a / distance) + b
		return width
	end


	# DEFINE SITES THAT ARE PLOTTED

	# sites to plot
	sites_to_plot = Array[]
	for index in 1:length(b_x)
        if site_labels == "LATTICE INDEX"
		    push!(sites_to_plot, [b_x[index], b_y[index], distances_to_camera[index], sum(positions[index].*c_vec), lattice.positions_indices[index], index])
        else
            push!(sites_to_plot, [b_x[index], b_y[index], distances_to_camera[index], sum(positions[index].*c_vec), lattice.positions_indices[index]])
        end
	end

	# sort the sites s.t. the z values arrange in descending order
	function sortfunctionSites(site)
		return -site[4]
	end
	sort!(sites_to_plot, by=sortfunctionSites)


	# DEFINE CONNECTIONS THAT ARE PLOTTED

	# connections to plot (sort out periodic ones)
	connections_to_plot = Array[]
	neutral_connection_wrap = (0,0,0)
	if length(lattice.connections[1][4]) == 2
		neutral_connection_wrap = (0,0)
	elseif length(lattice.connections[1][4]) == 1
		neutral_connection_wrap = 0
	end
	for c in connections
		if c[4] != neutral_connection_wrap && !visualize_periodic
			continue
		end
		# startpositions
		p1 = positions[Int(c[1])]
		p2 = positions[Int(c[2])]
		pm = (p1.+p2) * 0.5
		# calculate the camera distance
		pm_distance = vecabs(pm.-c_vec)
		# push the connection
		#push!(connections_to_plot, [c; pm_distance; pm[3]-c_z])
		push!(connections_to_plot, [c; pm_distance; min(sum(p1.*c_vec), sum(p2.*c_vec)) - 0.01])
	end
	# sort the connections s.t. the z values arrange in descending order
	function sortfunctionConnections(con)
		return -con[6] #max(con[5],con[6])
	end
	sort!(connections_to_plot, by=sortfunctionConnections)





	# overwrite positions from here on with projected positions
	for index in 1:length(b_x)
		positions[index] = [b_x[index],b_y[index]]
	end

	# generate the xvals and yvals only for size
	xvals = b_x
	yvals = b_y
	# define the width and height
	width_UC 	= (maximum(xvals) + border) - (minimum(xvals) - border)
	height_UC	= (maximum(yvals) + border) - (minimum(yvals) - border)
	min_x	= minimum(xvals) - border
	min_y	= minimum(yvals) - border
	if width_UC > height_UC
		width	= 1200
		conversion = width/width_UC
		height	= conversion*(height_UC)
	else
		height	= 1200
		conversion = height/height_UC
		width	= conversion*(width_UC)
	end

	# max und min distance
	max_distance = maximum(distances_to_camera)
	min_distance = minimum(distances_to_camera)

	# define the conversion functions for coordinates
	function X(x)
		return conversion*(x - min_x)
	end
	function Y(y)
		return + conversion*(y - min_y)
	end

	# open the file
	file = open(filename_output, "w")

	# write the headerstring
	write(file, getSVGHeaderString(round(Int64, width), round(Int64, height)))


	# write all strings
	i = 0
	while size(connections_to_plot,1) > 0 || size(sites_to_plot,1) > 0
		# if no more connections, only plot sites
		if size(connections_to_plot,1) == 0
			# pop the first site
			site_to_plot = pop!(sites_to_plot)
            site_color_basic = getProperSiteColor(site_to_plot[3],site_to_plot[5])
            site_color = color_hex(site_color_basic)
            label_color = color_hex([colorelment < 100 ? 255 : 0 for colorelment in site_color_basic])
            # write
            if site_labels == "POSITION INDEX" || site_labels == "LATTICE INDEX"
                write(file, getSVGStringEllipseStroked("el$(i)",
                    X(site_to_plot[1]), Y(site_to_plot[2]),
                    getProperSiteRadius(site_to_plot[3]), getProperSiteRadius(site_to_plot[3]),
                    site_color,
                    getProperSiteBorderColor(site_to_plot[3]), getProperSiteBorderRadius(site_to_plot[3]),
                    label=Int(site_to_plot[end]), labelcolor=label_color))
            else
                write(file, getSVGStringEllipseStroked("el$(i)",
                    X(site_to_plot[1]), Y(site_to_plot[2]),
                    getProperSiteRadius(site_to_plot[3]), getProperSiteRadius(site_to_plot[3]),
                    site_color,
                    getProperSiteBorderColor(site_to_plot[3]), getProperSiteBorderRadius(site_to_plot[3])))
            end
			# increase index
			i = i+1
		elseif size(sites_to_plot,1) == 0
			# pop the first connection
			connection_to_plot = pop!(connections_to_plot)
			c = connection_to_plot
			# write
			if visualize_periodic || c[4] == neutral_connection_wrap
                write(file,
                    getSVGStringLine("path$(i)",
                    [X(positions[Int(c[1])][1]), Y(positions[Int(c[1])][2])],
                    [X(positions[Int(c[2])][1]), Y(positions[Int(c[2])][2])],
                    getProperConnectionColor(c[5], c[3]),
                    getProperConnectionWidth(c[5]), dashed=(c[4] != neutral_connection_wrap)))
            end
			# increase index
			i = i+1
		else
			# decide which one to pop
			if sites_to_plot[end][4] < connections_to_plot[end][6] #min(connections_to_plot[end][5],connections_to_plot[end][6]) - 0.0001
				# pop the first site
				site_to_plot = pop!(sites_to_plot)
                site_color_basic = getProperSiteColor(site_to_plot[3],site_to_plot[5])
                site_color = color_hex(site_color_basic)
                label_color = color_hex([colorelment < 100 ? 255 : 0 for colorelment in site_color_basic])
                # write
                if site_labels == "POSITION INDEX" || site_labels == "LATTICE INDEX"
                    write(file, getSVGStringEllipseStroked("el$(i)",
                        X(site_to_plot[1]), Y(site_to_plot[2]),
                        getProperSiteRadius(site_to_plot[3]), getProperSiteRadius(site_to_plot[3]),
                        site_color,
                        getProperSiteBorderColor(site_to_plot[3]), getProperSiteBorderRadius(site_to_plot[3]),
                        label=Int(site_to_plot[end]), labelcolor=label_color))
                else
                    write(file, getSVGStringEllipseStroked("el$(i)",
                        X(site_to_plot[1]), Y(site_to_plot[2]),
                        getProperSiteRadius(site_to_plot[3]), getProperSiteRadius(site_to_plot[3]),
                        site_color,
                        getProperSiteBorderColor(site_to_plot[3]), getProperSiteBorderRadius(site_to_plot[3])))
                end
				# increase index
				i = i+1
			else
				# pop the first connection
				connection_to_plot = pop!(connections_to_plot)
				c = connection_to_plot
				# write
				if visualize_periodic || c[4] == neutral_connection_wrap
					write(file,
						getSVGStringLine("path$(i)",
						[X(positions[Int(c[1])][1]), Y(positions[Int(c[1])][2])],
						[X(positions[Int(c[2])][1]), Y(positions[Int(c[2])][2])],
						getProperConnectionColor(c[5], c[3]),
						getProperConnectionWidth(c[5]), dashed=(c[4] != neutral_connection_wrap)))
				end
			end
			# increase index
			i = i+1
		end
	end






	# write the footerstring to close the svg file
	write(file, getSVGFooterString())
	# close the file
	close(file)

	# convert to pdf
    if export_pdf
	    run(`inkscape $(filename_output) --export-pdf $(filename_output[1:end-4]).pdf`)
    end

	# if file shall be opened
	if openfile
        if is_linux()
		    run(`ristretto $(filename_output)`)
        elseif is_windows()
            try
                if export_pdf
                    run(`explorer $(filename_output[1:end-4]).pdf`)
                else
                    run(`explorer $(filename_output)`)
                end
            end
        elseif is_apple()
            run(`open $(filename_output)`)
        else
            println("wanted to show file but operating system not supported for showing file")
        end
	end

	# return the output filename
	return filename_output

end
export plotLattice3D


function plotLattice(
		lattice::Lattice;
		conversion = 160,
		border_percentage=0.1,
		filename_output::String="AUTO",
		site_radius=25,
		site_border_width_percentage::Float64=0.2,
        site_labels="OFF",
		bond_thickness::Int64=8,
		visualize_periodic=false,
		colorcode_sites = Dict(0 => [255,255,255], 1 => [255,255,255]),
		colorcode_bonds = Dict("0" => [0,0,0], "1.0" => [0,0,0]),
        colorcode_bonds_automation::String = "OFF",
		bond_color_BG = [100,100,100],
		site_color_BG = [125,125,125],
		lattice_rotation::Array{Float64}=[0.0,0.0,0.0],
		camera_position_relative::Array{Float64}=[0.0,0.0,5.0],
		DOF::Float64 = 0.8,
		openfile=false,
        export_pdf=true
		)

    # check which dimension the lattice has
    dimension = length(lattice.positions[1])

    # determine which function to call, depending on dimension
    if dimension == 2
        return plotLattice2D(
            lattice,
		    conversion=conversion,
            border_percentage=border_percentage,
            filename_output=filename_output,
            site_radius=site_radius,
            site_border_width_percentage=site_border_width_percentage,
            site_labels = site_labels,
            bond_thickness=bond_thickness,
            visualize_periodic=visualize_periodic,
            colorcode_sites = colorcode_sites,
            colorcode_bonds = colorcode_bonds,
            colorcode_bonds_automation = colorcode_bonds_automation,
            openfile=openfile,
            export_pdf=export_pdf)
    elseif dimension == 3
        return plotLattice3D(
            lattice,
            border_percentage=border_percentage,
            filename_output=filename_output,
            site_radius=site_radius,
            site_border_width_percentage=site_border_width_percentage,
            site_labels = site_labels,
            bond_thickness=bond_thickness,
            visualize_periodic=visualize_periodic,
            colorcode_sites = colorcode_sites,
            colorcode_bonds = colorcode_bonds,
            colorcode_bonds_automation = colorcode_bonds_automation,
            bond_color_BG = bond_color_BG,
            site_color_BG = site_color_BG,
            lattice_rotation=lattice_rotation,
            camera_position_relative=camera_position_relative,
            DOF=DOF,
            openfile=openfile,
            export_pdf=export_pdf)
    else
        println("plotting to SVG not implemented for dimension $(dimension) of the lattice")
        return
    end
end
export plotLattice






function plotPlaquettes2D(
		lattice::Lattice,
        plaquettes,
        plaquette_values;
		conversion = 160,
		border_percentage=0.1,
		filename_output::String="AUTO",
		site_radius=25,
		site_border_width_percentage::Float64=0.2,
        site_labels="OFF",
		bond_thickness::Int64=8,
		visualize_periodic=false,
		colorcode_sites = Dict(0 => [255,255,255], 1 => [255,255,255]),
		colorcode_bonds = Dict("0" => [0,0,0], "1.0" => [0,0,0]),
		colorcode_plaquettes = Dict(0.0 => [0,0,0], 1.0 => [255,0,0], -1.0 => [0,100,255]),
        colorcode_bonds_automation::String = "OFF",
		openfile=false,
        export_pdf=true
		)

    # define the filename_output if it is set to AUTO
    if filename_output=="AUTO"
        filename_output = "$(lattice.filename[1:end-4])_plq_plot.svg"
        filename_output = replace(filename_output, FOLDER_LATTICES, "")
    end

	# load positions and connections
	positions	= lattice.positions
	connections = lattice.connections

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
    colorcode_plaquettes[0.0] = get(colorcode_plaquettes, 0.0, [0,0,0])

	# define styles for the different sites, i.e. strings that are saved into the svg strings
	site_r 				= site_radius
	site_border			= "#000000"
	site_border_width	= "$(site_border_width_percentage*site_radius)px"

	# sites to plot
	pos_x = zeros(size(positions,1))
	pos_y = zeros(size(positions,1))
	for i in 1:size(positions,1)
		pos_x[i] = positions[i][1]
		pos_y[i] = positions[i][2]
	end
	sites_to_plot = positions
	indices_to_plot = lattice.positions_indices
	xvals = pos_x
	yvals = pos_y
	border 			= border_percentage * (maximum(pos_x) + maximum(pos_y) - minimum(pos_x) - minimum(pos_y))/4

	# connections to plot (all)
	connections_to_plot = copy(connections)
    # the neutral connection wrap, i.e. which wrap identifies a non-periodic connection
	neutral_connection_wrap = (0,0)
	if size(lattice.lattice_vectors,1) == 1
		neutral_connection_wrap = (0)
	end


	# define the width and height of the canvas
	width_UC 	= (maximum(xvals) + site_radius/conversion + border) - (minimum(xvals) - site_radius/conversion - border)
	height_UC	= (maximum(yvals) + site_radius/conversion + border) - (minimum(yvals) - site_radius/conversion - border)
	min_x	= minimum(xvals) - site_radius/conversion - border
	min_y	= minimum(yvals) - site_radius/conversion - border
	width	= conversion*(width_UC)
	height	= conversion*(height_UC)

	# define the conversion functions for coordinates
	function X(x)
		return conversion*(x - min_x)
	end
	function Y(y)
		return + conversion*(y - min_y)
	end

	# open the SVG file
	file = open(filename_output, "w")

	# write the headerstring
	write(file, getSVGHeaderString(round(Int64, width), round(Int64, height)))

    # write all plaquettes
    for (i,p) in enumerate(plaquettes)
        # get the color
        color_plaquette = color_hex(get(colorcode_plaquettes, plaquette_values[i], colorcode_plaquettes[0.0]))
        # get the plaquette point list
        points = Array[]
        for index in p
            push!(points, [X(positions[index][1]), Y(positions[index][2])])
        end
        # plot plaquette
        write(file, getSVGStringPlaquette("plaq$(i)", points, color_plaquette))
    end
	# write all connections
	for (i,c) in enumerate(connections_to_plot)
		if (c[4] == neutral_connection_wrap)
			connection_color = color_hex(get(colorcode_bonds, string(c[3]), colorcode_bonds["0"]))
			write(file, getSVGStringLine("path$(i)", [X(positions[Int(c[1])][1]), Y(positions[Int(c[1])][2])], [X(positions[Int(c[2])][1]), Y(positions[Int(c[2])][2])], connection_color, bond_thickness))
		elseif visualize_periodic
			connection_color = color_hex(get(colorcode_bonds, string(c[3]), colorcode_bonds["0"]))
			write(file, getSVGStringLine("path$(i)", [X(positions[Int(c[1])][1]), Y(positions[Int(c[1])][2])], [X(positions[Int(c[2])][1]), Y(positions[Int(c[2])][2])], connection_color, bond_thickness, dashed=true))
        end
	end

	# write all sites
	for (i,s) in enumerate(sites_to_plot)
        site_color_basic = get(colorcode_sites, indices_to_plot[i], colorcode_sites[0])
		site_color = color_hex(site_color_basic)
        label_color = color_hex([colorelment < 100 ? 255 : 0 for colorelment in site_color_basic])
        if site_labels == "POSITION INDEX"
		    write(file, getSVGStringEllipseStroked("el$(i)", X(s[1]), Y(s[2]), site_r, site_r, site_color, site_border, site_border_width, label=indices_to_plot[i], labelcolor=label_color))
        elseif site_labels == "LATTICE INDEX"
		    write(file, getSVGStringEllipseStroked("el$(i)", X(s[1]), Y(s[2]), site_r, site_r, site_color, site_border, site_border_width, label=i, labelcolor=label_color))
        else
		    write(file, getSVGStringEllipseStroked("el$(i)", X(s[1]), Y(s[2]), site_r, site_r, site_color, site_border, site_border_width))
        end
	end

	# write the footerstring to close the svg file
	write(file, getSVGFooterString())
	# close the file
	close(file)

	# convert to pdf
    if export_pdf
	    run(`inkscape $(filename_output) --export-pdf $(filename_output[1:end-4]).pdf`)
    end
	# if file shall be opened
	if openfile
        if is_linux()
		    run(`ristretto $(filename_output)`)
        elseif is_windows()
            try
                if export_pdf
                    run(`explorer $(filename_output[1:end-4]).pdf`)
                else
                    run(`explorer $(filename_output)`)
                end
            end
        elseif is_apple()
            run(`open $(filename_output)`)
        else
            println("wanted to show file but operating system not supported for showing file")
        end
	end

	# return the output filename
	return filename_output

end
export plotPlaquettes2D




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
#-----------------------------------------------------------------------------------------------------------------------------
#
#   CALCULATIONS OF BAND STRUCTURE AND SPIN GROUNDSTATES
#
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------






#-----------------------------------------------------------------------------------------------------------------------------
#
#   METHODS FOR CONSTRUCTION INTERACTION MATRICES FOR LATTICES
#
#-----------------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------------
#
#   Interaction matrix in REAL space
#
#   Parameters are
#   - lattice: The complete lattice of which the interaction matrix should be constructed
#   - enforce_hermitian (optional): if the matrix should be made hermitian by 0.5*(A + A_dagger)
#
#-----------------------------------------------------------------------------------------------------------------------------
function getInteractionMatrixRealSpace(lattice::Lattice; enforce_hermitian=false)
    # generate a new matrix
    matrix = zeros(size(lattice.positions,1),size(lattice.positions,1))
    # iterate over all connections
    for c in lattice.connections
        # get the indices
        index_from  = Int(c[1])
        index_to    = Int(c[2])
        strength    = c[3]
        # just add to the matrix
        matrix[index_from, index_to] += strength
    end
    # eventually ensure the hermitian nature of the matrix
    if enforce_hermitian
        matrix = 0.5*(matrix .+ transpose(conj(matrix)))
    end
    # return the matrix
    return matrix
end
function getInteractionMatrixRealSpace(unitcell::Unitcell; enforce_hermitian=false)
    # generate a new matrix
    matrix = zeros(size(unitcell.basis,1),size(unitcell.basis,1))
    # iterate over all connections
    for c in unitcell.connections
        # get the indices
        index_from  = Int(c[1])
        index_to    = Int(c[2])
        strength    = c[3]
        # just add to the matrix
        matrix[index_from, index_to] += strength
    end
    # eventually ensure the hermitian nature of the matrix
    if enforce_hermitian
        matrix = 0.5*(matrix .+ transpose(conj(matrix)))
    end
    # return the matrix
    return matrix
end
export getInteractionMatrixRealSpace

#-----------------------------------------------------------------------------------------------------------------------------
#
#   Interaction matrix in RECIPROCAL (K) space
#
#   Parameters are
#   - lattice: The complete lattice of which the interaction matrix should be constructed
#   - k_vector: The reciprocal point k at which the matrix is constructed
#   - enforce_hermitian (optional): if the matrix should be made hermitian by 0.5*(A + A_dagger)
#
#-----------------------------------------------------------------------------------------------------------------------------
function getInteractionMatrixKSpace(lattice::Lattice, k_vector::Array{Float64,1}; enforce_hermitian=false, majorana=false)
    # generate a new matrix
    matrix = zeros(Complex, size(lattice.positions,1),size(lattice.positions,1))
    # iterate over all connections
    for c in lattice.connections
        # get the indices
        index_from  = Int(c[1])
        index_to    = Int(c[2])
        strength    = c[3]
        wrap        = c[4]
        # get the difference vector
        pos_delta   = lattice.positions[index_to] .- lattice.positions[index_from]
        if size(lattice.lattice_vectors,1) > 0
            for pair in zip(wrap, lattice.lattice_vectors)
                pos_delta .+= pair[1].*pair[2]
            end
        end
        # take majorana fermionic statistics into account
        majorana_factor = 1
        if majorana && index_from < index_to
            majorana_factor = 1.0im
        elseif majorana && index_from > index_to
            majorana_factor = -1.0im
        end
        # just add to the matrix
        matrix[index_from, index_to] += strength * exp(-sum(pos_delta.*k_vector) * im) * majorana_factor
    end
    # eventually ensure the hermitian nature of the matrix
    if enforce_hermitian
        matrix = 0.5*(matrix .+ transpose(conj(matrix)))
    end
    # return the matrix
    return matrix
end
function getInteractionMatrixKSpace(unitcell::Unitcell, k_vector::Array{Float64,1}; enforce_hermitian=false, majorana=false)
    # generate a new matrix
    matrix = zeros(Complex, size(unitcell.basis,1),size(unitcell.basis,1))
    # iterate over all connections
    for c in unitcell.connections
        # get the indices
        index_from  = Int(c[1])
        index_to    = Int(c[2])
        strength    = c[3]
        wrap        = c[4]
        # get the difference vector
        pos_delta   = unitcell.basis[index_to] .- unitcell.basis[index_from]
        if size(unitcell.lattice_vectors,1) > 0
            for pair in zip(wrap, unitcell.lattice_vectors)
                pos_delta .+= pair[1].*pair[2]
            end
        end
        # take majorana fermionic statistics into account
        majorana_factor = 1
        if majorana && index_from < index_to
            majorana_factor = 1.0im
        elseif majorana && index_from > index_to
            majorana_factor = -1.0im
        end
        # just add to the matrix
        matrix[index_from, index_to] += strength * exp(-sum(pos_delta.*k_vector) * im) * majorana_factor
    end
    # eventually ensure the hermitian nature of the matrix
    if enforce_hermitian
        matrix = 0.5*(matrix .+ transpose(conj(matrix)))
    end
    # return the matrix
    return matrix
end
export getInteractionMatrixKSpace

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
