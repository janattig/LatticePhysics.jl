################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT LATTICE PLOTTING FUNCTIONS FOR
#   PLOTTING TO SVG IMAGES
#
#   STRUCTURE OF THE FILE
#
#   1) HELPER SVG FUNCTIONS
#
#   2) PLOTTING LATTICES
#
#   3) PLOTTING PLAQUETTES
#
################################################################################





################################################################################
#
#   HELPER METHODS FOR SVG CREATION
#   Get the correct SVG code lines for SVG beginning and end of document
#   Get the correct SVG code lines for certain objects
#
#   NOTE:   ALL OF THESE FUNCTIONS ARE NOT EXPORTED BUT ONLY USED!
#           ( THEREFORE NO DOCSTRINGS NECESSARY )
#
################################################################################







# Function to construct a HEADER STRING (must be first in every SVG file)
# width and height denote the dimensions of the image in px
function getSVGHeaderString(width::Int64, height::Int64)
	headerstring = """<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n<svg
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
	</metadata>\n\n"""
	return headerstring
end

# Function to construct the FOOTER STRING (must be end of every SVG file)
function getSVGFooterString()
    return """\n</svg>\n"""
end




# STRINGS FOR AN ELLIPSE
# Parameters are
# - id: The id which the object has later on
# - centerX / centerY: The coordinates of the ellipse center in the canvas
# - radiusX / radiusY: The radii of the ellipse
# - color: hex string of the ellipse color
# Optional:
# - label: The label of the ellipse
# - labelcolor: the color of the label
# - opacity: the opacity of filling the ellipse
function getSVGStringEllipse(
            id,
            centerX::Float64, centerY::Float64,
            radiusX::Float64, radiusY::Float64,
            color::String;
            label::String="NONE",
            labelcolor::String="#000000",
            opacity::Float64=1.0
        )
    # construct the String with only the ellipse
	es = """
	<ellipse
		style=\"color:$(color);fill:$(color);fill-opacity:$(opacity);fill-rule:nonzero\"
		id=\"$(id)\"
		cx=\"$(centerX)\"
		cy=\"$(centerY)\"
		rx=\"$(radiusX)px\"
		ry=\"$(radiusY)px\" />\n"""
    # check if to add a label
    if label != "NONE"
        es = """
    $(es)

    <text
        x=\"$(centerX)\"
        y=\"$(centerY+radiusY/2)\"
        style=\"text-anchor: middle; font-size: $(radiusX*1.2)px; fill:$(labelcolor)\"
    >
        $(label)
    </text>\n"""
    end
    # return the string
	return es
end

# OTHER VERSIONS OF ELLIPSE STRING
function getSVGStringEllipse(
            id,
            centerX::Float64, centerY::Float64,
            radius::Float64,
            color::String;
            label::String="NONE",
            labelcolor::String="#000000",
            opacity::Float64=1.0
        )
    return getSVGStringEllipse(id, centerX, centerY, radius, radius, color, label=label, labelcolor=labelcolor, opacity=opacity)
end


# STRING FOR A STROKED ELLIPSE
# Parameters are
# - id: The id which the object has later on
# - centerX / centerY: The coordinates of the ellipse center in the canvas
# - radiusX / radiusY: The radii of the ellipse
# - colorFill: hex string of the ellipse fill color
# - colorStroke: hex string of the ellipse stroke color
# - strokewidth: Float or Int giving the width of the surrounding stroke
# Optional:
# - label: The label of the ellipse
# - labelcolor: the color of the label
# - opacity: the opacity of filling the ellipse
function getSVGStringEllipseStroked(
            id,
            centerX::Float64, centerY::Float64,
            radiusX::Float64, radiusY::Float64,
            colorFill::String,
            colorStroke::String,
            strokewidth::Float64;
            label::String="NONE",
            labelcolor::String="#000000",
            opacityFill::Float64=1.0,
            opacityStroke::Float64=1.0
        )
    # construct the string for the ellipse
	es = """
	<ellipse
		style=\"color:$(colorFill);fill:$(colorFill);fill-opacity:$(opacityFill);fill-rule:nonzero;stroke:$(colorStroke);stroke-width:$(strokewidth);stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-dasharray:none;stroke-dashoffset:0;stroke-opacity:$(opacityStroke)\"
		id=\"$(id)\"
		cx=\"$(centerX)\"
		cy=\"$(centerY)\"
		rx=\"$(radiusX)px\"
		ry=\"$(radiusY)px\" />\n"""
    # check if to add a label
    if label != "NONE"
        es = """
    $(es)

    <text
        x=\"$(centerX)\"
        y=\"$(centerY+radiusY/2)\"
        style=\"text-anchor: middle; font-size: $(radiusX*1.2)px; fill:$(labelcolor)\"
    >
        $(label)
    </text>\n"""
    end
    # return the string
	return es
end


# OTHER VERSIONS OF ELLIPSE STRING
function getSVGStringEllipseStroked(
            id,
            centerX::Float64, centerY::Float64,
            radius::Float64,
            colorFill::String,
            colorStroke::String,
            strokewidth::Float64;
            label::String="NONE",
            labelcolor::String="#000000",
            opacityFill::Float64=1.0,
            opacityStroke::Float64=1.0
        )
    return getSVGStringEllipseStroked(id, centerX, centerY, radius, radius, colorFill, colorStroke, strokewidth, label=label, labelcolor=labelcolor, opacityFill=opacityFill, opacityStroke=opacityStroke)
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
#-----------------------------------------------------------------------------------------------------------------------------
#
#   PLOTTING TO SVG IMAGES
#
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------









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
