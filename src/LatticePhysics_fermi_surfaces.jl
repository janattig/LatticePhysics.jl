################################################################################
#
#   METHODS FOR CONSTRUCTION OF FERMI SURFACES (2D & 3D)
#
#   STRUCTURE OF THE FILE
#
#   1) CALCULATION OF FERMI SURFACE
#
#   2) PLOTTING OF FERMI SURFACE
#       - plotting from points
#       - plotting from unitcell
#
################################################################################




# Function for calculating Fermi surface in 2D (not exported)
function getFermiSurface2D(
            unitcell::Unitcell,
            N_points::Int64;
            fermi_energy::Float64=0.0,
            enforce_hermitian::Bool=false,
            epsilon::Float64=1e-10,
            epsilon_k::Float64=1e-10,
            bounds_lower::Array{Float64,1}=-2*pi.*ones(2),
            bounds_upper::Array{Float64,1}=2*pi.*ones(2),
            refold_to_first_BZ::Bool=true
        )

    # list of k values that contribute to the BZ
    k_values = zeros(Float64, N_points, 2)

    # funtion of energy
    function energy(kvector)
        eigenvalues = eigvals(getInteractionMatrixKSpace(unitcell, kvector, enforce_hermitian=enforce_hermitian)) .- fermi_energy
        eigenvalues = eigenvalues .* eigenvalues
        return minimum(eigenvalues)
    end

    # the current search index
    index = 1
    # search until there are enough points
    while index <= N_points

        # find a suitable starting point for the Newton algorithm
        k = Float64[rand(), rand()]
        for j in 1:length(k)
            k[j] = k[j] * bounds_lower[j] + (1-k[j]) * bounds_upper[j]
        end

        # start with the initial energy
        e0 = energy(k)
        # iterate i over 100 newton steps (maximum)
        for i in 1:100
            # check if the energy is already converged
            if e0 < epsilon
                # save the k vector
                k_values[index,:] = k
                # increment the index
                index = index+1
                # break the newton loop
                break
            end
            # the current energy
            H_0 = e0
            # the gradient of the energy
            H_eps =
            [
                energy(k .+ [epsilon_k, 0]),
                energy(k .+ [0, epsilon_k])
            ]
            dH = (H_eps - H_0) ./ epsilon_k
            # absolute value of the gradient
            dHdH = dot(dH, dH)
            # break the Newton loop if the gradient diverges or is flattened too much
            if abs(dHdH) < 1e-20 || abs(dHdH) > 1e20
                break
            end
            # increment k
            dk = dH * (H_0 / dHdH)
            k -= dk
            # calculate a new energy
            e0 = energy(k)
        end
    end

    # if refolding is enabled, refold all points
    if refold_to_first_BZ
        # get the lattice vectors
        a1 = unitcell.lattice_vectors[1]
        a2 = unitcell.lattice_vectors[2]
        # get the reciprocal lattice vectors
        b1 = [a2[2], -a2[1]]
        b2 = [a1[2], -a1[1]]
        # normalize the vectors
        b1 .*= 2*pi/sum(b1.*a1)
        b2 .*= 2*pi/sum(b2.*a2)

        # iterate over all k points
        for index in 1:N_points
            # assumend to be refolded to begin with
            refolded = true
            # iterate while still refolding possible
            while refolded
                # new iteration, not refolded
                refolded = false
                # try all possible refoldings
                for r in [b1, b2, -b1, -b2, b1+b2, -b2-b1, b1-b2, -b2+b1]
                    if sum((k_values[index,:].-r).*(k_values[index,:].-r)) < sum((k_values[index,:]).*(k_values[index,:]))
                        k_values[index,:] .-= r
                        refolded=true
                    end
                end
            end
        end
    end

    # return the list
    return k_values

end

# Function for calculating Fermi surface in 3D (not exported)
function getFermiSurface3D(
            unitcell::Unitcell,
            N_points::Int64;
            fermi_energy::Float64=0.0,
            enforce_hermitian::Bool=false,
            epsilon::Float64=1e-10,
            epsilon_k::Float64=1e-10,
            bounds_lower::Array{Float64,1}=-2*pi.*ones(3),
            bounds_upper::Array{Float64,1}=2*pi.*ones(3),
            refold_to_first_BZ::Bool=true
        )

    # list of k values that contribute to the BZ
    k_values = zeros(Float64, N_points, 3)

    # funtion of energy
    function energy(kvector)
        eigenvalues = eigvals(getInteractionMatrixKSpace(unitcell, kvector, enforce_hermitian=enforce_hermitian)) .- fermi_energy
        eigenvalues = eigenvalues .* eigenvalues
        return minimum(eigenvalues)
    end

    # the current search index
    index = 1
    # search until there are enough points
    while index <= N_points

        # find a suitable starting point for the Newton algorithm
        k = Float64[rand(), rand(), rand()]
        for j in 1:length(k)
            k[j] = k[j] * bounds_lower[j] + (1-k[j]) * bounds_upper[j]
        end

        # start with the initial energy
        e0 = energy(k)
        # iterate i over 100 newton steps (maximum)
        for i in 1:100
            # check if the energy is already converged
            if e0 < epsilon
                # save the k vector
                k_values[index,:] = k
                # increment the index
                index = index+1
                # break the newton loop
                break
            end
            # the current energy
            H_0 = e0
            # the gradient of the energy
            H_eps =
            [
                energy(k .+ [epsilon_k, 0, 0]),
                energy(k .+ [0, epsilon_k, 0]),
                energy(k .+ [0, 0, epsilon_k])
            ]
            dH = (H_eps - H_0) ./ epsilon_k
            # absolute value of the gradient
            dHdH = dot(dH, dH)
            # break the Newton loop if the gradient diverges or is flattened too much
            if abs(dHdH) < 1e-20 || abs(dHdH) > 1e20
                break
            end
            # increment k
            dk = dH * (H_0 / dHdH)
            k -= dk
            # calculate a new energy
            e0 = energy(k)
        end
    end

    # if refolding is enabled, refold all points
    if refold_to_first_BZ
        # get the lattice vectors
        a1 = unitcell.lattice_vectors[1]
        a2 = unitcell.lattice_vectors[2]
        a3 = unitcell.lattice_vectors[3]
        # get the reciprocal lattice vectors
        b1 = cross(a2, a3)
        b2 = cross(a3, a1)
        b3 = cross(a1, a2)
        # normalize the vectors
        b1 .*= 2*pi/sum(b1.*a1)
        b2 .*= 2*pi/sum(b2.*a2)
        b3 .*= 2*pi/sum(b3.*a3)
        # put together refolding
        refoldings = Array{Float64,1}[]
        for i in [-1,0,1]
        for j in [-1,0,1]
        for l in [-1,0,1]
            refolding = i*b1 + j*b2 + l*b3
            if sum(abs.(refolding))>1e-8
                push!(refoldings, refolding)
            end
        end
        end
        end
        # iterate over all k points
        for index in 1:N_points
            # assumend to be refolded to begin with
            refolded = true
            # iterate while still refolding possible
            while refolded
                # new iteration, not refolded
                refolded = false
                # try all possible refoldings
                for r in refoldings
                    if sum((k_values[index,:].-r).*(k_values[index,:].-r)) < sum((k_values[index,:]).*(k_values[index,:]))
                        k_values[index,:] .-= r
                        refolded=true
                    end
                end
            end
        end
    end

    # return the list
    return k_values

end





# obtain Fermi surface
"""
    getFermiSurface(
                unitcell::Unitcell,
                N_points::Int64
             [; fermi_energy::Float64=0.0,
                enforce_hermitian::Bool=false,
                epsilon::Float64=1e-10,
                epsilon_k::Float64=1e-10,
                bounds_lower::Array{Float64,1}=-2*pi.*ones(4),
                bounds_upper::Array{Float64,1}=2*pi.*ones(4),
                refold_to_first_BZ::Bool=true ]
            )

calculates `N_points` points which belong to the Fermi surface of the tight-binding model given by
the `Unitcell` object `unitcell`. The Fermi energy can be adjusted by using `fermi_energy=...`.

The points are located using Newton's method with random starting positions. The procedure can be
modified by using some of the optional keywords:
- `enforce_hermitian` determines if the matrix used for energy calculation is made hermitian by construction
- `epsilon` energy threshold for values to be considered close to the Fermi energy
- `epsilon_k` small k increment that is used for calculation of derivatives
- `bounds_upper` / `bounds_lower` the bounds between which the random starting location for the Newton method is searched

After the search has been completed, the values will optionally be refolded into the first brillouin zone by choosing `refold_to_first_BZ=true`.




# Examples

```julia-repl
julia> fermi_surace = getFermiSurface(unitcell, 100)    # 2D unitcell
100x2 Array{Float64}
...

julia> fermi_surace = getFermiSurface(unitcell, 100)    # 3D unitcell
100x3 Array{Float64}
...

julia> fermi_surace = getFermiSurface(unitcell, 100, fermi_energy=1.0)
100x2 Array{Float64}
...

julia> fermi_surace = getFermiSurface(unitcell, 100, refold_to_first_BZ=false)
100x2 Array{Float64}
...

```
"""
function getFermiSurface(
            unitcell::Unitcell,
            N_points::Int64;
            fermi_energy::Float64=0.0,
            enforce_hermitian::Bool=false,
            epsilon::Float64=1e-10,
            epsilon_k::Float64=1e-10,
            bounds_lower::Array{Float64,1}=-2*pi.*ones(3),
            bounds_upper::Array{Float64,1}=2*pi.*ones(3),
            refold_to_first_BZ::Bool=true
        )
    # check if the unitcell can be put in either method
    if length(unitcell.basis[1]) != length(unitcell.lattice_vectors)
        println("Unitcell has not the same number of lattice vectors as dimensions")
        return zeros(1,1)
    end
    # check which function to pass to
    if length(unitcell.lattice_vectors) == 2
        # 2D case
        return getFermiSurface2D(
                unitcell,
                N_points,
                fermi_energy=fermi_energy,
                enforce_hermitian=enforce_hermitian,
                epsilon=epsilon,
                epsilon_k=epsilon_k,
                bounds_lower=bounds_lower,
                bounds_upper=bounds_upper,
                refold_to_first_BZ=refold_to_first_BZ
            )
    else length(unitcell.lattice_vectors) == 3
        # 3D case
        return getFermiSurface3D(
                unitcell,
                N_points,
                fermi_energy=fermi_energy,
                enforce_hermitian=enforce_hermitian,
                epsilon=epsilon,
                epsilon_k=epsilon_k,
                bounds_lower=bounds_lower,
                bounds_upper=bounds_upper,
                refold_to_first_BZ=refold_to_first_BZ
            )
    end

end
export getFermiSurface














# plot the Fermi surface (2D)
function plotFermiSurface2D(
            k_values::Array{Float64, 2};
            brillouin_zone::BrillouinZone=BrillouinZone(),
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,6),
            showPlot::Bool=true,
            save_filename::String="NONE"
        )


    ###########################
    #   INITIAL SETTINGS
    ###########################

    # configure plot environment
    rc("font", family="serif")

    # create a new figure
    fig = figure(figsize=figsize)



    ###########################
    #   PLOT FERMI SURFACE
    ###########################

    # if brillouin zone not empty, plot it as well
    if length(brillouin_zone.points) > 0
        plotBrillouinZone(brillouin_zone, new_figure=false)
    end

    # scatter the points
    scatter(k_values[:,1], k_values[:,2], color=plot_color)



    ###########################
    #   CONFIGURE AXIS & TITLE
    ###########################

    # set the title
    if plot_title == "AUTO"
        # set the title to an automatically generated title
        title("Fermi surface")
    elseif plot_title == ""
        # do nothing title related
    else
        # set the title to the given title
        title(plot_title)
    end

    # equal axis
    gca()[:set_aspect]("equal")



    ###########################
    #   FINISH THE PLOT
    ###########################

    # tighten the layout
    tight_layout()

    # save the plot
    if save_filename != "NONE"
        # make sure the directory exists
        if contains(save_filename, "/")
    		# get the containing folder
    		folder = save_filename[1:findlast(save_filename, '/')]
    		# build the path to that folder
    		mkpath(folder)
    	end
        # save the plot
        savefig(save_filename)
    end

    # maybe show the plot
    if showPlot
        show()
    end

    # return the figure object
    return fig

end
function plotFermiSurface2D(
            unitcell::Unitcell,
            N_points::Int64;
            fermi_energy::Float64=0.0,
            enforce_hermitian::Bool=false,
            epsilon::Float64=1e-10,
            epsilon_k::Float64=1e-10,
            bounds_lower::Array{Float64,1}=-2*pi.*ones(4),
            bounds_upper::Array{Float64,1}=2*pi.*ones(4),
            refold_to_first_BZ::Bool=true,
            brillouin_zone::BrillouinZone=BrillouinZone(),
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,6),
            showPlot::Bool=true,
            save_filename::String="NONE"
        )

    # check if BZ should be calculated
    if length(brillouin_zone.points) == 0
        brillouin_zone = createBrillouinZone(unitcell)
    end

    # calculate the Fermi surface first
    fermi_surface = getFermiSurface2D(
            unitcell,
            N_points,
            fermi_energy=fermi_energy,
            enforce_hermitian=enforce_hermitian,
            epsilon=epsilon,
            epsilon_k=epsilon_k,
            bounds_lower=bounds_lower,
            bounds_upper=bounds_upper,
            refold_to_first_BZ=refold_to_first_BZ
        )

    # plot the fermi surface
    plotFermiSurface2D(
            fermi_surface,
            brillouin_zone=brillouin_zone,
            plot_title=plot_title,
            plot_color=plot_color,
            figsize=figsize,
            showPlot=showPlot,
            save_filename=save_filename
        )

end


# plot the Fermi surface (3D)
function plotFermiSurface3D(
            k_values::Array{Float64, 2};
            brillouin_zone::BrillouinZone=BrillouinZone(),
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,6),
            showPlot::Bool=true,
            save_filename::String="NONE"
        )


    ###########################
    #   INITIAL SETTINGS
    ###########################

    # configure plot environment
    rc("font", family="serif")

    # create a new figure
    fig = figure(figsize=figsize)



    ###########################
    #   PLOT FERMI SURFACE
    ###########################

    # if brillouin zone not empty, plot it as well
    if length(brillouin_zone.points) > 0
        plotBrillouinZone(brillouin_zone, new_figure=false)
    end

    # scatter the points
    scatter3D(k_values[:,1], k_values[:,2], k_values[:,3], color=plot_color)



    ###########################
    #   CONFIGURE AXIS & TITLE
    ###########################

    # set the title
    if plot_title == "AUTO"
        # set the title to an automatically generated title
        title("Fermi surface")
    elseif plot_title == ""
        # do nothing title related
    else
        # set the title to the given title
        title(plot_title)
    end

    # equal axis
    gca()[:set_aspect]("equal")



    ###########################
    #   FINISH THE PLOT
    ###########################

    # tighten the layout
    tight_layout()

    # save the plot
    if save_filename != "NONE"
        # make sure the directory exists
        if contains(save_filename, "/")
    		# get the containing folder
    		folder = save_filename[1:findlast(save_filename, '/')]
    		# build the path to that folder
    		mkpath(folder)
    	end
        # save the plot
        savefig(save_filename)
    end

    # maybe show the plot
    if showPlot
        show()
    end

    # return the figure object
    return fig

end
function plotFermiSurface3D(
            unitcell::Unitcell,
            N_points::Int64;
            fermi_energy::Float64=0.0,
            enforce_hermitian::Bool=false,
            epsilon::Float64=1e-10,
            epsilon_k::Float64=1e-10,
            bounds_lower::Array{Float64,1}=-2*pi.*ones(3),
            bounds_upper::Array{Float64,1}=2*pi.*ones(3),
            refold_to_first_BZ::Bool=true,
            brillouin_zone::BrillouinZone=BrillouinZone(),
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,6),
            showPlot::Bool=true,
            save_filename::String="NONE"
        )

    # check if BZ should be calculated
    if length(brillouin_zone.points) == 0
        brillouin_zone = createBrillouinZone(unitcell)
    end

    # calculate the Fermi surface first
    fermi_surface = getFermiSurface3D(
            unitcell,
            N_points,
            fermi_energy=fermi_energy,
            enforce_hermitian=enforce_hermitian,
            epsilon=epsilon,
            epsilon_k=epsilon_k,
            bounds_lower=bounds_lower,
            bounds_upper=bounds_upper,
            refold_to_first_BZ=refold_to_first_BZ
        )

    # plot the fermi surface
    plotFermiSurface3D(
            fermi_surface,
            brillouin_zone=brillouin_zone,
            plot_title=plot_title,
            plot_color=plot_color,
            figsize=figsize,
            showPlot=showPlot,
            save_filename=save_filename
        )

end





"""
    plotFermiSurface(
                unitcell::Unitcell,
                N_points::Int64
             [; fermi_energy::Float64=0.0,
                enforce_hermitian::Bool=false,
                epsilon::Float64=1e-10,
                epsilon_k::Float64=1e-10,
                bounds_lower::Array{Float64,1}=-2*pi.*ones(4),
                bounds_upper::Array{Float64,1}=2*pi.*ones(4),
                refold_to_first_BZ::Bool=true,
                brillouin_zone::BrillouinZone=BrillouinZone(),
                plot_title::String="",
                plot_color="b",
                figsize::Tuple=(6,6),
                showPlot::Bool=true,
                save_filename::String="NONE" ]
            )

calculates `N_points` points which belong to the Fermi surface of the tight-binding model given by
the `Unitcell` object `unitcell` and plots them using `PyPlot`. The Fermi energy can be adjusted by using `fermi_energy=...`.

The points are located using Newton's method with random starting positions. The procedure can be
modified by using some of the optional keywords:
- `enforce_hermitian` determines if the matrix used for energy calculation is made hermitian by construction
- `epsilon` energy threshold for values to be considered close to the Fermi energy
- `epsilon_k` small k increment that is used for calculation of derivatives
- `bounds_upper` / `bounds_lower` the bounds between which the random starting location for the Newton method is searched

After the search has been completed, the values will optionally be refolded into the first brillouin zone by choosing `refold_to_first_BZ=true`.

For plotting, several more options can be used.



# Examples

```julia-repl
julia> plotFermiSurface(unitcell, 100)
PyPlot.Figure(...)

julia> plotFermiSurface(unitcell, 100, fermi_energy=1.0, showPlot=false)
PyPlot.Figure(...)

```
"""
function plotFermiSurface(
            k_values::Array{Float64, 2};
            brillouin_zone::BrillouinZone=BrillouinZone(),
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,6),
            showPlot::Bool=true,
            save_filename::String="NONE"
        )

    # check which function to pass to
    if size(k_values,2) == 2
        # 2D case
        return plotFermiSurface2D(
                k_values,
                brillouin_zone=brillouin_zone,
                plot_title=plot_title,
                plot_color=plot_color,
                figsize=figsize,
                showPlot=showPlot,
                save_filename=save_filename
            )
    else size(k_values,2) == 3
        # 3D case
        return plotFermiSurface3D(
                k_values,
                brillouin_zone=brillouin_zone,
                plot_title=plot_title,
                plot_color=plot_color,
                figsize=figsize,
                showPlot=showPlot,
                save_filename=save_filename
            )
    end

end
function plotFermiSurface(
            unitcell::Unitcell,
            N_points::Int64;
            fermi_energy::Float64=0.0,
            enforce_hermitian::Bool=false,
            epsilon::Float64=1e-10,
            epsilon_k::Float64=1e-10,
            bounds_lower::Array{Float64,1}=-2*pi.*ones(4),
            bounds_upper::Array{Float64,1}=2*pi.*ones(4),
            refold_to_first_BZ::Bool=true,
            brillouin_zone::BrillouinZone=BrillouinZone(),
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,6),
            showPlot::Bool=true,
            save_filename::String="NONE"
        )

    # check if the unitcell can be put in either method
    if length(unitcell.basis[1]) != length(unitcell.lattice_vectors)
        println("Unitcell has not the same number of lattice vectors as dimensions")
        return zeros(1,1)
    end
    # check if BZ should be calculated
    if length(brillouin_zone.points) == 0
        brillouin_zone = createBrillouinZone(unitcell)
    end
    # check which function to pass to
    if length(unitcell.lattice_vectors) == 2
        # 2D case
        return plotFermiSurface2D(
                unitcell,
                N_points,
                fermi_energy=fermi_energy,
                enforce_hermitian=enforce_hermitian,
                epsilon=epsilon,
                epsilon_k=epsilon_k,
                bounds_lower=bounds_lower,
                bounds_upper=bounds_upper,
                refold_to_first_BZ=refold_to_first_BZ,
                brillouin_zone=brillouin_zone,
                plot_title=plot_title,
                plot_color=plot_color,
                figsize=figsize,
                showPlot=showPlot,
                save_filename=save_filename
            )
    else length(unitcell.lattice_vectors) == 3
        # 3D case
        return plotFermiSurface3D(
                unitcell,
                N_points,
                fermi_energy=fermi_energy,
                enforce_hermitian=enforce_hermitian,
                epsilon=epsilon,
                epsilon_k=epsilon_k,
                bounds_lower=bounds_lower,
                bounds_upper=bounds_upper,
                refold_to_first_BZ=refold_to_first_BZ,
                brillouin_zone=brillouin_zone,
                plot_title=plot_title,
                plot_color=plot_color,
                figsize=figsize,
                showPlot=showPlot,
                save_filename=save_filename
            )
    end

end
export plotFermiSurface
