# obtain Fermi surface
function getFermiSurface2D(
            unitcell::Unitcell,
            N_points::Int64;
            fermi_energy::Float64=0.0,
            enforce_hermitian::Bool=false,
            epsilon::Float64=1e-10,
            epsilon_k::Float64=1e-10,
            slowdown_factor::Float64=0.75,
            bounds_lower::Array{Float64,1}=-2*pi.*ones(4),
            bounds_upper::Array{Float64,1}=2*pi.*ones(4),
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
