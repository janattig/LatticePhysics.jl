
function getKPointsForZeroEnergy(lattice::Lattice, n::Int64; enforce_hermitian=false, epsilon=1e-10, epsilon_k=1e-10, slowdown_factor=0.75, bounds_lower=-2*pi.*ones(4), bounds_upper=2*pi.*ones(4), dumpEvery=-1, dumpFilename="NONE")

    # obtain the dimension
    dimension = size(lattice.positions[1],1)
    # list of k values
    k_values = Array[]
    # funtion of energy
    if enforce_hermitian
        function energy(kvector)
            evals = calculateEnergySingleKHermitian(lattice, kvector)
            evals = evals .* evals
            return minimum(evals)
        end
    else
        function energy(kvector)
            evals = calculateEnergySingleK(lattice, kvector)
            evals = evals .* evals
            return minimum(evals)
        end
    end

    # print the start
    #println("starting the calculation of $(n) points")

    # search until there are enough points
    while size(k_values,1) < n
        # maybe dump
        if dumpEvery>0 && size(k_values,1)%dumpEvery == 0 && dumpFilename!="NONE"
            dumpKPointListToFolder(lattice, k_values, "", dumpFilename)
        end
        # find a suitable starting point
        k = zeros(dimension)
        for j in 1:length(k)
            gamma = rand()
            k[j] = gamma * bounds_lower[j] + (1-gamma) * bounds_upper[j]
        end
        # small increment and derivatives
        e0 = energy(k)
        # iterate
        for i in 1:100
            if e0 < epsilon
                push!(k_values, k)
                # break the newton loop
                break
            end
            H_0 = e0
            if dimension == 3
                H_eps =
                [
                    energy(k .+ [epsilon_k, 0, 0]),
                    energy(k .+ [0, epsilon_k, 0]),
                    energy(k .+ [0, 0, epsilon_k])
                ]
            elseif dimension == 2
                H_eps =
                [
                    energy(k .+ [epsilon_k, 0]),
                    energy(k .+ [0, epsilon_k])
                ]
            end
            dH = (H_eps - H_0) ./ epsilon_k
            dHdH = dot(dH, dH)
            if abs(dHdH) < 1e-20 || abs(dHdH) > 1e20
                break
            end
            dk = dH * (H_0 / dHdH)
            k -= dk
            e0 = energy(k)
        end
    end

    #println("finished the calculation of $(n) points")

    # return the list
    return k_values

end
