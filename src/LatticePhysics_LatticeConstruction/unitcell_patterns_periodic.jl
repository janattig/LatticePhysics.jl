# CREATE LATTICES FROM UNITCELLS
# PUT TOGEHTER PATTERNS OF UNITCELLS PERIODICALLY



# Construction of periodic lattice (in 2D)
function getLatticePeriodic(
        unitcell        :: U,
        extent          :: NTuple{2,Int64},
        lattice_type    :: Type{LA} = Lattice{D,2,L,S,B}
    ) :: LA where {D,L,S<:AbstractSite{L,D},B<:AbstractBond{L,2},U<:AbstractUnitcell{D,2,L,S,B},LA<:AbstractLattice{D,2,L,S,B}}

    # extract the cardinal directions of the lattice from the array
    N_a1	= abs(extent[1])           :: Int64
    N_a2	= abs(extent[2])           :: Int64
	N_sites = length(sites(unitcell))  :: Int64
	N_bonds = length(bonds(unitcell))  :: Int64

    # GENERATE NEW POSITIONS
	site_list  = Array{S,1}(undef, N_a1*N_a2*N_sites)

	# set all positions to their correct values
	for i in 1:N_a1
	for j in 1:N_a2
	    @simd for a in 1:N_sites
			# calculate the site index
			site_index             = index(i,j,a, N_a1,N_a2,N_sites)
			# calculate the position
			site_position          = point(sites(unitcell)[a]) .+ (i.*latticeVectors(unitcell)[1]) .+ (j.*latticeVectors(unitcell)[2])
			# set the site in the list
	        site_list[site_index]  = newSite(site_position, label(sites(unitcell)[a]), S)
		end
	end
    end

    # GENERATE NEW CONNECTIONS
	bond_list  = Array{B,1}(undef, N_a1*N_a2*N_bonds)

	# iterate over all unit cells
	for i in 1:N_a1
	for j in 1:N_a2
		# add all connections for unitcell (i,j)
		@simd for b in 1:N_bonds
			# get the respective bond
			bond = bonds(unitcell)[b] :: B
			# calculate the site index from where the connection goes
			index_from = index(i,j,from(bond), N_a1,N_a2,N_sites)
			# calculate the aimed unitcell
			i_to = i + wrap(bond)[1]
			j_to = j + wrap(bond)[2]
			# check if the connection goes around in a1 direction
            (i_to, offset_a1) = getConstrainedIndex(i_to, N_a1)
			# check if the connection goes around in a2 direction
            (j_to, offset_a2) = getConstrainedIndex(j_to, N_a2)
			# get the site index to where the connection goes
			index_to = index(i_to, j_to, to(bond), N_a1,N_a2,N_sites)
            # generate a new connection and set it in the list
			bond_list[index(i,j,b, N_a1,N_a2,N_bonds)] = newBond(index_from, index_to, label(bond), (offset_a1, offset_a2), B)
		end
	end
	end

    # generate new lattice vectors, now spanning the WHOLE lattice
    lattice_vectors = Vector{Float64}[
			# multiply by extent of lattice
			latticeVectors(unitcell)[1] .* N_a1,
			latticeVectors(unitcell)[2] .* N_a2
		]

    # save everything to a Lattice object
    lattice = newLattice(
        	lattice_vectors,
        	site_list,
        	bond_list,
            LA
        )

    # return the lattice
    return lattice

end

# Construction of periodic lattice (in 3D)
function getLatticePeriodic(
        unitcell        :: U,
        extent          :: NTuple{3,Int64},
        lattice_type    :: Type{LA} = Lattice{D,3,L,S,B}
    ) :: LA where {D,L,S<:AbstractSite{L,D},B<:AbstractBond{L,3},U<:AbstractUnitcell{D,3,L,S,B},LA<:AbstractLattice{D,3,L,S,B}}

    # extract the cardinal directions of the lattice from the array
    N_a1	= abs(extent[1])           :: Int64
    N_a2	= abs(extent[2])           :: Int64
    N_a3	= abs(extent[3])           :: Int64
	N_sites = length(sites(unitcell))  :: Int64
	N_bonds = length(bonds(unitcell))  :: Int64

    # GENERATE NEW POSITIONS
	site_list  = Array{S,1}(undef, N_a1*N_a2*N_a3*N_sites)

	# set all positions to their correct values
	for i in 1:N_a1
	for j in 1:N_a2
	for k in 1:N_a3
	    @simd for a in 1:N_sites
			# calculate the site index
			site_index             = index(i,j,k,a, N_a1,N_a2,N_a3,N_sites)
			# calculate the position
			site_position          = point(sites(unitcell)[a]) .+ (i.*latticeVectors(unitcell)[1]) .+ (j.*latticeVectors(unitcell)[2]) .+ (k.*latticeVectors(unitcell)[3])
			# set the site in the list
	        site_list[site_index]  = newSite(site_position, label(sites(unitcell)[a]), S)
		end
	end
    end
    end

    # GENERATE NEW CONNECTIONS
	bond_list  = Array{B,1}(undef, N_a1*N_a2*N_a3*N_bonds)

	# iterate over all unit cells
	for i in 1:N_a1
	for j in 1:N_a2
	for k in 1:N_a3
		# add all connections for unitcell (i,j)
		@simd for b in 1:N_bonds
			# get the respective bond
			bond = bonds(unitcell)[b] :: B
			# calculate the site index from where the connection goes
			index_from = index(i,j,k, from(bond), N_a1,N_a2,N_a3,N_sites)
			# calculate the aimed unitcell
			i_to = i + wrap(bond)[1]
			j_to = j + wrap(bond)[2]
			k_to = k + wrap(bond)[3]
			# check if the connection goes around in a1 direction
            (i_to, offset_a1) = getConstrainedIndex(i_to, N_a1)
			# check if the connection goes around in a2 direction
            (j_to, offset_a2) = getConstrainedIndex(j_to, N_a2)
			# check if the connection goes around in a3 direction
            (k_to, offset_a3) = getConstrainedIndex(k_to, N_a3)
			# get the site index to where the connection goes
			index_to = index(i_to,j_to,k_to, to(bond), N_a1,N_a2,N_a3,N_sites)
            # generate a new connection and set it in the list
			bond_list[index(i,j,k,b, N_a1,N_a2,N_a3,N_bonds)] = newBond(index_from, index_to, label(bond), (offset_a1, offset_a2, offset_a3), B)
		end
	end
	end
	end

    # generate new lattice vectors, now spanning the WHOLE lattice
    lattice_vectors = Vector{Float64}[
			# multiply by extent of lattice
			latticeVectors(unitcell)[1] .* N_a1,
			latticeVectors(unitcell)[2] .* N_a2,
			latticeVectors(unitcell)[3] .* N_a3
		]

    # save everything to a Lattice object
    lattice = newLattice(
        	lattice_vectors,
        	site_list,
        	bond_list,
            LA
        )

    # return the lattice
    return lattice

end




# short hand notation (only pass N instead of (N1,N2,N3))

# 2d
function getLatticePeriodic(
        unitcell        :: U,
        extent          :: Int64,
        lattice_type    :: Type{LA} = Lattice{D,2,L,S,B}
    ) :: LA where {D,L,S<:AbstractSite{L,D},B<:AbstractBond{L,2},U<:AbstractUnitcell{D,2,L,S,B},LA<:AbstractLattice{D,2,L,S,B}}

    # return the suitable function
    return getLatticePeriodic(unitcell, (extent, extent), lattice_type)
end

# 3d
function getLatticePeriodic(
        unitcell        :: U,
        extent          :: Int64,
        lattice_type    :: Type{LA} = Lattice{D,3,L,S,B}
    ) :: LA where {D,L,S<:AbstractSite{L,D},B<:AbstractBond{L,3},U<:AbstractUnitcell{D,3,L,S,B},LA<:AbstractLattice{D,3,L,S,B}}

    # return the suitable function
    return getLatticePeriodic(unitcell, (extent, extent, extent), lattice_type)
end
