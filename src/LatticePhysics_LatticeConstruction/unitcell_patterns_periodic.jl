# CREATE LATTICES FROM UNITCELLS
# PUT TOGEHTER PATTERNS OF UNITCELLS PERIODICALLY



# Construction of periodic lattice (in 2D)
function getLatticePeriodic(
        unitcell        :: U,
        extent          :: NTuple{2,Int64},
        lattice_type    :: Type{L} = Lattice{S,B,U}
    ) :: L where {
        D,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,2},
        U<:AbstractUnitcell{S,B},
        DL,LLS,LLB,SL<:AbstractSite{LLS,DL},BL<:AbstractBond{LLB,2},
        L<:AbstractLattice{SL,BL}
    }

    # extract the cardinal directions of the lattice from the array
    N_a1	= abs(extent[1])       :: Int64
    N_a2	= abs(extent[2])       :: Int64
	N_sites = numSites(unitcell)   :: Int64
	N_bonds = numBonds(unitcell)   :: Int64

    # GENERATE NEW POSITIONS
	site_list  = Array{SL,1}(undef, N_a1*N_a2*N_sites)

	# set all positions to their correct values
	for i in 1:N_a1
	for j in 1:N_a2
	    for a in 1:N_sites
			# calculate the site index
			site_index             = index(i,j,a, N_a1,N_a2,N_sites)
			# calculate the position
			site_position          = point(site(unitcell,a)) .+ (i.*a1(unitcell)) .+ (j.*a2(unitcell))
			# set the site in the list
	        site_list[site_index]  = newSite(SL, site_position, label(site(unitcell,a)))
		end
	end
    end

    # GENERATE NEW CONNECTIONS
	bond_list  = Array{BL,1}(undef, N_a1*N_a2*N_bonds)

	# iterate over all unit cells
	for i in 1:N_a1
	for j in 1:N_a2
		# add all connections for unitcell (i,j)
		for b in 1:N_bonds
			# get the respective bond
			current_bond = bond(unitcell,b) :: BL
			# calculate the site index from where the connection goes
			index_from = index(i,j,from(current_bond), N_a1,N_a2,N_sites)
			# calculate the aimed unitcell
			i_to = i + wrap(current_bond)[1]
			j_to = j + wrap(current_bond)[2]
			# check if the connection goes around in a1 direction
            (i_to, offset_a1) = getConstrainedIndex(i_to, N_a1)
			# check if the connection goes around in a2 direction
            (j_to, offset_a2) = getConstrainedIndex(j_to, N_a2)
			# get the site index to where the connection goes
			index_to = index(i_to, j_to, to(current_bond), N_a1,N_a2,N_sites)
            # generate a new connection and set it in the list
			bond_list[index(i,j,b, N_a1,N_a2,N_bonds)] = newBond(BL, index_from, index_to, label(current_bond), (offset_a1, offset_a2))
		end
	end
	end

    # generate new lattice vectors, now spanning the WHOLE lattice
    lattice_vectors = Vector{Float64}[
			# multiply by extent of lattice
			a1(unitcell) .* N_a1,
			a2(unitcell) .* N_a2
		]

    # save everything to a Lattice object
    lattice = newLattice(
            L,
        	lattice_vectors,
        	site_list,
        	bond_list,
            unitcell
        )

    # return the lattice
    return lattice

end






# Construction of periodic lattice (in 3D)
function getLatticePeriodic(
        unitcell        :: U,
        extent          :: NTuple{3,Int64},
        lattice_type    :: Type{L} = Lattice{S,B,U}
    ) :: L where {
        D,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,3},
        U<:AbstractUnitcell{S,B},
        DL,LLS,LLB,SL<:AbstractSite{LLS,DL},BL<:AbstractBond{LLB,3},
        L<:AbstractLattice{SL,BL}
    }

    # extract the cardinal directions of the lattice from the array
    N_a1	= abs(extent[1])       :: Int64
    N_a2	= abs(extent[2])       :: Int64
    N_a3	= abs(extent[3])       :: Int64
	N_sites = numSites(unitcell)   :: Int64
	N_bonds = numBonds(unitcell)   :: Int64

    # GENERATE NEW POSITIONS
	site_list  = Array{SL,1}(undef, N_a1*N_a2*N_a3*N_sites)

	# set all positions to their correct values
	for i in 1:N_a1
	for j in 1:N_a2
	for k in 1:N_a3
	    for a in 1:N_sites
			# calculate the site index
			site_index             = index(i,j,k,a, N_a1,N_a2,N_a3,N_sites)
			# calculate the position
			site_position          = point(site(unitcell,a)) .+ (i.*a1(unitcell)) .+ (j.*a2(unitcell)) .+ (k.*a3(unitcell))
			# set the site in the list
	        site_list[site_index]  = newSite(SL, site_position, label(site(unitcell,a)))
		end
	end
    end
    end

    # GENERATE NEW CONNECTIONS
	bond_list  = Array{BL,1}(undef, N_a1*N_a2*N_a3*N_bonds)

	# iterate over all unit cells
	for i in 1:N_a1
	for j in 1:N_a2
	for k in 1:N_a3
		# add all connections for unitcell (i,j)
		for b in 1:N_bonds
			# get the respective bond
			current_bond = bond(unitcell, b) :: BL
			# calculate the site index from where the connection goes
			index_from = index(i,j,k, from(current_bond), N_a1,N_a2,N_a3,N_sites)
			# calculate the aimed unitcell
			i_to = i + wrap(current_bond)[1]
			j_to = j + wrap(current_bond)[2]
			k_to = k + wrap(current_bond)[3]
			# check if the connection goes around in a1 direction
            (i_to, offset_a1) = getConstrainedIndex(i_to, N_a1)
			# check if the connection goes around in a2 direction
            (j_to, offset_a2) = getConstrainedIndex(j_to, N_a2)
			# check if the connection goes around in a3 direction
            (k_to, offset_a3) = getConstrainedIndex(k_to, N_a3)
			# get the site index to where the connection goes
			index_to = index(i_to,j_to,k_to, to(current_bond), N_a1,N_a2,N_a3,N_sites)
            # generate a new connection and set it in the list
			bond_list[index(i,j,k,b, N_a1,N_a2,N_a3,N_bonds)] = newBond(BL, index_from, index_to, label(current_bond), (offset_a1, offset_a2, offset_a3))
		end
	end
	end
	end

    # generate new lattice vectors, now spanning the WHOLE lattice
    lattice_vectors = Vector{Float64}[
			# multiply by extent of lattice
			a1(unitcell) .* N_a1,
			a2(unitcell) .* N_a2,
			a3(unitcell) .* N_a3
		]

    # save everything to a Lattice object
    lattice = newLattice(
            L,
        	lattice_vectors,
        	site_list,
        	bond_list,
            unitcell
        )

    # return the lattice
    return lattice

end




# short hand notation (only pass N instead of (N1,N2,N3))

# 2d
function getLatticePeriodic(
        unitcell        :: U,
        extent          :: Int64,
        lattice_type    :: Type{L} = Lattice{S,B,U}
    ) :: L where {D,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,2},U<:AbstractUnitcell{S,B},L<:AbstractLattice{S,B,U}}

    # return the suitable function
    return getLatticePeriodic(unitcell, (extent, extent), lattice_type)
end

# 3d
function getLatticePeriodic(
        unitcell        :: U,
        extent          :: Int64,
        lattice_type    :: Type{L} = Lattice{S,B,U}
    ) :: L where {D,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,3},U<:AbstractUnitcell{S,B},L<:AbstractLattice{S,B,U}}

    # return the suitable function
    return getLatticePeriodic(unitcell, (extent, extent, extent), lattice_type)
end
