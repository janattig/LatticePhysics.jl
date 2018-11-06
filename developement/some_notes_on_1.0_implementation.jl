################################################################################
#
#	NEW DEFINITIONS OF LATTICES AND UNITCELLS
#
################################################################################

# usings
using LinearAlgebra





################################################################################
#
#	ABSTRACT TYPES
#
################################################################################

# Bonds
abstract type AbstractBond{L,N} end

# Sites
abstract type AbstractSite{L,D} end

# Unitcells
abstract type AbstractUnitcell{D,N,L,B<:AbstractBond{L,N},S<:AbstractSite{L,D}} end

# Lattices
abstract type AbstractLattice{D,N,L,B<:AbstractBond{L,N},S<:AbstractSite{L,D},U<:AbstractUnitcell{D,N,L,B,S}} end




################################################################################
#
#	INTERFACING / ACCESSING
#	ABSTRACT TYPES
#	(have to be overwritten by concrete types)
#
################################################################################

# BONDS

# default constructor interface
function newBond(index_from::Int64, index_to::Int64, label::L, wrap::NTuple{N,Int64}, ::Type{B}) :: B where {L,N,B<:AbstractBond{L,N}} error("not implemented function 'newBond' for bond type $(B) with label type $(L) and wrap length $(N)") end

# from and to index
function indexFrom(b::AbstractBond{L,N}) :: Int64 where {L,N} error("not implemented interface function 'indexFrom' for bond type $(typeof(b))") end
function indexTo(b::AbstractBond{L,N})   :: Int64 where {L,N} error("not implemented interface function 'indexTo' for bond type $(typeof(b))") end

# label
function label(b::AbstractBond{L,N}) :: L where {L,N} error("not implemented interface function 'label' for bond type $(typeof(b))") end

# wrap
function wrap(b::AbstractBond{L,N}) :: NTuple{N,Int64} where {L,N} error("not implemented interface function 'wrap' for bond type $(typeof(b))") end


# TESTING THE BOND INTERFACE
function testInterface(::Type{T}) where {T<:AbstractBond}
	# get the parameterless constructor
	B = Base.typename(T).wrapper
	# get a new bond
	bond_1 = newBond(1, 1, 1.0, (1,1), B{Float64, 2})
	bond_2 = newBond(1, 1, "tx", (1,1), B{String, 2})
	bond_3 = newBond(1, 1, 2, (1,1,1), B{Int64, 3})
	# get the indices and labels and wraps
	for b in [bond_1, bond_2, bond_3]
		indexFrom(b)
		indexTo(b)
		label(b)
		wrap(b)
	end
	# return nothing
	return nothing
end




# SITES

# default constructor interface
function newSite(position::Vector{Float64}, label::L, ::Type{S}) :: S where {L,S<:AbstractSite{L,D} where D} error("not implemented function 'newSite' for site type $(S) with label type $(L)") end

# position
function point(s::AbstractSite{L,D}) :: Vector{Float64} where {L,D} error("not implemented interface function 'point' for site type $(typeof(s))") end

# label
function label(s::AbstractSite{L,D}) :: L where {L,D} error("not implemented interface function 'label' for site type $(typeof(s))") end


# TESTING THE SITE INTERFACE
function testInterface(::Type{T}) where {T<:AbstractSite}
	# get the parameterless constructor
	S = Base.typename(T).wrapper
	# get a new site
	site_1 = newSite([1.0, 0.0], "site_1", S{String,2})
	site_2 = newSite([1.0, 0.0, 0.0], 1.0, S{Float64,3})
	site_3 = newSite([1.0, 0.0], 1, S{Int64,2})
	# get the points and labels
	for s in [site_1, site_2, site_3]
		point(s)
		label(s)
	end
	# return nothing
	return nothing
end











# DEFINITION OF A BOND
# L - Label type
# N - dimension of the Bravais lattice (i.e. number of primitive directions to wrap)
mutable struct Bond{L,N} <: AbstractBond{L,N}

	# indices that the bond connects
	from	::Int64
	to		::Int64

	# label
	label	::L

	# wrap
	wrap	::NTuple{N, Int64}

end


# IMPLEMENTATION OF BOND RELATED FUNCTIONS

# creation
function newBond(index_from::Int64, index_to::Int64, label::L, wrap::NTuple{N,Int64}, ::Type{Bond{L,N}}) :: Bond{L,N} where {L,N} return Bond{L,N}(index_from, index_to, label, wrap) end

# from and to index
function indexFrom(b::Bond{L,N}) :: Int64 where {N,L} return b.from end
function indexTo(b::Bond{L,N})   :: Int64 where {N,L} return b.to end

# label
function label(b::Bond{L,N}) :: L where {N,L} return b.label end

# wrap
function wrap(b::Bond{L,N}) :: NTuple{N,Int64} where {N,L} return b.wrap end







# DEFINITION OF A SITE
# L - Label type
# D - dimension of the real space in which the site is located
mutable struct Site{L,D} <: AbstractSite{L,D}

	# label
	label	:: L

	# point
	point	:: Vector{Float64}

end



# IMPLEMENTATION OF SITE RELATED FUNCTIONS

# default constructor interface
function newSite(position::Vector{Float64}, label::L, ::Type{Site{L,D}}) :: Site{L,D} where {L,D} return Site{L,D}(label, position) end

# position
function point(s::Site{L,D}) :: Vector{Float64} where {L,D} return s.point end

# label
function label(s::Site{L,D}) :: L where {L,D} return s.label end


















# DEFINITION OF A UNITCELL
# D - space dimension
# N - dimension of the Bravais lattice
# L - label type
mutable struct Unitcell{D,N,L,B,S} <: AbstractUnitcell{D,N,L,B,S}

    # FIELDS

    # basis vectors of the Bravais lattice
    lattice_vectors	::Vector{Vector{Float64}}

    # basis sites within the unitcell
    sites			::Vector{S}

    # list of bonds
    bonds			::Vector{B}

end



# Unitcell obtaining (2D - Square)
function getUnitcellSquare_New()::Unitcell_New{2,2,Float64}
	# return the unitcell
	return Unitcell_New{2,2,Float64}(
		# Bravais lattice vectors
	 	Vector{Float64}[
	    	[1.0, 0.0],
	    	[0.0, 1.0]
		],
	    # Basis Definition
	    Vector{Float64}[
	        [0.0, 0.0]
	    ],
	    # Bond Definition
	    Bond{Float64,2}[
	        Bond(1,1, 1.0, ( 0, 1) ),
	        Bond(1,1, 1.0, ( 0,-1) ),
	        Bond(1,1, 1.0, ( 1, 0) ),
	        Bond(1,1, 1.0, (-1, 0) )
	    ]
	)
end

# Unitcell obtaining (3D - SC)
function getUnitcellSC_New()::Unitcell_New{3,3,Float64}
	# return the unitcell
	return Unitcell_New{3,3,Float64}(
		# Bravais lattice vectors
		Vector{Float64}[
	        [1.0, 0.0, 0.0],
	        [0.0, 1.0, 0.0],
	        [0.0, 0.0, 1.0]
		],
        # Basis Definition
        Vector{Float64}[
            [0.0, 0.0, 0.0]
        ],
		# Bond Definition
	    Bond{Float64,3}[
            Bond(1, 1, 1.0, ( 1, 0, 0)),
            Bond(1, 1, 1.0, (-1, 0, 0)),
            Bond(1, 1, 1.0, (0,  1, 0)),
            Bond(1, 1, 1.0, (0, -1, 0)),
            Bond(1, 1, 1.0, (0, 0,  1)),
            Bond(1, 1, 1.0, (0, 0, -1)),
        ]
	)
end









# DEFINITION OF A LATTICE
# D - space dimension
# N - dimension of the Bravais lattice
# L - label type
mutable struct Lattice_New{D,N,L}

    # FIELDS

	# The unitcell of the lattice if applicable
    unitcell			::Unitcell_New{D,N,L}
	# how often the unitcell is repeated in the elementary directions
    unitcellRepetitions	::Array{Int64, 1}


    # basis vectors of the Bravais lattice
    lattice_vectors		::Array{Vector{Float64}, 1}

    # basis sites within the unitcell
    sites				::Array{Vector{Float64}, 1}
	# array of ints that give the index of the site in the original unit cell
    site_indices		::Array{Int64, 1}

    # list of bonds
    bonds				::Array{Bond{L,N}, 1}

end




# define the index function to get the index of an element in the array
function index(
			i1::Int64, i2::Int64,
			N1::Int64, N2::Int64
		) :: Int64
	return (i1-1)*N2 + i2
end
function index(
			i1::Int64, i2::Int64, i3::Int64,
			N1::Int64, N2::Int64, N3::Int64
		) :: Int64
	return N3*((i1-1)*N2+i2-1) + i3
end
function index(
			i1::Int64, i2::Int64, i3::Int64, i4::Int64,
			N1::Int64, N2::Int64, N3::Int64, N4::Int64
		) :: Int64
	return N4*(N3*((i1-1)*N2+i2-1) + i3 -1) + i4
end

# function to constrain an index
# returns (index_constrained, offset)
# s.t. index = index_constraint + offset*index_max
# To be replaced by ic = div(.,.) and offset = rem(.,.) ???
# apparently not faster...
function getConstrainedIndex(
		index		::Int64,
		index_max	::Int64
	) :: Tuple{Int64, Int64}
	# specify the initial offset as 0
	offset = 0 :: Int64
	# check that the index is bigger than the minimum value
	while index < 1
		index += index_max
		offset += -1
	end
	# check that the index is smaller than the maximum value
	while index > index_max
		index -= index_max
		offset += 1
	end
	# return both index and offset
	return (index, offset)
end




# Construction of periodic lattice (in 2D)
function getLatticePeriodic_New(
			unitcell			::Unitcell_New{D,2,L},
			repetition_array	::Array{Int64, 1}
		) :: Lattice_New{D,2,L} where {D,L}

    # extract the cardinal directions of the lattice from the array
    N_a1	= abs(repetition_array[1]) 	:: Int64
    N_a2	= abs(repetition_array[2]) 	:: Int64
	N_sites = length(unitcell.sites)	:: Int64
	N_bonds = length(unitcell.bonds)	:: Int64

    # GENERATE NEW POSITIONS
	sites 				= Array{Vector{Float64},1}(undef, N_a1*N_a2*N_sites)
    site_indices 		= Array{Int64,1}(undef, N_a1*N_a2*N_sites)

	# set all positions to their correct values
	for i in 1:N_a1
	for j in 1:N_a2
	    @simd for a in 1:N_sites
			# calculate the site index
			site_index 					= index(i,j,a, N_a1,N_a2,N_sites)
			# calculate the position
			sites[site_index]			= unitcell.sites[a] .+ (i.*unitcell.lattice_vectors[1]) .+ (j.*unitcell.lattice_vectors[2])
			# set the site index
	        site_indices[site_index]	= a
		end
	end
    end

    # GENERATE NEW CONNECTIONS
	bonds = Array{Bond{L,2},1}(undef, N_a1*N_a2*N_bonds)

	# iterate over all unit cells
	for i in 1:N_a1
	for j in 1:N_a2
		# add all connections for unitcell (i,j)
		@simd for b in 1:N_bonds
			# get the respective bond
			bond = unitcell.bonds[b] :: Bond{L,2}
			# calculate the index from where the connection goes
			index_from = index(i,j,bond.from, N_a1,N_a2,N_sites)
			# calculate the aimed unitcell
			i_to = i + bond.wrap[1]
			j_to = j + bond.wrap[2]
			# check if the connection goes around in a1 direction
            (i_to, offset_a1) = getConstrainedIndex(i_to, N_a1)
			# check if the connection goes around in a2 direction
            (j_to, offset_a2) = getConstrainedIndex(j_to, N_a2)
			# get the index to where the connection goes
			index_to = index(i_to, j_to, bond.to, N_a1,N_a2,N_sites)
            # generate a new connection and set it in the list
			bonds[index(i,j,b, N_a1,N_a2,N_bonds)] = Bond{L,2}(index_from, index_to, bond.label, (offset_a1, offset_a2))
		end
	end
	end

    # generate new lattice vectors, now spanning the WHOLE lattice
    lattice_vectors = Vector{Float64}[
			# multiply by extent of lattice
			unitcell.lattice_vectors[1] .* N_a1,
			unitcell.lattice_vectors[2] .* N_a2
		]

    # save everything to a Lattice object
    lattice = Lattice_New{D,2,L}(
        	unitcell,
        	Int64[N_a1, N_a2],
        	lattice_vectors,
        	sites,
        	site_indices,
        	bonds
        )

    # return the lattice
    return lattice

end

# Construction of periodic lattice (in 3D)
function getLatticePeriodic_New(
			unitcell			::Unitcell_New{D,3,L},
			repetition_array	::Array{Int64, 1}
		) :: Lattice_New{D,3,L} where {D,L}

    # extract the cardinal directions of the lattice from the array
    N_a1	= abs(repetition_array[1]) 	:: Int64
    N_a2	= abs(repetition_array[2]) 	:: Int64
    N_a3	= abs(repetition_array[3]) 	:: Int64
	N_sites = length(unitcell.sites)	:: Int64
	N_bonds = length(unitcell.bonds)	:: Int64

    # GENERATE NEW POSITIONS
	sites 				= Array{Vector{Float64},1}(undef, N_a1*N_a2*N_a3*N_sites)
    site_indices 		= Array{Int64,1}(undef, N_a1*N_a2*N_a3*N_sites)

	# set all positions to their correct values
	for i in 1:N_a1
	for j in 1:N_a2
	for k in 1:N_a3
	    @simd for a in 1:N_sites
			# calculate the site index
			site_index 					= index(i,j,k,a, N_a1,N_a2,N_a3,N_sites)
			# calculate the position
			sites[site_index]			= unitcell.sites[a] .+ (i.*unitcell.lattice_vectors[1]) .+ (j.*unitcell.lattice_vectors[2]) .+ (k.*unitcell.lattice_vectors[3])
			# set the site index
	        site_indices[site_index]	= a
		end
	end
	end
    end

    # GENERATE NEW CONNECTIONS
	bonds = Array{Bond{L,3},1}(undef, N_a1*N_a2*N_a3*N_bonds)

	# iterate over all unit cells
	for i in 1:N_a1
	for j in 1:N_a2
	for k in 1:N_a3
		# add all connections for unitcell (i,j)
		@simd for b in 1:N_bonds
			# get the respective bond
			bond = unitcell.bonds[b] :: Bond{L,3}
			# calculate the index from where the connection goes
			index_from = index(i,j,k,bond.from, N_a1,N_a2,N_a3,N_sites)
			# calculate the aimed unitcell
			i_to = i + bond.wrap[1] ::Int64
			j_to = j + bond.wrap[2] ::Int64
			k_to = k + bond.wrap[3] ::Int64
			# check if the connection goes around in a1 direction
            (i_to, offset_a1) = getConstrainedIndex(i_to, N_a1)
			# check if the connection goes around in a2 direction
            (j_to, offset_a2) = getConstrainedIndex(j_to, N_a2)
			# check if the connection goes around in a3 direction
            (k_to, offset_a3) = getConstrainedIndex(k_to, N_a3)
			# get the index to where the connection goes
			index_to = index(i_to, j_to, k_to, bond.to, N_a1,N_a2,N_a3,N_sites)
            # generate a new connection and set it in the list
			bonds[index(i,j,k,b, N_a1,N_a2,N_a3,N_bonds)] = Bond{L,3}(index_from, index_to, bond.label, (offset_a1, offset_a2, offset_a3))
		end
	end
	end
	end

	# generate new lattice vectors, now spanning the WHOLE lattice
    lattice_vectors = Vector{Float64}[
			# multiply by extent of lattice
			unitcell.lattice_vectors[1] .* N_a1,
			unitcell.lattice_vectors[2] .* N_a2,
			unitcell.lattice_vectors[3] .* N_a3
		]


    # save everything to a Lattice object
    lattice = Lattice_New{D,3,L}(
        	unitcell,
        	Int64[N_a1, N_a2, N_a3],
        	lattice_vectors,
        	sites,
        	site_indices,
        	bonds
        )

    # return the lattice
    return lattice

end
