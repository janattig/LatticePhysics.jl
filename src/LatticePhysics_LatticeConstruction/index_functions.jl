# FUNCTIONS THAT DEFINE INDICES

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
