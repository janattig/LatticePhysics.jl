################################################################################
#
#	ABSTRACT TYPE
#
#   Bond{L,N}
#   --> L is the label type
#   --> N is the dimension of the Bravais lattice that the bond is located in
#
#   FILE CONTAINS
#       - abstract type definition
#       - interface definition
#       - interface testing
#
################################################################################





################################################################################
#
#   ABSTRACT TYPE DEFINITION
#
################################################################################
abstract type AbstractBond{L,N} end




################################################################################
#
#	INTERFACING / ACCESSING BONDS
#	(functions have to be overwritten by concrete types)
#
################################################################################

# default constructor interface
# used for creation of new bonds
function newBond(
            from    :: Int64,
            to      :: Int64,
            label   :: L,
            wrap    :: NTuple{N,Int64},
            :: Type{B}
        ) :: B where {L,N,B<:AbstractBond{L,N}}

    # print an error because implementation for concrete type is missing
    error("not implemented function 'newBond' for concrete bond type $(B) with label type $(L) and wrap length $(N)")
end



# UIDs of sites between which the bond is located

# from index / UID (Int64)
function from(
            b :: AbstractBond{L,N}
        ) :: Int64 where {L,N}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'indexFrom' for bond type $(typeof(b))")
end

# to index / UID (Int64)
function to(
            b :: AbstractBond{L,N}
        ) :: Int64 where {L,N}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'indexTo' for bond type $(typeof(b))")
end







# label
function label(
            b :: AbstractBond{L,N}
        ) :: L where {L,N}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'label' for bond type $(typeof(b))")
end



# wrap
function wrap(
            b :: AbstractBond{L,N}
        ) :: NTuple{N,Int64} where {L,N}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'wrap' for bond type $(typeof(b))")
end




################################################################################
#
#	INTERFACING / ACCESSING BONDS
#	(should never be overwritten by concrete types)
#
################################################################################

# TESTING THE BOND INTERFACE
function testInterface(
            ::Type{T}
        ) :: Bool where {T<:AbstractBond}

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

	# return true to indicate the test passed
	return true
end
