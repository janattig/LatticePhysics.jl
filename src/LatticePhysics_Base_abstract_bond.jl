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
    error(  "not implemented function 'newBond' for concrete bond type " *
            string(B) * " with label type " * string(L) *
            " and wrap length " * string(N)   )
end



# UIDs of sites between which the bond is located

# from index / UID (Int64)
function from(
            b :: AbstractBond{L,N}
        ) :: Int64 where {L,N}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'indexFrom' for bond type " * string(typeof(b)))
end

# to index / UID (Int64)
function to(
            b :: AbstractBond{L,N}
        ) :: Int64 where {L,N}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'indexTo' for bond type " * string(typeof(b)))
end







# label
function label(
            b :: AbstractBond{L,N}
        ) :: L where {L,N}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'label' for bond type " * string(typeof(b)))
end



# wrap
function wrap(
            b :: AbstractBond{L,N}
        ) :: NTuple{N,Int64} where {L,N}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'wrap' for bond type " * string(typeof(b)))
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

    # iterate over some standard wraps
    for w in [(1,), (0,0), (0,0,0), (1,2,0,0)]
    # iterate over some standard labels
    for l in ["t", 1, 1.0]
        # create a new bond
        bond = newBond(1, 1, l, w, B{typeof(l), length(w)})
        # test the interface
		indexFrom(b)
		indexTo(b)
		label(b)
		wrap(b)
    end
    end
	# return true to indicate the test passed
	return true
end
