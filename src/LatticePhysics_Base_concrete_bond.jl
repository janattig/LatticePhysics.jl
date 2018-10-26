################################################################################
#
#	CONCRETE TYPE
#
#   Bond{L,N} <: AbstractBond{L,N}
#   --> L is the label type
#   --> N is the dimension of the Bravais lattice that the bond is located in
#
#   FILE CONTAINS
#       - concrete struct definition
#       - interface implementation
#
################################################################################






################################################################################
#
#   CONCRETE STRUCT DEFINITION
#
################################################################################
mutable struct Bond{L,N} <: AbstractBond{L,N}

	# indices that the bond connects
	from	:: Int64
	to		:: Int64

	# label
	label	:: L

	# wrap
	wrap	:: NTuple{N, Int64}

end




################################################################################
#
#	IMPLEMENTATION OF INTERFACE FOR CONCRETE BOND TYPE
#	(functions that had to be overwritten by concrete type)
#
################################################################################

# default constructor interface
# used for creation of new bonds
function newBond(
            from    :: Int64,
            to      :: Int64,
            label   :: L,
            wrap    :: NTuple{N,Int64},
            :: Type{Bond{L,N}}
        ) :: Bond{L,N} where {L,N}

    # return the constructed Bond
    return Bond{L,N}(from, to, label, wrap)
end



# from index / UID (Int64)
function from(
            b :: Bond{L,N}
        ) :: Int64 where {L,N}

    # return the index
    return b.from
end

# to index / UID (Int64)
function to(
            b :: Bond{L,N}
        ) :: Int64 where {L,N}

    # return the index
    return b.to
end



# label
function label(
            b :: Bond{L,N}
        ) :: L where {L,N}

    # return the label
    return b.label
end


# wrap
function wrap(
            b :: Bond{L,N}
        ) :: NTuple{N,Int64} where {L,N}

    # return the wrap
    return b.wrap
end
