################################################################################
#
#	CONCRETE TYPE
#
#   Site{L,D} <: AbstractSite{L,D}
#   --> L is the label type
#   --> D is the dimension of embedding space
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
mutable struct Site{L,D} <: AbstractSite{L,D}

	# point
	point	:: Vector{Float64}

	# label
	label	:: L

end





################################################################################
#
#	IMPLEMENTATION OF INTERFACE FOR CONCRETE SITE TYPE
#	(functions that had to be overwritten by concrete type)
#
################################################################################

# default constructor interface
# used for creation of new sites
function newSite(
            point   :: Vector{<:Real},
            label   :: L,
            :: Type{Site{L,D}}
        ) :: Site{L,D} where {L,D}

    # check if correct dimension is given
    @assert length(point) == D

    # return a new object
    return Site{L,D}(point, label)
end





# label
function label(
            s :: Site{L,D}
        ) :: L where {L,D}

    # return the label
    return s.label
end


# point
function point(
            s :: Site{L,D}
        ) :: Vector{Float64} where {L,D}

    # return the point
    return s.point
end
