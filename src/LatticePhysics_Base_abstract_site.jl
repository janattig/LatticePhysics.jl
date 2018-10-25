################################################################################
#
#	ABSTRACT TYPE
#
#   Site{L,D}
#   --> L is the label type
#   --> D is the dimension of embedding space
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
abstract type AbstractSite{L,D} end




################################################################################
#
#	INTERFACING / ACCESSING BONDS
#	(functions have to be overwritten by concrete types)
#
################################################################################

# default constructor interface
# used for creation of new sites
function newSite(
            point   :: Vector{<:Real},
            label   :: L,
            :: Type{S}
        ) :: S where {L,S<:AbstractSite{L,D} where D}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'newSite' for concrete site type " *
            string(S) * " with label type " * string(L) )
end





# label
function label(
            s :: AbstractSite{L,D}
        ) :: L where {L,D}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'label' for site type " * string(typeof(s)))
end



# point
function point(
            s :: AbstractSite{L,D}
        ) :: Vector{Float64} where {L,D}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'point' for site type " * string(typeof(s)))
end




################################################################################
#
#	TESTING THE INTERFACE OF SITES
#	(should never be overwritten by concrete types)
#
################################################################################

# TESTING THE SITE INTERFACE
function testInterface(
            ::Type{T}
        ) :: Bool where {T<:AbstractSite}

	# get the parameterless constructor
	S = Base.typename(T).wrapper

    # iterate over some standard points
    for p in Vector{Float64}[[1.0,], [1.0, 1.0], [0.0, 0.0, 0.0, 0.0]]
    # iterate over some standard labels
    for l in ["t", 1, 1.0]
        # create a new site
        s = newSite(p,l, S{typeof(l), length(p)})
        # test the interface
		label(s)
		point(s)
    end
    end

	# return true to indicate the test passed
	return true
end
