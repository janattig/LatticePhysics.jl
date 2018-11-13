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
#	INTERFACING / ACCESSING SITES
#	(functions have to be overwritten by concrete types)
#
################################################################################

# default constructor interface
# used for creation of new sites
function newSite(
            :: Type{S},
            point   :: Vector{<:Real},
            label   :: L
        ) :: S where {L,S<:AbstractSite{L,D} where D}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'newSite' for concrete site type " *
            string(S) * " with label type " * string(L) )
end





# get label
function label(
            s :: AbstractSite{L,D}
        ) :: L where {L,D}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'label' for site type " * string(typeof(s)))
end
# set label
function label!(
            s :: AbstractSite{L,D},
            l :: L
        ) where {L,D}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'label!' for site type " * string(typeof(s)))
end



# get point
function point(
            s :: AbstractSite{L,D}
        ) :: Vector{Float64} where {L,D}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'point' for site type " * string(typeof(s)))
end
# set point
function point!(
            s :: AbstractSite{L,D},
            p :: Vector{<:Real}
        ) where {L,D}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'point!' for site type " * string(typeof(s)))
end





# SIMILAR FUNCTION (can be overwritten but does not have to be overwritten)

# without new parameters
function similar(
            s :: S
        ) :: S where {L,D,S<:AbstractSite{L,D}}

    # return a new site object
    return newSite(S, deepcopy(point(s)), deepcopy(label(s)))
end
# with new parameters
function similar(
            s :: S,
            p :: Vector{<:Real},
            l :: L
        ) :: S where {L,D,S<:AbstractSite{L,D}}

    # create a new site object
    s_new = similar(s)
    # set parameters
    point!(s_new, p)
    label!(s_new, l)
    # return the new object
    return s_new
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
        s = newSite(S{typeof(l), length(p)}, p,l)
        # test the interface by getting label and point
		label(s)
		point(s)
        # set both label and point
        label!(s, l)
        point!(s, p)
    end
    end

	# return true to indicate the test passed
	return true
end
