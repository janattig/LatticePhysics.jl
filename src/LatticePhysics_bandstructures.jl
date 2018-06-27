################################################################################
#
#   METHODS FOR CONSTRUCTION OF BANDSTRUCTURES (ALONG PATHS)
#
#   STRUCTURE OF THE FILE
#
#   TYPE BANDSTRUCTURE
#       - type definition
#       - printInfo function
#
################################################################################




################################################################################
#
#   TYPE BANDSTRUCTURE
#       - type definition
#       - printInfo function
#
################################################################################
struct Bandstructure

    # the path along which the band structure is calcualted
    path::Path

    # bands for each segment
    # bands[i] gives all bands of segment i
    # bands[i][j] gives all energy values for band j of segment i
    # bands[i][j][k] gives the energy value at kpoint index k of band j in segment i
    bands::Array{Array{Array{Float64, 1}, 1}, 1}

    # ONLY DEFAULT CONSTRUCTOR

end
