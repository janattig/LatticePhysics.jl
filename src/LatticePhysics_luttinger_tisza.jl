










################################################################################
#
#   TYPE LTBANDSTRUCTURE
#       - type definition
#       - TODO printInfo function
#
################################################################################
struct LTBandstructure

    # the path along which the band structure is calcualted
    path::Path

    # bands for each segment
    # bands[i] gives all bands of segment i
    # bands[i][j] gives all energy values for band j of segment i
    # bands[i][j][k] gives the energy value at kpoint index k of band j in segment i
    bands::Array{Array{Array{Float64, 1}, 1}, 1}

    # constraint value for all bands
    # value is the minimum of all sum(|s_i - 1.0|^2) for constructed s_i
    constraint_value::Array{Array{Array{Float64, 1}, 1}, 1}

    # ONLY DEFAULT CONSTRUCTOR

end


# export the type
export LTBandstructure
