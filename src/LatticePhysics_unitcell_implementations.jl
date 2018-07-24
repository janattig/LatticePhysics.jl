################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT UNITCELLS in 2D and 3D
#
#   STRUCTURE OF THE FILE
#
#   2D UNITCELLS
#   - SQUARE / RECTANGLE
#       - getUnitcellSquare
#       - getUnitcellExtendedSquare
#       - getUnitcellCheckerboard
#       - getUnitcellShastrySutherland
#       - getUnitcellAdvancedSquare
#       - getUnitcellSquareOctagon
#       - getUnitcellBCC2D
#       - getUnitcellFullyConnectedSquare
#   - TRIANGULAR
#       - getUnitcellTriangular
#       - getUnitcellHoneycomb
#       - getUnitcellKagome
#       - getUnitcellKagomeMinus
#       - getUnitcellHoneycombXXX
#
#   3D UNITCELLS
#   - CUBIC / FCC
#       - getUnitcellDiamond
#       - getUnitcellFCC
#       - getUnitcellBCC
#       - getUnitcellPyrochlore
#   - (X,3)y FAMILY
#       - getUnitcell_8_3_a
#       - getUnitcell_8_3_b
#       - getUnitcell_8_3_c
#       - getUnitcell_8_3_n
#       - getUnitcell_9_3_a
#       - getUnitcell_10_3_a / getUnitcellHyperoctagon
#       - getUnitcell_10_3_b / getUnitcellHyperhoneycomb
#       - getUnitcell_10_3_c
#       - getUnitcell_10_3_d
#
################################################################################





#-----------------------------------------------------------------------------------------------------------------------------
#
#   Individual FUNCTIONS FOR 2D UNITCELLS
#
#   They all generate the special filenames for the unitcells and return the unitcells as objects
#   Functions can be given a version integer to distinguish several implementations of the same lattice
#   Functions can be specified to already save the unitcell to file
#
#-----------------------------------------------------------------------------------------------------------------------------



#----------------------------
#
#   2D SQUARE BASIS
#
#----------------------------

"""
    getUnitcellSquare([version::Int64=1; save::Bool=false])

get the implementation of the *2D square lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 - simple (DEFAULT)

Bravais lattice vectors are

    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]

1 site per unitcell, located at

    r1 = [0.0, 0.0]

All couplings have strength `1.0`



#### 2 - BCC in 2D

Bravais lattice vectors are

    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]

2 sites per unitcell, located at

    r1 = [0.0, 0.0]
    r2 = [0.5, 0.5]

All couplings have strength `1.0`



#### 3 - advanced

Bravais lattice vectors are

    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]

4 sites per unitcell, located at

    r1 = [0.0, 0.0]
    r2 = [0.5, 0.0]
    r3 = [0.5, 0.5]
    r4 = [0.0, 0.5]

couplings are `"t1"` and `"t2"` which are in a staggered square pattern.



# Examples

```julia-repl
julia> unitcell = getUnitcellSquare()
LatticePhysics.Unitcell(...)

julia> unitcell = getUnitcellSquare(2)
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellSquare(version::Int64=1; save::Bool=false)
    # SIMPLE SQUARE LATTICE
    if version == 1
        # the lattice vectors
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 1; 1.0; (0, 1)],
            [1; 1; 1.0; (0, -1)],
            [1; 1; 1.0; (1, 0)],
            [1; 1; 1.0; (-1, 0)],
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_square_unitcell.jld"
    # SQUARE LATTICE THAT LOOKS LIKE BCC in 2D
    elseif version == 2
        # the lattice vectors
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0],
            [0.5, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0, 0)],
            [1; 2; 1.0; (-1, 0)],
            [1; 2; 1.0; (0, -1)],
            [1; 2; 1.0; (-1, -1)],

            [2; 1; 1.0; (0, 0)],
            [2; 1; 1.0; (1, 0)],
            [2; 1; 1.0; (0, 1)],
            [2; 1; 1.0; (1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_unitcell_bcc_2d.jld"
    # ADVANCED SQUARE LATTICE
    elseif version == 3
        # the lattice vectors
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0],
            [0.5, 0.0],
            [0.5, 0.5],
            [0.0, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; "t2"; (-1,0)],
            [1; 4; "t2"; (0,-1)],
            [1; 2; "t1"; (0, 0)],
            [1; 4; "t1"; (0, 0)],

            [2; 1; "t2"; (1, 0)],
            [2; 3; "t1"; (0, 0)],
            [2; 1; "t1"; (0, 0)],
            [2; 3; "t2"; (0,-1)],

            [3; 4; "t1"; (0, 0)],
            [3; 2; "t1"; (0, 0)],
            [3; 4; "t2"; (1, 0)],
            [3; 2; "t2"; (0, 1)],

            [4; 3; "t1"; (0, 0)],
            [4; 1; "t2"; (0, 1)],
            [4; 3; "t2"; (-1,0)],
            [4; 1; "t1"; (0, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_advanced_square_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellSquare




"""
    getUnitcellExtendedSquare([version::Int64=1; save::Bool=false])

get the implementation of the *2D extended square lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 - simple (DEFAULT)

Bravais lattice vectors are

    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]

3 sites per unitcell, located at

    r1 = [0.0, 0.0]
    r2 = [0.5, 0.0]
    r3 = [0.0, 0.5]

All couplings have strength `1.0`



# Examples

```julia-repl
julia> unitcell = getUnitcellExtendedSquare()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellExtendedSquare(version::Int64=1; save::Bool=false)
    # SIMPLE EXTENDED SQUARE LATTICE
    if version == 1
        # the lattice vectors
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0],
            [0.5, 0.0],
            [0.0, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0, 0)],
            [1; 2; 1.0; (-1, 0)],
            [1; 3; 1.0; (0, 0)],
            [1; 3; 1.0; (0, -1)],

            [2; 1; 1.0; (0, 0)],
            [2; 1; 1.0; (1, 0)],

            [3; 1; 1.0; (0, 0)],
            [3; 1; 1.0; (0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_extended_square_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellExtendedSquare



"""
    getUnitcellCheckerboard([version::Int64=1; save::Bool=false])

get the implementation of the *2D checkerboard lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 - simple (DEFAULT)

Bravais lattice vectors are

    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]

2 sites per unitcell, located at

    r1 = [0.0, 0.0]
    r2 = [0.5, 0.5]

All couplings have strength `1.0`



# Examples

```julia-repl
julia> unitcell = getUnitcellCheckerboard()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellCheckerboard(version::Int64=1; save::Bool=false)
    # SIMPLE CHECKERBOARD LATTICE
    if version == 1
        # the lattice vectors
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0],
            [0.5, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 1; 1.0; (1, 0)],
            [1; 1; 1.0; (-1, 0)],
            [1; 2; 1.0; (0, 0)],
            [1; 2; 1.0; (0, -1)],
            [1; 2; 1.0; (-1, 0)],
            [1; 2; 1.0; (-1, -1)],

            [2; 2; 1.0; (0, 1)],
            [2; 2; 1.0; (0, -1)],
            [2; 1; 1.0; (0, 0)],
            [2; 1; 1.0; (0, 1)],
            [2; 1; 1.0; (1, 0)],
            [2; 1; 1.0; (1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_checkerboard_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellCheckerboard


"""
    getUnitcellShastrySutherland([version::Int64=1; save::Bool=false])

get the implementation of the *2D Shastry-Sutherland lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 4 - simple (DEFAULT) & simple (extended) Kitaev

Bravais lattice vectors are

    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]

4 sites per unitcell, located at

    r1 = [0.0, 0.0]
    r2 = [0.5, 0.0]
    r3 = [0.5, 0.5]
    r4 = [0.0, 0.5]

For version `1`, the couplings all have strength `1.0`.

For version `4`, the coupling values are `["t1", ..., "t5"]` for the distinct bonds.



# Examples

```julia-repl
julia> unitcell = getUnitcellShastrySutherland()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellShastrySutherland(version::Int64=1; save::Bool=false)
    # SIMPLE SHASTRY SUTHERLAND LATTICE
    if version == 1
        # the lattice vectors
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0],
            [0.5, 0.0],
            [0.5, 0.5],
            [0.0, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 3; 1.0; (0, 0)],
            [1; 2; 1.0; (-1,0)],
            [1; 4; 1.0; (0,-1)],
            [1; 2; 1.0; (0, 0)],
            [1; 4; 1.0; (0, 0)],

            [2; 4; 1.0; (1,-1)],
            [2; 1; 1.0; (1, 0)],
            [2; 3; 1.0; (0, 0)],
            [2; 1; 1.0; (0, 0)],
            [2; 3; 1.0; (0,-1)],

            [3; 1; 1.0; (0, 0)],
            [3; 4; 1.0; (0, 0)],
            [3; 2; 1.0; (0, 0)],
            [3; 4; 1.0; (1, 0)],
            [3; 2; 1.0; (0, 1)],

            [4; 2; 1.0; (-1,1)],
            [4; 3; 1.0; (0, 0)],
            [4; 1; 1.0; (0, 1)],
            [4; 3; 1.0; (-1,0)],
            [4; 1; 1.0; (0, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_shastry_sutherland_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0],
            [0.5, 0.0],
            [0.5, 0.5],
            [0.0, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 3; "t1"; (0, 0)],
            [1; 2; "t2"; (-1,0)],
            [1; 4; "t3"; (0,-1)],
            [1; 2; "t4"; (0, 0)],
            [1; 4; "t5"; (0, 0)],

            [2; 4; "t1"; (1,-1)],
            [2; 1; "t2"; (1, 0)],
            [2; 3; "t3"; (0, 0)],
            [2; 1; "t4"; (0, 0)],
            [2; 3; "t5"; (0,-1)],

            [3; 1; "t1"; (0, 0)],
            [3; 4; "t2"; (0, 0)],
            [3; 2; "t3"; (0, 0)],
            [3; 4; "t4"; (1, 0)],
            [3; 2; "t5"; (0, 1)],

            [4; 2; "t1"; (-1,1)],
            [4; 3; "t2"; (0, 0)],
            [4; 1; "t3"; (0, 1)],
            [4; 3; "t4"; (-1,0)],
            [4; 1; "t5"; (0, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_shastry_sutherland_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellShastrySutherland




"""
    getUnitcellAdvancedSquare([version::Int64=1; save::Bool=false])

get the implementation of the *2D advanced square lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.

NOTE: This function redirects to `getUnitcellSquare(3)` and is only kept for
historical reasons.



# Versions

#### 1 - simple (DEFAULT)

Bravais lattice vectors are

    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]

4 sites per unitcell, located at

    r1 = [0.0, 0.0]
    r2 = [0.5, 0.0]
    r3 = [0.5, 0.5]
    r4 = [0.0, 0.5]

couplings are `"t1"` and `"t2"` which are in a staggered square pattern.



# Examples

```julia-repl
julia> unitcell = getUnitcellAdvancedSquare()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellAdvancedSquare(version::Int64=1; save::Bool=false)
    # redirect to the square implmentation with version 2
    if version == 1
        return getUnitcellSquare(3, save=save)
    end
end
export getUnitcellAdvancedSquare




"""
    getUnitcellSquareOctagon([version::Int64=1; save::Bool=false])

get the implementation of the *2D square-octagon lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 4 - simple (DEFAULT) & simple Kitaev

Bravais lattice vectors are

    a1 = [3*sqrt(3.0)/4., -3*sqrt(3.0)/4.]
    a2 = [3*sqrt(3.0)/4.,  3*sqrt(3.0)/4.]

4 sites per unitcell, located at

    r1 = [0.0, 0.0]
    r2 = [0.0, 1/sqrt(3.0)]
    r3 = [-1/sqrt(3.0), 0.0]
    r4 = [-1/sqrt(3.0), 1/sqrt(3.0)]

For version `1`, the coupling values are `1.0` for all bonds.

For version `4`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`



#### 2 & 5 - scaled shifted simple & scaled shifted simple Kitaev

Bravais lattice vectors are

    a1 = [1.0 + 1./sqrt(2.0), -(1.0 + 1./sqrt(2.0))]
    a2 = [1.0 + 1./sqrt(2.0),  (1.0 + 1./sqrt(2.0))]

4 sites per unitcell, located at

    r1 = [ 0.5, -0.5]
    r2 = [ 0.5,  0.5]
    r3 = [-0.5, -0.5]
    r4 = [-0.5,  0.5]

For version `2`, the coupling values are `1.0` for all bonds.

For version `5`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`




# Examples

```julia-repl
julia> unitcell = getUnitcellSquareOctagon()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellSquareOctagon(version::Int64=1; save::Bool=false)
    # SQUARE OCTAGON LATTICE
    if version == 1
        # the lattice vectors
        a1 = [3*sqrt(3.0)/4., -3*sqrt(3.0)/4.]
        a2 = [3*sqrt(3.0)/4., 3*sqrt(3.0)/4.]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0],
            [0.0, 1/sqrt(3.0)],
			[-1/sqrt(3.0), 0.0],
			[-1/sqrt(3.0), 1/sqrt(3.0)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0, 0)],
			[2; 4; 1.0; (0, 0)],
			[3; 4; 1.0; (0, 0)],
			[3; 1; 1.0; (0, 0)],
			[2; 1; 1.0; (0, 0)],
			[4; 2; 1.0; (0, 0)],
			[4; 3; 1.0; (0, 0)],
			[1; 3; 1.0; (0, 0)],

            [3; 2; 1.0; (0, -1)],
            [2; 3; 1.0; (0, 1)],
            [4; 1; 1.0; (-1, 0)],
            [1; 4; 1.0; (1, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_squareoctagon_unitcell.jld"
    elseif version == 2
        # the lattice vectors
        a1 = [1.0 + 1./sqrt(2.0), -(1.0 + 1./sqrt(2.0))]
        a2 = [1.0 + 1./sqrt(2.0),  (1.0 + 1./sqrt(2.0))]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [ 0.5, -0.5],
            [ 0.5,  0.5],
            [-0.5, -0.5],
            [-0.5,  0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0, 0)],
            [2; 4; 1.0; (0, 0)],
            [3; 4; 1.0; (0, 0)],
            [3; 1; 1.0; (0, 0)],
            [2; 1; 1.0; (0, 0)],
            [4; 2; 1.0; (0, 0)],
            [4; 3; 1.0; (0, 0)],
            [1; 3; 1.0; (0, 0)],
            [3; 2; 1.0; (0, -1)],
            [2; 3; 1.0; (0, 1)],
            [4; 1; 1.0; (-1, 0)],
            [1; 4; 1.0; (1, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_squareoctagon_alternative_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [3*sqrt(3.0)/4., -3*sqrt(3.0)/4.]
        a2 = [3*sqrt(3.0)/4., 3*sqrt(3.0)/4.]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0],
            [0.0, 1/sqrt(3.0)],
			[-1/sqrt(3.0), 0.0],
			[-1/sqrt(3.0), 1/sqrt(3.0)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; "tx"; (0,0)],
            [2; 4; "ty"; (0,0)],
            [3; 4; "tx"; (0,0)],
            [3; 1; "ty"; (0,0)],
            [2; 1; "tx"; (0,0)],
            [4; 2; "ty"; (0,0)],
            [4; 3; "tx"; (0,0)],
            [1; 3; "ty"; (0,0)],
            [3; 2; "tz"; (0,-1)],
            [2; 3; "tz"; (0,1)],
            [4; 1; "tz"; (-1,0)],
            [1; 4; "tz"; (1,0)],
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_squareoctagon_kitaev_unitcell.jld"
    elseif version == 5
        # the lattice vectors
        a1 = [1.0 + 1./sqrt(2.0), -(1.0 + 1./sqrt(2.0))]
        a2 = [1.0 + 1./sqrt(2.0),  (1.0 + 1./sqrt(2.0))]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [ 0.5, -0.5],
            [ 0.5,  0.5],
            [-0.5, -0.5],
            [-0.5,  0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; "tx"; (0,0)],
            [2; 4; "ty"; (0,0)],
            [3; 4; "tx"; (0,0)],
            [3; 1; "ty"; (0,0)],
            [2; 1; "tx"; (0,0)],
            [4; 2; "ty"; (0,0)],
            [4; 3; "tx"; (0,0)],
            [1; 3; "ty"; (0,0)],
            [3; 2; "tz"; (0,-1)],
            [2; 3; "tz"; (0,1)],
            [4; 1; "tz"; (-1,0)],
            [1; 4; "tz"; (1,0)],
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_squareoctagon_alternative_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellSquareOctagon



"""
    getUnitcellBCC2D([version::Int64=1; save::Bool=false])

get the implementation of the *2D bcc-like square lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.

NOTE: This function redirects to `getUnitcellSquare(2)` and is only kept for
historical reasons.



# Versions

#### 1 - simple (DEFAULT)

Bravais lattice vectors are

    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]

2 sites per unitcell, located at

    r1 = [0.0, 0.0]
    r2 = [0.5, 0.5]

All couplings have strength 1.0.


# Examples

```julia-repl
julia> unitcell = getUnitcellAdvancedSquare()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellBCC2D(version::Int64=1; save::Bool=false)

    # only return the redirected Function
    return getUnitcellSquare(2, save=save)

end
export getUnitcellBCC2D



"""
    getUnitcellFullyConnectedSquare([version::Int64=1; save::Bool=false, J1=1.0, J1X=0.5])

get the implementation of the *2D (fully connected) square lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 - simple (DEFAULT)

Bravais lattice vectors are

    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]

1 site per unitcell, located at

    r1 = [0.0, 0.0]

All couplings have strength `J1` and `J1X` with default values

    J1  = 1.0
    J1X = 0.5

but can be passed exactly as well.




# Examples

```julia-repl
julia> unitcell = getUnitcellFullyConnectedSquare()
LatticePhysics.Unitcell(...)

julia> unitcell = getUnitcellFullyConnectedSquare(J1=1.0, J1X=10.0)
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellFullyConnectedSquare(version::Int64=1; save::Bool=false, J1=1.0, J1X=0.5)
    if version == 1
        if J1==1.0 && J1X==0.5
            # the lattice vectors
            a1 = [1.0, 0.0]
            a2 = [0.0, 1.0]
            lattice_vectors = Array{Float64, 1}[]
            push!(lattice_vectors, a1)
            push!(lattice_vectors, a2)
            # Basis Definition
            basis = Array{Float64, 1}[
                [0.0, 0.0]
            ]
            # Connection Definition
            # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            connections = Array{Any, 1}[
                [1; 1; 1.0; (0, 1)],
                [1; 1; 1.0; (0, -1)],
                [1; 1; 1.0; (1, 0)],
                [1; 1; 1.0; (-1, 0)],
                [1; 1; 0.5; (1, 1)],
                [1; 1; 0.5; (-1, -1)],
                [1; 1; 0.5; (1, -1)],
                [1; 1; 0.5; (-1, 1)]
            ]
            # filename
            filename = "$(FOLDER_UNITCELLS)2d_fully_connected_square_unitcell.jld"
        else
            # the lattice vectors
            a1 = [1.0, 0.0]
            a2 = [0.0, 1.0]
            lattice_vectors = Array{Float64, 1}[]
            push!(lattice_vectors, a1)
            push!(lattice_vectors, a2)
            # Basis Definition
            basis = Array{Float64, 1}[
                [0.0, 0.0]
            ]
            # Connection Definition
            # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            connections = Array{Any, 1}[
                [1; 1; J1; (0, 1)],
                [1; 1; J1; (0, -1)],
                [1; 1; J1; (1, 0)],
                [1; 1; J1; (-1, 0)],
                [1; 1; J1X; (1, 1)],
                [1; 1; J1X; (-1, -1)],
                [1; 1; J1X; (1, -1)],
                [1; 1; J1X; (-1, 1)]
            ]
            # filename
            filename = "$(FOLDER_UNITCELLS)2d_fully_connected_square_$(J1)_$(J1X)_unitcell.jld"
        end
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellFullyConnectedSquare

















#----------------------------
#
#   2D TRIANGULAR BASIS
#
#----------------------------




"""
    getUnitcellTriangular([version::Int64=1; save::Bool=false])

get the implementation of the *2D triangular lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 - simple (DEFAULT)

Bravais lattice vectors are

    a1 = [sqrt(3.0)/2, -0.5]
    a2 = [sqrt(3.0)/2, +0.5]

1 site per unitcell, located at

    r1 = [0.0, 0.0]

All couplings have strength `1.0`.




#### 3 - anisotropic couplings

Bravais lattice vectors are

    a1 = [sqrt(3.0)/2, -0.5]
    a2 = [sqrt(3.0)/2, +0.5]

1 site per unitcell, located at

    r1 = [0.0, 0.0]

Couplings have strengths `"1"` and `"2"`.



# Examples

```julia-repl
julia> unitcell = getUnitcellTriangular()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellTriangular(version::Int64=1; save::Bool=false)
    # SIMPLE TRIANGULAR LATTICE
    if version == 1
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 1; 1.0; (0, 1)],
            [1; 1; 1.0; (0, -1)],
            [1; 1; 1.0; (1, 0)],
            [1; 1; 1.0; (-1, 0)],
            [1; 1; 1.0; (1, -1)],
            [1; 1; 1.0; (-1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_triangular_unitcell.jld"
    elseif version == 3
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 1; "1"; (0, 1)],
            [1; 1; "1"; (0, -1)],
            [1; 1; "1"; (1, 0)],
            [1; 1; "1"; (-1, 0)],
            [1; 1; "2"; (1, -1)],
            [1; 1; "2"; (-1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_triangular_aniso_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellTriangular




"""
    getUnitcellHoneycomb([version::Int64=1; save::Bool=false])

get the implementation of the *2D honeycomb lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 4 - simple ZZ (DEFAULT) & simple ZZ Kitaev

Bravais lattice vectors are

    a1 = [sqrt(3.0)/2, -0.5]
    a2 = [sqrt(3.0)/2, +0.5]

2 sites per unitcell, located at

    r1 = [0.0, 0.0]
    r2 = [1/sqrt(3.0), 0.0]

In a stripe, these definitions lead to a zigzag boundary.

For version `1`, the coupling values are `1.0` for all bonds.

For version `4`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`.





#### 2 & 5 - simple AC & simple AC Kitaev

Bravais lattice vectors are

    a1 = [sqrt(3.0)/2, -0.5]
    a2 = [sqrt(3.0),    0.0]

2 sites per unitcell, located at

    r1 = [0.0, 0.0]
    r2 = [1/sqrt(3.0), 0.0]

In a stripe, these definitions lead to an armchair boundary.

For version `2`, the coupling values are `1.0` for all bonds.

For version `5`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`.




#### 3 - anisotropic interaction ZZ

Bravais lattice vectors are

    a1 = [sqrt(3.0)/2, -0.5]
    a2 = [sqrt(3.0)/2, +0.5]

2 sites per unitcell, located at

    r1 = [0.0, 0.0]
    r2 = [1/sqrt(3.0), 0.0]

In a stripe, these definitions lead to a zigzag boundary.

Couplings have strengths `"1"` and `"2"`.





# Examples

```julia-repl
julia> unitcell = getUnitcellHoneycomb()
LatticePhysics.Unitcell(...)

julia> unitcell = getUnitcellHoneycomb(2)
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellHoneycomb(version::Int64=1; save::Bool=false)
    # SIMPLE HONEYCOMB LATTICE
    if version == 1
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0, 0)],
            [1; 2; 1.0; (-1, 0)],
            [1; 2; 1.0; (0, -1)],
            [2; 1; 1.0; (0, 0)],
            [2; 1; 1.0; (1, 0)],
            [2; 1; 1.0; (0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_honeycomb_ZZ_unitcell.jld"
    elseif version == 2
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0), 0.0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0, 0)],
            [1; 2; 1.0; (-1, 0)],
            [1; 2; 1.0; (1, -1)],
            [2; 1; 1.0; (0, 0)],
            [2; 1; 1.0; (1, 0)],
            [2; 1; 1.0; (-1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_honeycomb_AC_unitcell.jld"
    elseif version == 3
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; "1"; (0, 0)],
            [1; 2; "2"; (-1, 0)],
            [1; 2; "2"; (0, -1)],
            [2; 1; "1"; (0, 0)],
            [2; 1; "2"; (1, 0)],
            [2; 1; "2"; (0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_honeycomb_aniso_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; "tx"; (0, 0)],
            [1; 2; "ty"; (-1, 0)],
            [1; 2; "tz"; (0, -1)],
            [2; 1; "tx"; (0, 0)],
            [2; 1; "ty"; (1, 0)],
            [2; 1; "tz"; (0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_honeycomb_ZZ_kitaev_unitcell.jld"
    elseif version == 5
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0), 0.0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; "tx"; (0, 0)],
            [1; 2; "ty"; (-1, 0)],
            [1; 2; "tz"; (1, -1)],
            [2; 1; "tx"; (0, 0)],
            [2; 1; "ty"; (1, 0)],
            [2; 1; "tz"; (-1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_honeycomb_AC_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellHoneycomb


"""
    getUnitcellKagome([version::Int64=1; save::Bool=false])

get the implementation of the *2D Kagome lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 - simple (DEFAULT)

Bravais lattice vectors are

    a1 = [sqrt(3.0)/2, -0.5]
    a2 = [sqrt(3.0)/2, +0.5]

3 sites per unitcell, located at

    r1 = [0.0, 0.0]
    r2 = [sqrt(3.0)/4, -0.25]
    r3 = [sqrt(3.0)/4, +0.25]

All couplings have strength 1.0.





# Examples

```julia-repl
julia> unitcell = getUnitcellKagome()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellKagome(version::Int64=1; save::Bool=false)
    # SIMPLE KAGOME LATTICE
    if version == 1
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0],
            [sqrt(3.0)/4, -0.25],
            [sqrt(3.0)/4, +0.25]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0, 0)],
            [1; 3; 1.0; (0, 0)],
            [1; 2; 1.0; (-1, 0)],
            [1; 3; 1.0; (0, -1)],

            [2; 1; 1.0; (0, 0)],
            [2; 3; 1.0; (0, 0)],
            [2; 1; 1.0; (1, 0)],
            [2; 3; 1.0; (1, -1)],

            [3; 1; 1.0; (0, 0)],
            [3; 2; 1.0; (0, 0)],
            [3; 1; 1.0; (0, 1)],
            [3; 2; 1.0; (-1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_kagome_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellKagome



"""
    getUnitcellKagomeMinus([version::Int64=1; save::Bool=false])

get the implementation of the *2D "Kagome^-" lattice* unitcell
(similar looking to the Kagome but with additional bonds between triangles). The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 4 - simple (DEFAULT) & simple Kitaev

Bravais lattice vectors are

    a1 = [sqrt(3.0)/2, -0.5]
    a2 = [sqrt(3.0)/2, +0.5]

6 sites per unitcell, located at

    r1  =  b1  +  1/4 * (b2 - b1)
    r2  =  b1  +  1/4 * (b2 - a2)
    r3  =  b1  +  1/4 * (b2 - a1)
    r4  =  b2  +  1/4 * (b1 - b2)
    r5  =  b2  +  1/4 * (a2 - b2)
    r6  =  b2  +  1/4 * (a1 - b2)

with the additional definition of

    b1 = [0.0, 0.0]
    b2 = [1/sqrt(3.0), 0.0]


For version `1`, all couplings have strength 1.0.

For version `4`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`
as well as `"txp"`, `"typ"` and `"tzp"` for the second triangle.





# Examples

```julia-repl
julia> unitcell = getUnitcellKagomeMinus()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellKagomeMinus(version::Int64=1; save::Bool=false)
    # SIMPLE KAGOME MINUS LATTICE
    if version == 1
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        b1 = [0.0, 0.0]
        b2 = [1/sqrt(3.0), 0.0]
        basis = Array{Float64, 1}[
            b1 .+ 1/4 * (b2 .- b1),
            b1 .+ 1/4 * (b2 .- a2),
            b1 .+ 1/4 * (b2 .- a1),
            b2 .+ 1/4 * (b1 .- b2),
            b2 .+ 1/4 * (a2 .- b2),
            b2 .+ 1/4 * (a1 .- b2)
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0,0)],
            [2; 1; 1.0; (0,0)],
            [1; 3; 1.0; (0,0)],
            [3; 1; 1.0; (0,0)],
            [2; 3; 1.0; (0,0)],
            [3; 2; 1.0; (0,0)],
            [4; 5; 1.0; (0,0)],
            [5; 4; 1.0; (0,0)],
            [4; 6; 1.0; (0,0)],
            [6; 4; 1.0; (0,0)],
            [5; 6; 1.0; (0,0)],
            [6; 5; 1.0; (0,0)],
            [1; 4; 1.0; (0,0)],
            [4; 1; 1.0; (0,0)],
            [2; 5; 1.0; (0,-1)],
            [5; 2; 1.0; (0,1)],
            [3; 6; 1.0; (-1,0)],
            [6; 3; 1.0; (1,0)],
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_kagome_minus_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        b1 = [0.0, 0.0]
        b2 = [1/sqrt(3.0), 0.0]
        basis = Array{Float64, 1}[
            b1 .+ 1/4 * (b2 .- b1),
            b1 .+ 1/4 * (b2 .- a2),
            b1 .+ 1/4 * (b2 .- a1),
            b2 .+ 1/4 * (b1 .- b2),
            b2 .+ 1/4 * (a2 .- b2),
            b2 .+ 1/4 * (a1 .- b2)
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; "tzp"; (0,0)],
            [2; 1; "tzp"; (0,0)],
            [1; 3; "txp"; (0,0)],
            [3; 1; "txp"; (0,0)],
            [2; 3; "typ"; (0,0)],
            [3; 2; "typ"; (0,0)],
            [4; 5; "tzp"; (0,0)],
            [5; 4; "tzp"; (0,0)],
            [4; 6; "txp"; (0,0)],
            [6; 4; "txp"; (0,0)],
            [5; 6; "typ"; (0,0)],
            [6; 5; "typ"; (0,0)],
            [1; 4; "ty"; (0,0)],
            [4; 1; "ty"; (0,0)],
            [2; 5; "tx"; (0,-1)],
            [5; 2; "tx"; (0,1)],
            [3; 6; "tz"; (-1,0)],
            [6; 3; "tz"; (1,0)],
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_kagome_minus_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellKagomeMinus




"""
    getUnitcellHoneycombXXX([version::Int64=1; save::Bool=false])

get the implementation of the *2D honeycomb-XXX lattice* unitcell
(similar looking to the honeycomb but with 3 additional sites on every bond). The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 2 - simple (DEFAULT) & fine-tuned simple

Bravais lattice vectors are

    a1 = [sqrt(3.0)/2, -0.5]
    a2 = [sqrt(3.0)/2, +0.5]

11 sites per unitcell, located at

    r1  = [0.0,        0.0]
    r2  = [0.57735,    0.0]
    r3  = [0.288675,   0.0]
    r4  = [-0.144338,  0.25]
    r5  = [-0.144338, -0.25]
    r6  = [0.144338,   0.0]
    r7  = [0.433013,   0.0]
    r8  = [-0.0721688, 0.125]
    r9  = [-0.216506,  0.375]
    r10 = [-0.0721688,-0.125]
    r11 = [-0.216506, -0.375]

For version `1`, all couplings have strength 1.0.

For version `2`, the couplings are fine-tuned to be a perfect squareroot
as `sqrt(sqrt(3.0/2.0))` and `sqrt(sqrt(2.0/3.0))`.





# Examples

```julia-repl
julia> unitcell = getUnitcellHoneycombXXX()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellHoneycombXXX(version::Int64=1; save::Bool=false)
    # distinguish by version
    if version == 1
       # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0,0.0],
            [0.57735,0.0],
            [0.288675,0.0],
            [-0.144338,0.25],
            [-0.144338,-0.25],
            [0.144338,0.0],
            [0.433013,0.0],
            [-0.0721688,0.125],
            [-0.216506,0.375],
            [-0.0721688,-0.125],
            [-0.216506,-0.375]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 6; 1.0; (0,0)],
            [6; 1; 1.0; (0,0)],
            [6; 3; 1.0; (0,0)],
            [3; 6; 1.0; (0,0)],
            [3; 7; 1.0; (0,0)],
            [7; 3; 1.0; (0,0)],
            [7; 2; 1.0; (0,0)],
            [2; 7; 1.0; (0,0)],
            [1; 8; 1.0; (0,0)],
            [8; 1; 1.0; (0,0)],
            [8; 4; 1.0; (0,0)],
            [4; 8; 1.0; (0,0)],
            [4; 9; 1.0; (0,0)],
            [9; 4; 1.0; (0,0)],
            [9; 2; 1.0; (-1,0)],
            [2; 9; 1.0; (1,0)],
            [1; 10; 1.0; (0,0)],
            [10; 1; 1.0; (0,0)],
            [10; 5; 1.0; (0,0)],
            [5; 10; 1.0; (0,0)],
            [5; 11; 1.0; (0,0)],
            [11; 5; 1.0; (0,0)],
            [11; 2; 1.0; (0,-1)],
            [2; 11; 1.0; (0,1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_honeycomb-XXX_unitcell.jld"
    elseif version == 2
       # the lattice vectors
        a1 = [sqrt(3.0)/2, -0.5]
        a2 = [sqrt(3.0)/2, +0.5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0,0.0],
            [0.57735,0.0],
            [0.288675,0.0],
            [-0.144338,0.25],
            [-0.144338,-0.25],
            [0.144338,0.0],
            [0.433013,0.0],
            [-0.0721688,0.125],
            [-0.216506,0.375],
            [-0.0721688,-0.125],
            [-0.216506,-0.375]
        ]
        # Connection Definition
        a = sqrt(sqrt(2.0/3.0))
        b = sqrt(sqrt(3.0/2.0))
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 6; a; (0,0)],
            [6; 1; a; (0,0)],
            [6; 3; b; (0,0)],
            [3; 6; b; (0,0)],
            [3; 7; b; (0,0)],
            [7; 3; b; (0,0)],
            [7; 2; a; (0,0)],
            [2; 7; a; (0,0)],
            [1; 8; a; (0,0)],
            [8; 1; a; (0,0)],
            [8; 4; b; (0,0)],
            [4; 8; b; (0,0)],
            [4; 9; b; (0,0)],
            [9; 4; b; (0,0)],
            [9; 2; a; (-1,0)],
            [2; 9; a; (1,0)],
            [1; 10; a; (0,0)],
            [10; 1; a; (0,0)],
            [10; 5; b; (0,0)],
            [5; 10; b; (0,0)],
            [5; 11; b; (0,0)],
            [11; 5; b; (0,0)],
            [11; 2; a; (0,-1)],
            [2; 11; a; (0,1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)2d_honeycomb-XXX_a_b_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellHoneycombXXX

































#-----------------------------------------------------------------------------------------------------------------------------
#
#   Individual FUNCTIONS FOR 3D UNITCELLS
#
#   They all generate the special filenames for the unitcells and return the unitcells as objects
#   Functions can be given a version integer to distinguish several implementations of the same lattice
#   Functions can be specified to already save the unitcell to file
#
#-----------------------------------------------------------------------------------------------------------------------------




#----------------------------
#
#   3D CUBIC (FCC) BASIS
#
#----------------------------


"""
    getUnitcellDiamond([version::Int64=1; save::Bool=false])

get the implementation of the *3D diamond lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 2 & 3 - simple (DEFAULT) & simple with NNN & simple with NNN and anisotropy

Bravais lattice vectors are

    a1 = 0.5 .* [0, 1, 1]
    a2 = 0.5 .* [1, 0, 1]
    a3 = 0.5 .* [1, 1, 0]

2 sites per unitcell, located at

    r1 = [0.0,  0.0,  0.0]
    r2 = [0.25, 0.25, 0.25]

For version `1`, all couplings have strength `1.0`.

For version `2`, the nearest neighbor couplings have strength `"J1"` and the lattice
includes next-nearest neighbor couplings with strength `"J2"`.

For version `3`, the next-nearest neighbor couplings are split between those with
strength `"J21"` and those with strength `"J22"`.





# Examples

```julia-repl
julia> unitcell = getUnitcellDiamond()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellDiamond(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = 0.5 .* [0, 1, 1]
        a2 = 0.5 .* [1, 0, 1]
        a3 = 0.5 .* [1, 1, 0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0, 0.0],
            [0.25, 0.25, 0.25]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0, 0, 0)],
            [1; 2; 1.0; (-1, 0, 0)],
            [1; 2; 1.0; (0, -1, 0)],
            [1; 2; 1.0; (0, 0, -1)],

            [2; 1; 1.0; (0, 0, 0)],
            [2; 1; 1.0; (1, 0, 0)],
            [2; 1; 1.0; (0, 1, 0)],
            [2; 1; 1.0; (0, 0, 1)],
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_diamond_unitcell.jld"
    elseif version == 2
        # the lattice vectors
        a1 = 0.5 .* [0, 1, 1]
        a2 = 0.5 .* [1, 0, 1]
        a3 = 0.5 .* [1, 1, 0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0, 0.0],
            [0.25, 0.25, 0.25]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; "J1"; (0, 0, 0)],
            [1; 2; "J1"; (-1, 0, 0)],
            [1; 2; "J1"; (0, -1, 0)],
            [1; 2; "J1"; (0, 0, -1)],

            [2; 1; "J1"; (0, 0, 0)],
            [2; 1; "J1"; (1, 0, 0)],
            [2; 1; "J1"; (0, 1, 0)],
            [2; 1; "J1"; (0, 0, 1)],


            [1; 1; "J2"; (1, 0, 0)],
            [1; 1; "J2"; (-1, 0, 0)],
            [1; 1; "J2"; (0, 1, 0)],
            [1; 1; "J2"; (0, -1, 0)],
            [1; 1; "J2"; (0, 0, 1)],
            [1; 1; "J2"; (0, 0, -1)],

            [1; 1; "J2"; (1, -1, 0)],
            [1; 1; "J2"; (-1, 1, 0)],
            [1; 1; "J2"; (0, 1, -1)],
            [1; 1; "J2"; (0, -1, 1)],
            [1; 1; "J2"; (-1, 0, 1)],
            [1; 1; "J2"; (1, 0, -1)],


            [2; 2; "J2"; (1, 0, 0)],
            [2; 2; "J2"; (-1, 0, 0)],
            [2; 2; "J2"; (0, 1, 0)],
            [2; 2; "J2"; (0, -1, 0)],
            [2; 2; "J2"; (0, 0, 1)],
            [2; 2; "J2"; (0, 0, -1)],

            [2; 2; "J2"; (1, -1, 0)],
            [2; 2; "J2"; (-1, 1, 0)],
            [2; 2; "J2"; (0, 1, -1)],
            [2; 2; "J2"; (0, -1, 1)],
            [2; 2; "J2"; (-1, 0, 1)],
            [2; 2; "J2"; (1, 0, -1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_diamond_NN_NNN_unitcell.jld"
    elseif version == 3
        # the lattice vectors
        a1 = 0.5 .* [0, 1, 1]
        a2 = 0.5 .* [1, 0, 1]
        a3 = 0.5 .* [1, 1, 0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0, 0.0],
            [0.25, 0.25, 0.25]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; "J1"; (0, 0, 0)],
            [1; 2; "J1"; (-1, 0, 0)],
            [1; 2; "J1"; (0, -1, 0)],
            [1; 2; "J1"; (0, 0, -1)],

            [2; 1; "J1"; (0, 0, 0)],
            [2; 1; "J1"; (1, 0, 0)],
            [2; 1; "J1"; (0, 1, 0)],
            [2; 1; "J1"; (0, 0, 1)],


            [1; 1; "J21"; (1, 0, 0)],
            [1; 1; "J21"; (-1, 0, 0)],
            [1; 1; "J21"; (0, 1, 0)],
            [1; 1; "J21"; (0, -1, 0)],
            [1; 1; "J22"; (0, 0, 1)],
            [1; 1; "J22"; (0, 0, -1)],

            [1; 1; "J22"; (1, -1, 0)],
            [1; 1; "J22"; (-1, 1, 0)],
            [1; 1; "J21"; (0, 1, -1)],
            [1; 1; "J21"; (0, -1, 1)],
            [1; 1; "J21"; (-1, 0, 1)],
            [1; 1; "J21"; (1, 0, -1)],


            [2; 2; "J21"; (1, 0, 0)],
            [2; 2; "J21"; (-1, 0, 0)],
            [2; 2; "J21"; (0, 1, 0)],
            [2; 2; "J21"; (0, -1, 0)],
            [2; 2; "J22"; (0, 0, 1)],
            [2; 2; "J22"; (0, 0, -1)],

            [2; 2; "J22"; (1, -1, 0)],
            [2; 2; "J22"; (-1, 1, 0)],
            [2; 2; "J21"; (0, 1, -1)],
            [2; 2; "J21"; (0, -1, 1)],
            [2; 2; "J21"; (-1, 0, 1)],
            [2; 2; "J21"; (1, 0, -1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_diamond_NN_aniso_NNN_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellDiamond


"""
    getUnitcellFCC([version::Int64=1; save::Bool=false])

get the implementation of the *3D fcc lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 4 - simple (DEFAULT) & simple Kitaev

Bravais lattice vectors are

    a1 = 0.5 .* [0, 1, 1]
    a2 = 0.5 .* [1, 0, 1]
    a3 = 0.5 .* [1, 1, 0]

1 site per unitcell, located at

    r1 = [0.0,  0.0,  0.0]

Connections are only between nearest neighbors.

For version `1`, all couplings have strength `1.0`.

For version `4`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`.





#### 2 & 5 - NNN & Kitaev + NNN

Bravais lattice vectors are

    a1 = 0.5 .* [0, 1, 1]
    a2 = 0.5 .* [1, 0, 1]
    a3 = 0.5 .* [1, 1, 0]

1 site per unitcell, located at

    r1 = [0.0,  0.0,  0.0]

Connections are between nearest neighbors and next-nearest neighbors (which are chosen by real-space distance).

For version `1`, nearest neighbors have strength `"J1"` and next-nearest neighbors have strength `"J2"`.

For version `4`, the nearest-neighbor couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`
whereas next-nearest neighbors have strength `"J2"`.




# Examples

```julia-repl
julia> unitcell = getUnitcellFCC()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellFCC(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = 0.5 .* [0, 1, 1]
        a2 = 0.5 .* [1, 0, 1]
        a3 = 0.5 .* [1, 1, 0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 1; 1.0; (1, 0, 0)],
            [1; 1; 1.0; (-1, 0, 0)],
            [1; 1; 1.0; (0, 1, -1)],
            [1; 1; 1.0; (0, -1, 1)],
            [1; 1; 1.0; (0, 1, 0)],
            [1; 1; 1.0; (0, -1, 0)],
            [1; 1; 1.0; (-1, 0, 1)],
            [1; 1; 1.0; (1, 0, -1)],
            [1; 1; 1.0; (0, 0, 1)],
            [1; 1; 1.0; (0, 0, -1)],
            [1; 1; 1.0; (-1, 1, 0)],
            [1; 1; 1.0; (1, -1, 0)],
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_fcc_unitcell.jld"
    elseif version == 2
        # the lattice vectors
        a1 = 0.5 .* [0, 1, 1]
        a2 = 0.5 .* [1, 0, 1]
        a3 = 0.5 .* [1, 1, 0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 1; "J1"; (1, 0, 0)],
            [1; 1; "J1"; (-1, 0, 0)],
            [1; 1; "J1"; (0, 1, -1)],
            [1; 1; "J1"; (0, -1, 1)],
            [1; 1; "J1"; (0, 1, 0)],
            [1; 1; "J1"; (0, -1, 0)],
            [1; 1; "J1"; (-1, 0, 1)],
            [1; 1; "J1"; (1, 0, -1)],
            [1; 1; "J1"; (0, 0, 1)],
            [1; 1; "J1"; (0, 0, -1)],
            [1; 1; "J1"; (-1, 1, 0)],
            [1; 1; "J1"; (1, -1, 0)],
            [1; 1; "J2"; (1, 1, -1)],
            [1; 1; "J2"; (1, -1, 1)],
            [1; 1; "J2"; (-1, 1, 1)],
            [1; 1; "J2"; (-1, -1, 1)],
            [1; 1; "J2"; (1, -1, -1)],
            [1; 1; "J2"; (-1, 1, -1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_fcc_NN_NNN_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = 0.5 .* [0, 1, 1]
        a2 = 0.5 .* [1, 0, 1]
        a3 = 0.5 .* [1, 1, 0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 1; "tx"; (1, 0, 0)],
            [1; 1; "tx"; (-1, 0, 0)],
            [1; 1; "tx"; (0, 1, -1)],
            [1; 1; "tx"; (0, -1, 1)],
            [1; 1; "ty"; (0, 1, 0)],
            [1; 1; "ty"; (0, -1, 0)],
            [1; 1; "ty"; (-1, 0, 1)],
            [1; 1; "ty"; (1, 0, -1)],
            [1; 1; "tz"; (0, 0, 1)],
            [1; 1; "tz"; (0, 0, -1)],
            [1; 1; "tz"; (-1, 1, 0)],
            [1; 1; "tz"; (1, -1, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_fcc_kitaev_unitcell.jld"
    elseif version == 5
        # the lattice vectors
        a1 = 0.5 .* [0, 1, 1]
        a2 = 0.5 .* [1, 0, 1]
        a3 = 0.5 .* [1, 1, 0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 1; "tx"; (1, 0, 0)],
            [1; 1; "tx"; (-1, 0, 0)],
            [1; 1; "tx"; (0, 1, -1)],
            [1; 1; "tx"; (0, -1, 1)],
            [1; 1; "ty"; (0, 1, 0)],
            [1; 1; "ty"; (0, -1, 0)],
            [1; 1; "ty"; (-1, 0, 1)],
            [1; 1; "ty"; (1, 0, -1)],
            [1; 1; "tz"; (0, 0, 1)],
            [1; 1; "tz"; (0, 0, -1)],
            [1; 1; "tz"; (-1, 1, 0)],
            [1; 1; "tz"; (1, -1, 0)],
            [1; 1; "J2"; (1, 1, -1)],
            [1; 1; "J2"; (1, -1, 1)],
            [1; 1; "J2"; (-1, 1, 1)],
            [1; 1; "J2"; (-1, -1, 1)],
            [1; 1; "J2"; (1, -1, -1)],
            [1; 1; "J2"; (-1, 1, -1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_fcc_J1_K1_J2_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellFCC


"""
    getUnitcellBCC([version::Int64=1; save::Bool=false])

get the implementation of the *3D body centered cubic (BCC) lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 - simple (DEFAULT)

Bravais lattice vectors are

    a1 = [1, 0, 0]
    a2 = [0, 1, 0]
    a3 = [0, 0, 1]

2 sites per unitcell, located at

    r1 = [0.0, 0.0, 0.0]
    r2 = [0.5, 0.5, 0.5]

All couplings have strength `1.0`.




# Examples

```julia-repl
julia> unitcell = getUnitcellBCC()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellBCC(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = [1, 0, 0]
        a2 = [0, 1, 0]
        a3 = [0, 0, 1]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; ( 0,  0,  0)],
            [1; 2; 1.0; (-1,  0,  0)],
            [1; 2; 1.0; ( 0, -1,  0)],
            [1; 2; 1.0; (-1, -1,  0)],
            [1; 2; 1.0; ( 0,  0, -1)],
            [1; 2; 1.0; (-1,  0, -1)],
            [1; 2; 1.0; ( 0, -1, -1)],
            [1; 2; 1.0; (-1, -1, -1)],

            [2; 1; 1.0; (0, 0, 0)],
            [2; 1; 1.0; (1, 0, 0)],
            [2; 1; 1.0; (0, 1, 0)],
            [2; 1; 1.0; (1, 1, 0)],
            [2; 1; 1.0; (0, 0, 1)],
            [2; 1; 1.0; (1, 0, 1)],
            [2; 1; 1.0; (0, 1, 1)],
            [2; 1; 1.0; (1, 1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_bcc_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellBCC




"""
    getUnitcellPyrochlore([version::Int64=1; save::Bool=false])

get the implementation of the *3D pyrochlore lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 - simple (DEFAULT)

Bravais lattice vectors are

    a1 = 0.5 .* [0, 1, 1]
    a2 = 0.5 .* [1, 0, 1]
    a3 = 0.5 .* [1, 1, 0]

4 sites per unitcell, located at

    r1 = [0.0,  0.0,  0.0]
    r2 = [0.0,  0.25, 0.25]
    r3 = [0.25, 0.0,  0.25]
    r4 = [0.25, 0.25, 0.0]

All couplings have strength `1.0`.





# Examples

```julia-repl
julia> unitcell = getUnitcellPyrochlore()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcellPyrochlore(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = [0, 0.5, 0.5]
        a2 = [0.5, 0, 0.5]
        a3 = [0.5, 0.5, 0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0., 0., 0.],
            [0., 0.25, 0.25],
            [0.25, 0., 0.25],
            [0.25, 0.25, 0.]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0, 0, 0)],
            [1; 3; 1.0; (0, 0, 0)],
            [1; 4; 1.0; (0, 0, 0)],
            [2; 1; 1.0; (0, 0, 0)],
            [2; 3; 1.0; (0, 0, 0)],
            [2; 4; 1.0; (0, 0, 0)],
            [3; 1; 1.0; (0, 0, 0)],
            [3; 2; 1.0; (0, 0, 0)],
            [3; 4; 1.0; (0, 0, 0)],
            [4; 1; 1.0; (0, 0, 0)],
            [4; 2; 1.0; (0, 0, 0)],
            [4; 3; 1.0; (0, 0, 0)],

            [1; 4; 1.0; (0, 0, -1)],
            [4; 1; 1.0; (0, 0, 1)],
            [1; 2; 1.0; (-1, 0, 0)],
            [2; 1; 1.0; (1, 0, 0)],
            [1; 3; 1.0; (0, -1, 0)],
            [3; 1; 1.0; (0, 1, 0)],

            [2; 3; 1.0; (1, -1, 0)],
            [3; 2; 1.0; (-1, 1, 0)],
            [2; 4; 1.0; (1, 0, -1)],
            [4; 2; 1.0; (-1, 0, 1)],

            [3; 4; 1.0; (0, 1, -1)],
            [4; 3; 1.0; (0, -1, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_pyrochlore_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcellPyrochlore



















#-----------------------------------------
#
#   (X,3)y FAMILY OF LATTICES (in 3D)
#
#-----------------------------------------





"""
    getUnitcell_8_3_a([version::Int64=1; save::Bool=false])

get the implementation of the *3D (8,3)a lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 4 - simple (DEFAULT) & simple Kitaev

Bravais lattice vectors are

    a1 = [ 1.0,         0.0,            0.0]
    a2 = [-0.5,  sqrt(3)/2.,            0.0]
    a3 = [ 0.0,         0.0, (3*sqrt(2))/5.]

6 sites per unitcell, located at

    r1 = [ 0.5,     sqrt(3)/10.,            0.0]
    r2 = [3/5.,      sqrt(3)/5., (2*sqrt(2))/5.]
    r3 = [ 0.1, (3*sqrt(3))/10.,     sqrt(2)/5.]
    r4 = [ 0.4,      sqrt(3)/5.,     sqrt(2)/5.]
    r5 = [ 0.0,  (2*sqrt(3))/5.,            0.0]
    r6 = [-0.1, (3*sqrt(3))/10., (2*sqrt(2))/5.]

For version `1`, all couplings have strength 1.0.

For version `4`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`.





# Examples

```julia-repl
julia> unitcell = getUnitcell_8_3_a()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcell_8_3_a(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = [1.0, 0.0, 0.0]
        a2 = [-0.5, sqrt(3)/2., 0.0]
        a3 = [0.0, 0.0, (3*sqrt(2))/5.]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.5, sqrt(3)/10., 0.0],
            [3/5., sqrt(3)/5., (2*sqrt(2))/5.],
            [0.1, (3*sqrt(3))/10., sqrt(2)/5.],
            [0.4, sqrt(3)/5., sqrt(2)/5.],
            [0.0, (2*sqrt(3))/5., 0.0],
            [-0.1, (3*sqrt(3))/10., (2*sqrt(2))/5.]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 4; 1.0; (0, 0, 0)],
            [4; 2; 1.0; (0, 0, 0)], # zz
            [4; 3; 1.0; (0, 0, 0)],
            [5; 3; 1.0; (0, 0, 0)],
            [3; 6; 1.0; (0, 0, 0)], # zz

            [4; 1; 1.0; (0, 0, 0)],
            [2; 4; 1.0; (0, 0, 0)], # zz
            [3; 4; 1.0; (0, 0, 0)],
            [3; 5; 1.0; (0, 0, 0)],
            [6; 3; 1.0; (0, 0, 0)], # zz

            [5; 1; 1.0; (0, 1, 0)], # zz
            [1; 5; 1.0; (0, -1, 0)], # zz

            [2; 6; 1.0; (1, 0, 0)],
            [6; 2; 1.0; (-1, 0, 0)],

            [1; 2; 1.0; (0, 0, -1)],
            [2; 1; 1.0; (0, 0, 1)],

            [5; 6; 1.0; (0, 0, -1)],
            [6; 5; 1.0; (0, 0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_8_3_a_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [1.0, 0.0, 0.0]
        a2 = [-0.5, sqrt(3)/2., 0.0]
        a3 = [0.0, 0.0, (3*sqrt(2))/5.]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.5, sqrt(3)/10., 0.0],
            [3/5., sqrt(3)/5., (2*sqrt(2))/5.],
            [0.1, (3*sqrt(3))/10., sqrt(2)/5.],
            [0.4, sqrt(3)/5., sqrt(2)/5.],
            [0.0, (2*sqrt(3))/5., 0.0],
            [-0.1, (3*sqrt(3))/10., (2*sqrt(2))/5.]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 4; "ty"; (0, 0, 0)],
            [4; 2; "tz"; (0, 0, 0)], # zz
            [4; 3; "tx"; (0, 0, 0)],
            [5; 3; "ty"; (0, 0, 0)],
            [3; 6; "tz"; (0, 0, 0)], # zz

            [4; 1; "ty"; (0, 0, 0)],
            [2; 4; "tz"; (0, 0, 0)], # zz
            [3; 4; "tx"; (0, 0, 0)],
            [3; 5; "ty"; (0, 0, 0)],
            [6; 3; "tz"; (0, 0, 0)], # zz

            [5; 1; "tz"; (0, 1, 0)], # zz
            [1; 5; "tz"; (0, -1, 0)], # zz

            [2; 6; "ty"; (1, 0, 0)],
            [6; 2; "ty"; (-1, 0, 0)],

            [1; 2; "tx"; (0, 0, -1)],
            [2; 1; "tx"; (0, 0, 1)],

            [5; 6; "tx"; (0, 0, -1)],
            [6; 5; "tx"; (0, 0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_8_3_a_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcell_8_3_a




"""
    getUnitcell_8_3_b([version::Int64=1; save::Bool=false])

get the implementation of the *3D (8,3)b lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 4 - simple (DEFAULT) & simple Kitaev

Bravais lattice vectors are

    a1 = [1/2., 1/(2*sqrt(3)),     sqrt(2)/(5*sqrt(3))]
    a2 = [   0,     1/sqrt(3), (2*sqrt(2))/(5*sqrt(3))]
    a3 = [   0,             0,               sqrt(6)/5]

6 sites per unitcell, located at

    r1 = [1/10.,   1/(2*sqrt(3)),     sqrt(2)/(5*sqrt(3))]
    r2 = [ 1/5.,       sqrt(3)/5,               sqrt(6)/5]
    r3 = [3/10., 11/(10*sqrt(3)), (4*sqrt(2))/(5*sqrt(3))]
    r4 = [ 1/5.,   2/(5*sqrt(3)), (2*sqrt(2))/(5*sqrt(3))]
    r5 = [3/10., (3*sqrt(3))/10.,               sqrt(6)/5]
    r6 = [ 2/5.,       1/sqrt(3),         sqrt(2)/sqrt(3)]

For version `1`, all couplings have strength 1.0.

For version `4`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`.





# Examples

```julia-repl
julia> unitcell = getUnitcell_8_3_b()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcell_8_3_b(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = [1/2., 1/(2*sqrt(3)), sqrt(2)/(5*sqrt(3))]
        a2 = [0, 1/sqrt(3), (2*sqrt(2))/(5*sqrt(3))]
        a3 = [0, 0, sqrt(6)/5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [1/10., 1/(2*sqrt(3)), sqrt(2)/(5*sqrt(3))],
            [1/5., sqrt(3)/5, sqrt(6)/5],
            [3/10., 11/(10*sqrt(3)), (4*sqrt(2))/(5*sqrt(3))],
            [1/5., 2/(5*sqrt(3)), (2*sqrt(2))/(5*sqrt(3))],
            [3/10., (3*sqrt(3))/10., sqrt(6)/5],
            [2/5., 1/sqrt(3), sqrt(2)/sqrt(3)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 4; 1.0; (0, 0, 0)], # zz
            [4; 2; 1.0; (0, 0, 0)],
            [2; 5; 1.0; (0, 0, 0)], # zz
            [5; 3; 1.0; (0, 0, 0)],
            [3; 6; 1.0; (0, 0, 0)], # zz

            [4; 1; 1.0; (0, 0, 0)], # zz
            [2; 4; 1.0; (0, 0, 0)],
            [5; 2; 1.0; (0, 0, 0)], # zz
            [3; 5; 1.0; (0, 0, 0)],
            [6; 3; 1.0; (0, 0, 0)], # zz

            [6; 1; 1.0; (1, 0, 1)],
            [1; 6; 1.0; (-1, 0, -1)],

            [4; 3; 1.0; (0, -1, 0)],
            [3; 4; 1.0; (0, 1, 0)],

            [1; 2; 1.0; (0, 0, -1)],
            [2; 1; 1.0; (0, 0, 1)],

            [5; 6; 1.0; (0, 0, -1)],
            [6; 5; 1.0; (0, 0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_8_3_b_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [1/2., 1/(2*sqrt(3)), sqrt(2)/(5*sqrt(3))]
        a2 = [0, 1/sqrt(3), (2*sqrt(2))/(5*sqrt(3))]
        a3 = [0, 0, sqrt(6)/5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [1/10., 1/(2*sqrt(3)), sqrt(2)/(5*sqrt(3))],
            [1/5., sqrt(3)/5, sqrt(6)/5],
            [3/10., 11/(10*sqrt(3)), (4*sqrt(2))/(5*sqrt(3))],
            [1/5., 2/(5*sqrt(3)), (2*sqrt(2))/(5*sqrt(3))],
            [3/10., (3*sqrt(3))/10., sqrt(6)/5],
            [2/5., 1/sqrt(3), sqrt(2)/sqrt(3)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 4; "tz"; (0, 0, 0)], # zz
            [4; 2; "ty"; (0, 0, 0)],
            [2; 5; "tz"; (0, 0, 0)], # zz
            [5; 3; "ty"; (0, 0, 0)],
            [3; 6; "tz"; (0, 0, 0)], # zz

            [4; 1; "tz"; (0, 0, 0)], # zz
            [2; 4; "ty"; (0, 0, 0)],
            [5; 2; "tz"; (0, 0, 0)], # zz
            [3; 5; "ty"; (0, 0, 0)],
            [6; 3; "tz"; (0, 0, 0)], # zz

            [6; 1; "ty"; (1, 0, 1)],
            [1; 6; "ty"; (-1, 0, -1)],

            [4; 3; "tx"; (0, -1, 0)],
            [3; 4; "tx"; (0, 1, 0)],

            [1; 2; "tx"; (0, 0, -1)],
            [2; 1; "tx"; (0, 0, 1)],

            [5; 6; "tx"; (0, 0, -1)],
            [6; 5; "tx"; (0, 0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_8_3_b_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcell_8_3_b



"""
    getUnitcell_8_3_c([version::Int64=1; save::Bool=false])

get the implementation of the *3D (8,3)c lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 4 - simple (DEFAULT) & simple Kitaev

Bravais lattice vectors are

    a1 = [   1.,         0.,   0.]
    a2 = [-1/2., sqrt(3)/2.,   0.]
    a3 = [   0.,         0., 2/5.]

8 sites per unitcell, located at

    r1 = [ -1/5.,  4/(5*sqrt(3)), 1/10.]
    r2 = [    0.,  7/(5*sqrt(3)), 1/10.]
    r3 = [  1/5.,  4/(5*sqrt(3)), 1/10.]
    r4 = [  1/2.,  1/(2*sqrt(3)), 3/10.]
    r5 = [    0.,      1/sqrt(3), 1/10.]
    r6 = [ 3/10., 7/(10*sqrt(3)), 3/10.]
    r7 = [  1/2., 1/(10*sqrt(3)), 3/10.]
    r8 = [ 7/10., 7/(10*sqrt(3)), 3/10.]

For version `1`, all couplings have strength 1.0.

For version `4`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`.





# Examples

```julia-repl
julia> unitcell = getUnitcell_8_3_c()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcell_8_3_c(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = [1., 0., 0.]
        a2 = [-1/2., sqrt(3)/2., 0.]
        a3 = [0., 0., 2/5.]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [-1/5., 4/(5*sqrt(3)), 1/10.],  #1
            [ 0., 7/(5*sqrt(3)), 1/10.],    #2
            [ 1/5., 4/(5*sqrt(3)), 1/10.],  #3
            [ 1/2., 1/(2*sqrt(3)), 3/10.],  #4
            [ 0., 1/sqrt(3), 1/10.],        #5
            [ 3/10., 7/(10*sqrt(3)), 3/10.],#6
            [ 1/2., 1/(10*sqrt(3)), 3/10.], #7
            [ 7/10., 7/(10*sqrt(3)), 3/10.],#8
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 5; 1.0; (0, 0, 0)],
            [2; 5; 1.0; (0, 0, 0)], # zz
            [5; 3; 1.0; (0, 0, 0)],
            [3; 6; 1.0; (0, 0, 0)], # zz
            [6; 4; 1.0; (0, 0, 0)],
            [4; 7; 1.0; (0, 0, 0)], # zz
            [4; 8; 1.0; (0, 0, 0)],

            [5; 1; 1.0; (0, 0, 0)],
            [5; 2; 1.0; (0, 0, 0)], # zz
            [3; 5; 1.0; (0, 0, 0)],
            [6; 3; 1.0; (0, 0, 0)], # zz
            [4; 6; 1.0; (0, 0, 0)],
            [7; 4; 1.0; (0, 0, 0)], # zz
            [8; 4; 1.0; (0, 0, 0)],

            [8; 1; 1.0; (1, 0, 0)],
            [1; 8; 1.0; (-1, 0, 0)],

            [8; 1; 1.0; (1, 0, 1)],   # zz
            [1; 8; 1.0; (-1, 0, -1)], # zz

            [2; 7; 1.0; (0, 1, 0)],
            [7; 2; 1.0; (0, -1, 0)],

            [2; 7; 1.0; (0, 1, -1)],
            [7; 2; 1.0; (0, -1, 1)],

            [3; 6; 1.0; (0, 0, -1)],
            [6; 3; 1.0; (0, 0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_8_3_c_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [1., 0., 0.]
        a2 = [-1/2., sqrt(3)/2., 0.]
        a3 = [0., 0., 2/5.]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [-1/5., 4/(5*sqrt(3)), 1/10.], #1
            [ 0., 7/(5*sqrt(3)), 1/10.], #2
            [ 1/5., 4/(5*sqrt(3)), 1/10.], #3
            [ 1/2., 1/(2*sqrt(3)), 3/10.], #4
            [ 0., 1/sqrt(3), 1/10.], #5
            [ 3/10., 7/(10*sqrt(3)), 3/10.], #6
            [ 1/2., 1/(10*sqrt(3)), 3/10.], #7
            [ 7/10., 7/(10*sqrt(3)), 3/10.], #8
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 5; "tx"; (0, 0, 0)],
            [2; 5; "tz"; (0, 0, 0)], # zz
            [5; 3; "ty"; (0, 0, 0)],
            [3; 6; "tz"; (0, 0, 0)], # zz
            [6; 4; "ty"; (0, 0, 0)],
            [4; 7; "tz"; (0, 0, 0)], # zz
            [4; 8; "tx"; (0, 0, 0)],

            [5; 1; "tx"; (0, 0, 0)],
            [5; 2; "tz"; (0, 0, 0)], # zz
            [3; 5; "ty"; (0, 0, 0)],
            [6; 3; "tz"; (0, 0, 0)], # zz
            [4; 6; "ty"; (0, 0, 0)],
            [7; 4; "tz"; (0, 0, 0)], # zz
            [8; 4; "tx"; (0, 0, 0)],

            [8; 1; "ty"; (1, 0, 0)],
            [1; 8; "ty"; (-1, 0, 0)],

            [8; 1; "tz"; (1, 0, 1)],   # zz
            [1; 8; "tz"; (-1, 0, -1)], # zz

            [2; 7; "tx"; (0, 1, 0)],
            [7; 2; "tx"; (0, -1, 0)],

            [2; 7; "ty"; (0, 1, -1)],
            [7; 2; "ty"; (0, -1, 1)],

            [3; 6; "tx"; (0, 0, -1)],
            [6; 3; "tx"; (0, 0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_8_3_c_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcell_8_3_c



"""
    getUnitcell_8_3_n([version::Int64=1; save::Bool=false])

get the implementation of the *3D (8,3)n lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 4 - simple (DEFAULT) & simple Kitaev

Bravais lattice vectors are

    a1 = [1.0, 0.0, 0.0]
    a2 = [0.0, 1.0, 0.0]
    a3 = [0.5, 0.5, 2/(2*sqrt(3) + sqrt(2))]

16 sites per unitcell, located at

    r01 =         x*a + (0.5 - x)*b +      0.25*c
    r02 =     (1-x)*a + (0.5 - x)*b +      0.25*c
    r03 = (0.5 + x)*a +       0.5*b + (0.5 - z)*c
    r04 =     (1-x)*a + (0.5 + x)*b +      0.25*c
    r05 =         x*a + (0.5 + x)*b +      0.25*c
    r06 = (0.5 - x)*a +       0.5*b + (0.5 - z)*c
    r07 =                   (1-x)*b +         z*c
    r08 =                       x*b +         z*c
    r09 = (0.5 - x)*a +         x*b +      0.25*c
    r10 =       0.5*a + (0.5 - x)*b + (0.5 - z)*c
    r11 = (0.5 + x)*a +         x*b +      0.25*c
    r12 = (0.5 + x)*a +   (1 - x)*b +      0.25*c
    r13 =       0.5*a + (0.5 + x)*b + (0.5 - z)*c
    r14 = (0.5 - x)*a +   (1 - x)*b +      0.25*c
    r15 =         x*a               +         z*c
    r16 =     (1-x)*a               +         z*c

with the additional definition of

    a = [1.0, 0.0, 0.0]
    b = [0.0, 1.0, 0.0]
    c = [0.0, 0.0, 4/(2*sqrt(3) + sqrt(2))]

    x = (sqrt(3) + sqrt(2))/(2*(2*sqrt(3) + sqrt(2)))
    z = 0.125

For version `1`, all couplings have strength 1.0.

For version `4`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`.





# Examples

```julia-repl
julia> unitcell = getUnitcell_8_3_n()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcell_8_3_n(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a = [1.0, 0.0, 0.0]
        b = [0.0, 1.0, 0.0]
        c = [0.0, 0.0, 4/(2*sqrt(3) + sqrt(2))]
        a1 = a
        a2 = b
        a3 = 0.5*(a+b+c)
        x = (sqrt(3) + sqrt(2))/(2*(2*sqrt(3) + sqrt(2)))
        z = 0.125
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
                    x.*a + (0.5 - x).*b +      0.25.*c,
                (1-x).*a + (0.5 - x).*b +      0.25.*c,
            (0.5 + x).*a +       0.5.*b + (0.5 - z).*c,
                (1-x).*a + (0.5 + x).*b +      0.25.*c,
                    x.*a + (0.5 + x).*b +      0.25.*c,
            (0.5 - x).*a +       0.5.*b + (0.5 - z).*c,
                               (1-x).*b +         z.*c,
                                   x.*b +         z.*c,
            (0.5 - x).*a +         x.*b +      0.25.*c,
                  0.5.*a + (0.5 - x).*b + (0.5 - z).*c,
            (0.5 + x).*a +         x.*b +      0.25.*c,
            (0.5 + x).*a +   (1 - x).*b +      0.25.*c,
                  0.5.*a + (0.5 + x).*b + (0.5 - z).*c,
            (0.5 - x).*a +   (1 - x).*b +      0.25.*c,
                    x.*a                +         z.*c,
                (1-x).*a                +         z.*c
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 10; 1.0; (0, 0, 0)],
            [10; 2; 1.0; (0, 0, 0)],
            [2; 11; 1.0; (0, 0, 0)], # zz
            [11; 3; 1.0; (0, 0, 0)],
            [3; 12; 1.0; (0, 0, 0)],
            [12; 4; 1.0; (0, 0, 0)], # zz
            [4; 13; 1.0; (0, 0, 0)],
            [13; 5; 1.0; (0, 0, 0)],
            [5; 14; 1.0; (0, 0, 0)], # zz
            [14; 6; 1.0; (0, 0, 0)],
            [6; 9; 1.0; (0, 0, 0)],
            [9; 8; 1.0; (0, 0, 0)],
            [9; 1; 1.0; (0, 0, 0)], # zz
            [1; 15; 1.0; (0, 0, 0)],
            [2; 16; 1.0; (0, 0, 0)],
            [14; 7; 1.0; (0, 0, 0)],

            [10; 1; 1.0; (0, 0, 0)],
            [2; 10; 1.0; (0, 0, 0)],
            [11; 2; 1.0; (0, 0, 0)], # zz
            [3; 11; 1.0; (0, 0, 0)],
            [12; 3; 1.0; (0, 0, 0)],
            [4; 12; 1.0; (0, 0, 0)], # zz
            [13; 4; 1.0; (0, 0, 0)],
            [5; 13; 1.0; (0, 0, 0)],
            [14; 5; 1.0; (0, 0, 0)], # zz
            [6; 14; 1.0; (0, 0, 0)],
            [9; 6; 1.0; (0, 0, 0)],
            [8; 9; 1.0; (0, 0, 0)],
            [1; 9; 1.0; (0, 0, 0)], # zz
            [15; 1; 1.0; (0, 0, 0)],
            [16; 2; 1.0; (0, 0, 0)],
            [7; 14; 1.0; (0, 0, 0)],

            [11; 8; 1.0; (1, 0, 0)],
            [8; 11; 1.0; (-1, 0, 0)],

            [12; 7; 1.0; (1, 0, 0)],
            [7; 12; 1.0; (-1, 0, 0)],

            [7; 10; 1.0; (0, 1, -1)], # zz
            [10; 7; 1.0; (0, -1, 1)], # zz

            [8; 13; 1.0; (0, 0, -1)], # zz
            [13; 8; 1.0; (0, 0, 1)],  # zz

            [16; 4; 1.0; (0, -1, 0)],
            [4; 16; 1.0; (0, 1, 0)],

            [15; 5; 1.0; (0, -1, 0)],
            [5; 15; 1.0; (0, 1, 0)],

            [15; 3; 1.0; (0, 0, -1)], # zz
            [3; 15; 1.0; (0, 0, 1)],  # zz

            [16; 6; 1.0; (1, 0, -1)], # zz
            [6; 16; 1.0; (-1, 0, 1)]  # zz
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_8_3_n_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a = [1.0, 0.0, 0.0]
        b = [0.0, 1.0, 0.0]
        c = [0.0, 0.0, 4/(2*sqrt(3) + sqrt(2))]
        a1 = a
        a2 = b
        a3 = 0.5*(a+b+c)
        x = (sqrt(3) + sqrt(2))/(2*(2*sqrt(3) + sqrt(2)))
        z = 0.125
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
                    x.*a + (0.5 - x).*b +      0.25.*c,
                (1-x).*a + (0.5 - x).*b +      0.25.*c,
            (0.5 + x).*a +       0.5.*b + (0.5 - z).*c,
                (1-x).*a + (0.5 + x).*b +      0.25.*c,
                    x.*a + (0.5 + x).*b +      0.25.*c,
            (0.5 - x).*a +       0.5.*b + (0.5 - z).*c,
                               (1-x).*b +         z.*c,
                                   x.*b +         z.*c,
            (0.5 - x).*a +         x.*b +      0.25.*c,
                  0.5.*a + (0.5 - x).*b + (0.5 - z).*c,
            (0.5 + x).*a +         x.*b +      0.25.*c,
            (0.5 + x).*a +   (1 - x).*b +      0.25.*c,
                  0.5.*a + (0.5 + x).*b + (0.5 - z).*c,
            (0.5 - x).*a +   (1 - x).*b +      0.25.*c,
                    x.*a                +         z.*c,
                (1-x).*a                +         z.*c
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 10; "tx"; (0, 0, 0)],
            [10; 2; "ty"; (0, 0, 0)],
            [2; 11; "tz"; (0, 0, 0)], # zz
            [11; 3; "tx"; (0, 0, 0)],
            [3; 12; "ty"; (0, 0, 0)],
            [12; 4; "tz"; (0, 0, 0)], # zz
            [4; 13; "tx"; (0, 0, 0)],
            [13; 5; "ty"; (0, 0, 0)],
            [5; 14; "tz"; (0, 0, 0)], # zz
            [14; 6; "tx"; (0, 0, 0)],
            [6; 9; "ty"; (0, 0, 0)],
            [9; 8; "tx"; (0, 0, 0)],
            [9; 1; "tz"; (0, 0, 0)], # zz
            [1; 15; "ty"; (0, 0, 0)],
            [2; 16; "tx"; (0, 0, 0)],
            [14; 7; "ty"; (0, 0, 0)],

            [10; 1; "tx"; (0, 0, 0)],
            [2; 10; "ty"; (0, 0, 0)],
            [11; 2; "tz"; (0, 0, 0)], # zz
            [3; 11; "tx"; (0, 0, 0)],
            [12; 3; "ty"; (0, 0, 0)],
            [4; 12; "tz"; (0, 0, 0)], # zz
            [13; 4; "tx"; (0, 0, 0)],
            [5; 13; "ty"; (0, 0, 0)],
            [14; 5; "tz"; (0, 0, 0)], # zz
            [6; 14; "tx"; (0, 0, 0)],
            [9; 6; "ty"; (0, 0, 0)],
            [8; 9; "tx"; (0, 0, 0)],
            [1; 9; "tz"; (0, 0, 0)], # zz
            [15; 1; "ty"; (0, 0, 0)],
            [16; 2; "tx"; (0, 0, 0)],
            [7; 14; "ty"; (0, 0, 0)],

            [11; 8; "ty"; (1, 0, 0)],
            [8; 11; "ty"; (-1, 0, 0)],

            [12; 7; "tx"; (1, 0, 0)],
            [7; 12; "tx"; (-1, 0, 0)],

            [7; 10; "tz"; (0, 1, -1)], # zz
            [10; 7; "tz"; (0, -1, 1)], # zz

            [8; 13; "tz"; (0, 0, -1)], # zz
            [13; 8; "tz"; (0, 0, 1)],  # zz

            [16; 4; "ty"; (0, -1, 0)],
            [4; 16; "ty"; (0, 1, 0)],

            [15; 5; "tx"; (0, -1, 0)],
            [5; 15; "tx"; (0, 1, 0)],

            [15; 3; "tz"; (0, 0, -1)], # zz
            [3; 15; "tz"; (0, 0, 1)],  # zz

            [16; 6; "tz"; (1, 0, -1)], # zz
            [6; 16; "tz"; (-1, 0, 1)] # zz
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_8_3_n_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcell_8_3_n





"""
    getUnitcell_9_3_a([version::Int64=1; save::Bool=false])

get the implementation of the *3D (9,3)a lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 4 - simple (DEFAULT) & simple Kitaev

Bravais lattice vectors are

    a1 = (-1/3)*a + (1/3)*b + (1/3)*c
    a2 = (-1/3)*a + (2/3)*b + (1/3)*c
    a3 =  (2/3)*a + (1/3)*b + (1/3)*c

12 sites per unitcell, located at

    r01 =    d_f*a
    r02 =  2*d_h*a +   d_h*b + (1/12)*c
    r03 =    d_f*a +   d_f*b
    r04 =    d_h*a + 2*d_h*b - (1/12)*c
    r05 =              d_f*b
    r06 =   -d_h*a +   d_h*b + (1/12)*c
    r07 =   -d_f*a
    r08 = -2*d_h*a -   d_h*b - (1/12)*c
    r09 =   -d_f*a -   d_f*b,
    r10 =   -d_h*a - 2*d_h*b + (1/12)*c
    r11 =          -   d_f*b
    r12 =    d_h*a -   d_h*b - (1/12)*c

with the additional definition of

    a = [ 1.0,        0.0, 0.0]
    b = [-0.5, sqrt(3)/2., 0.0]
    c = [ 0.0,        0.0, sqrt(6*(4 + sqrt(3)))/(1 + 2*sqrt(3))]

    d_f = sqrt(3)/(1+2*sqrt(3))
    d_h = (29 - 3*sqrt(3))/132

For version `1`, all couplings have strength 1.0.

For version `4`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`.




#### 2 & 5 - alternative & alternative Kitaev

Bravais lattice vectors are

    a1 = [-sqrt(3)/2, 1/2, 1/sqrt(3)]
    a2 = [         0,  -1, 1/sqrt(3)]
    a3 = [ sqrt(3)/2, 1/2, 1/sqrt(3)]

12 sites per unitcell, located at

    r01 = [0, 0, 0],
    r02 = (1/6).*a1 - (1/6).*a2,
    r03 = (2/6).*a1 - (2/6).*a2,
    r04 = (3/6).*a1 - (2/6).*a2 - (1/6).*a3,
    r05 = (4/6).*a1 - (2/6).*a2 - (2/6).*a3,
    r06 = (4/6).*a1 - (1/6).*a2 - (3/6).*a3,
    r07 = (4/6).*a1             - (4/6).*a3,
    r08 = (3/6).*a1 + (1/6).*a2 - (4/6).*a3,
    r09 = (2/6).*a1 + (2/6).*a2 - (4/6).*a3,
    r10 = (1/6).*a1 + (2/6).*a2 - (3/6).*a3,
    r11 = (2/6).*a1             - (2/6).*a3,
    r12 = (1/6).*a1             - (1/6).*a3

For version `2`, all couplings have strength 1.0.

For version `5`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`.





# Examples

```julia-repl
julia> unitcell = getUnitcell_9_3_a()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcell_9_3_a(version::Int64=1; save::Bool=false)
    if version==1
        # the lattice vectors
        a = [ 1.0,        0.0, 0.0]
        b = [-0.5, sqrt(3)/2., 0.0]
        c = [ 0.0,        0.0, sqrt(6*(4 + sqrt(3)))/(1 + 2*sqrt(3))]
        a1 = (-1/3.).*a + (1/3.).*b + (1/3.).*c
        a2 = (-1/3.).*a + (2/3.).*b + (1/3.).*c
        a3 =  (2/3.).*a + (1/3.).*b + (1/3.).*c
        d_f = sqrt(3)/(1+2*sqrt(3))
        d_h = (29 - 3*sqrt(3))/132.
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
               d_f.*a,
             2*d_h.*a +   d_h.*b + (1/12).*c,
               d_f.*a +   d_f.*b,
               d_h.*a + 2*d_h.*b - (1/12).*c,
                          d_f.*b,
              -d_h.*a +   d_h.*b + (1/12).*c,
              -d_f.*a,
            -2*d_h.*a -   d_h.*b - (1/12).*c,
              -d_f.*a -   d_f.*b,
              -d_h.*a - 2*d_h.*b + (1/12).*c,
                      -   d_f.*b,
               d_h.*a -   d_h.*b - (1/12).*c
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0, 0, 0)],
            [2; 3; 1.0; (0, 0, 0)],
            [3; 4; 1.0; (0, 0, 0)],
            [4; 5; 1.0; (0, 0, 0)],
            [5; 6; 1.0; (0, 0, 0)],
            [6; 7; 1.0; (0, 0, 0)],
            [7; 8; 1.0; (0, 0, 0)],
            [8; 9; 1.0; (0, 0, 0)],
            [9; 10; 1.0; (0, 0, 0)],
            [10; 11; 1.0; (0, 0, 0)],
            [11; 12; 1.0; (0, 0, 0)],
            [12; 1; 1.0; (0, 0, 0)],

            [2; 1; 1.0; (0, 0, 0)],
            [3; 2; 1.0; (0, 0, 0)],
            [4; 3; 1.0; (0, 0, 0)],
            [5; 4; 1.0; (0, 0, 0)],
            [6; 5; 1.0; (0, 0, 0)],
            [7; 6; 1.0; (0, 0, 0)],
            [8; 7; 1.0; (0, 0, 0)],
            [9; 8; 1.0; (0, 0, 0)],
            [10; 9; 1.0; (0, 0, 0)],
            [11; 10; 1.0; (0, 0, 0)],
            [12; 11; 1.0; (0, 0, 0)],
            [1; 12; 1.0; (0, 0, 0)],


            [3; 9; 1.0; (0, -1, 1)], # zz
            [9; 3; 1.0; (0, 1, -1)], # zz

            [1; 7; 1.0; (-1, 0, 1)], # zz
            [7; 1; 1.0; (1, 0, -1)], # zz

            [5; 11; 1.0; (1, -1, 0)], # zz
            [11; 5; 1.0; (-1, 1, 0)], # zz

            [12; 6; 1.0; (-1, 0, 0)], # zz
            [6; 12; 1.0; (1, 0, 0)], # zz

            [8; 2; 1.0; (0, 0, -1)], # zz
            [2; 8; 1.0; (0, 0, 1)], # zz

            [4; 10; 1.0; (0, -1, 0)], # zz
            [10; 4; 1.0; (0, 1, 0)]   # zz

        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_9_3_a_unitcell.jld"
    elseif version == 2
        # the lattice vectors
        a1 = [-sqrt(3)/2., 1/2., 1/sqrt(3)]
        a2 = [         0.,  -1., 1/sqrt(3)]
        a3 = [ sqrt(3)/2., 1/2., 1/sqrt(3)]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0, 0, 0],
            (1/6).*a1 - (1/6).*a2,
            (2/6).*a1 - (2/6).*a2,
            (3/6).*a1 - (2/6).*a2 - (1/6).*a3,
            (4/6).*a1 - (2/6).*a2 - (2/6).*a3,
            (4/6).*a1 - (1/6).*a2 - (3/6).*a3,
            (4/6).*a1             - (4/6).*a3,
            (3/6).*a1 + (1/6).*a2 - (4/6).*a3,
            (2/6).*a1 + (2/6).*a2 - (4/6).*a3,
            (1/6).*a1 + (2/6).*a2 - (3/6).*a3,
            (2/6).*a1             - (2/6).*a3,
            (1/6).*a1             - (1/6).*a3
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0, 0, 0)],
            [2; 3; 1.0; (0, 0, 0)],
            [3; 4; 1.0; (0, 0, 0)],
            [4; 5; 1.0; (0, 0, 0)],
            [5; 6; 1.0; (0, 0, 0)],
            [6; 7; 1.0; (0, 0, 0)],
            [7; 8; 1.0; (0, 0, 0)],
            [8; 9; 1.0; (0, 0, 0)],
            [9; 10; 1.0; (0, 0, 0)],
            [10; 11; 1.0; (0, 0, 0)],
            [11; 12; 1.0; (0, 0, 0)],
            [12; 1; 1.0; (0, 0, 0)],

            [2; 1; 1.0; (0, 0, 0)],
            [3; 2; 1.0; (0, 0, 0)],
            [4; 3; 1.0; (0, 0, 0)],
            [5; 4; 1.0; (0, 0, 0)],
            [6; 5; 1.0; (0, 0, 0)],
            [7; 6; 1.0; (0, 0, 0)],
            [8; 7; 1.0; (0, 0, 0)],
            [9; 8; 1.0; (0, 0, 0)],
            [10; 9; 1.0; (0, 0, 0)],
            [11; 10; 1.0; (0, 0, 0)],
            [12; 11; 1.0; (0, 0, 0)],
            [1; 12; 1.0; (0, 0, 0)],

            [1; 7; 1.0; (-1, 0, 1)],
            [7; 1; 1.0; (1, 0, -1)],
            [2; 8; 1.0; (0, 0, 1)],
            [8; 2; 1.0; (0, 0, -1)],
            [3; 9; 1.0; (0, -1, 1)],
            [9; 3; 1.0; (0, 1, -1)],
            [4; 10; 1.0; (0, -1, 0)],
            [10; 4; 1.0; (0, 1, 0)],
            [5; 11; 1.0; (1, -1, 0)],
            [11; 5; 1.0; (-1, 1, 0)],
            [6; 12; 1.0; (1, 0, 0)],
            [12; 6; 1.0; (-1, 0, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_9_3_a_V2_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a = [ 1.0,        0.0, 0.0]
        b = [-0.5, sqrt(3)/2., 0.0]
        c = [ 0.0,        0.0, sqrt(6*(4 + sqrt(3)))/(1 + 2*sqrt(3))]
        a1 = (-1/3.).*a + (1/3.).*b + (1/3.).*c
        a2 = (-1/3.).*a + (2/3.).*b + (1/3.).*c
        a3 =  (2/3.).*a + (1/3.).*b + (1/3.).*c
        d_f = sqrt(3)/(1+2*sqrt(3))
        d_h = (29 - 3*sqrt(3))/132.
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
               d_f.*a,
             2*d_h.*a +   d_h.*b + (1/12).*c,
               d_f.*a +   d_f.*b,
               d_h.*a + 2*d_h.*b - (1/12).*c,
                          d_f.*b,
              -d_h.*a +   d_h.*b + (1/12).*c,
              -d_f.*a,
            -2*d_h.*a -   d_h.*b - (1/12).*c,
              -d_f.*a -   d_f.*b,
              -d_h.*a - 2*d_h.*b + (1/12).*c,
                      -   d_f.*b,
               d_h.*a -   d_h.*b - (1/12).*c
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; "tx"; (0, 0, 0)],
            [2; 3; "ty"; (0, 0, 0)],
            [3; 4; "tx"; (0, 0, 0)],
            [4; 5; "ty"; (0, 0, 0)],
            [5; 6; "tx"; (0, 0, 0)],
            [6; 7; "ty"; (0, 0, 0)],
            [7; 8; "tx"; (0, 0, 0)],
            [8; 9; "ty"; (0, 0, 0)],
            [9; 10; "tx"; (0, 0, 0)],
            [10; 11; "ty"; (0, 0, 0)],
            [11; 12; "tx"; (0, 0, 0)],
            [12; 1; "ty"; (0, 0, 0)],

            [2; 1; "tx"; (0, 0, 0)],
            [3; 2; "ty"; (0, 0, 0)],
            [4; 3; "tx"; (0, 0, 0)],
            [5; 4; "ty"; (0, 0, 0)],
            [6; 5; "tx"; (0, 0, 0)],
            [7; 6; "ty"; (0, 0, 0)],
            [8; 7; "tx"; (0, 0, 0)],
            [9; 8; "ty"; (0, 0, 0)],
            [10; 9; "tx"; (0, 0, 0)],
            [11; 10; "ty"; (0, 0, 0)],
            [12; 11; "tx"; (0, 0, 0)],
            [1; 12; "ty"; (0, 0, 0)],


            [3; 9; "tz"; (0, -1, 1)], # zz
            [9; 3; "tz"; (0, 1, -1)], # zz

            [1; 7; "tz"; (-1, 0, 1)], # zz
            [7; 1; "tz"; (1, 0, -1)], # zz

            [5; 11; "tz"; (1, -1, 0)], # zz
            [11; 5; "tz"; (-1, 1, 0)], # zz

            [12; 6; "tz"; (-1, 0, 0)], # zz
            [6; 12; "tz"; (1, 0, 0)], # zz

            [8; 2; "tz"; (0, 0, -1)], # zz
            [2; 8; "tz"; (0, 0, 1)], # zz

            [4; 10; "tz"; (0, -1, 0)], # zz
            [10; 4; "tz"; (0, 1, 0)]   # zz

        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_9_3_a_kitaev_unitcell.jld"
    elseif version == 5
        # the lattice vectors
        a1 = [-sqrt(3)/2., 1/2., 1/sqrt(3)]
        a2 = [         0.,  -1., 1/sqrt(3)]
        a3 = [ sqrt(3)/2., 1/2., 1/sqrt(3)]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0, 0, 0],
            (1/6).*a1 - (1/6).*a2,
            (2/6).*a1 - (2/6).*a2,
            (3/6).*a1 - (2/6).*a2 - (1/6).*a3,
            (4/6).*a1 - (2/6).*a2 - (2/6).*a3,
            (4/6).*a1 - (1/6).*a2 - (3/6).*a3,
            (4/6).*a1             - (4/6).*a3,
            (3/6).*a1 + (1/6).*a2 - (4/6).*a3,
            (2/6).*a1 + (2/6).*a2 - (4/6).*a3,
            (1/6).*a1 + (2/6).*a2 - (3/6).*a3,
            (2/6).*a1             - (2/6).*a3,
            (1/6).*a1             - (1/6).*a3
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; "tx"; (0, 0, 0)],
            [2; 3; "ty"; (0, 0, 0)],
            [3; 4; "tx"; (0, 0, 0)],
            [4; 5; "ty"; (0, 0, 0)],
            [5; 6; "tx"; (0, 0, 0)],
            [6; 7; "ty"; (0, 0, 0)],
            [7; 8; "tx"; (0, 0, 0)],
            [8; 9; "ty"; (0, 0, 0)],
            [9; 10; "tx"; (0, 0, 0)],
            [10; 11; "ty"; (0, 0, 0)],
            [11; 12; "tx"; (0, 0, 0)],
            [12; 1; "ty"; (0, 0, 0)],

            [2; 1; "tx"; (0, 0, 0)],
            [3; 2; "ty"; (0, 0, 0)],
            [4; 3; "tx"; (0, 0, 0)],
            [5; 4; "ty"; (0, 0, 0)],
            [6; 5; "tx"; (0, 0, 0)],
            [7; 6; "ty"; (0, 0, 0)],
            [8; 7; "tx"; (0, 0, 0)],
            [9; 8; "ty"; (0, 0, 0)],
            [10; 9; "tx"; (0, 0, 0)],
            [11; 10; "ty"; (0, 0, 0)],
            [12; 11; "tx"; (0, 0, 0)],
            [1; 12; "ty"; (0, 0, 0)],

            [1; 7; "tz"; (-1, 0, 1)],
            [7; 1; "tz"; (1, 0, -1)],
            [2; 8; "tz"; (0, 0, 1)],
            [8; 2; "tz"; (0, 0, -1)],
            [3; 9; "tz"; (0, -1, 1)],
            [9; 3; "tz"; (0, 1, -1)],
            [4; 10; "tz"; (0, -1, 0)],
            [10; 4; "tz"; (0, 1, 0)],
            [5; 11; "tz"; (1, -1, 0)],
            [11; 5; "tz"; (-1, 1, 0)],
            [6; 12; "tz"; (1, 0, 0)],
            [12; 6; "tz"; (-1, 0, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_9_3_a_V2_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcell_9_3_a



"""
    getUnitcell_10_3_a([version::Int64=1; save::Bool=false])
    getUnitcellHyperoctagon([version::Int64=1; save::Bool=false])

get the implementation of the *3D (10,3)a / hyperoctagon lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 4 - simple (DEFAULT) & simple Kitaev

Bravais lattice vectors are

    a1 = [1.0,  0.0,  0.0]
    a2 = [0.5,  0.5, -0.5]
    a3 = [0.5,  0.5,  0.5]

4 sites per unitcell, located at

    r1 = [1/8,  1/8,  1/8]
    r2 = [5/8,  3/8, -1/8]
    r3 = [3/8,  1/8, -1/8]
    r4 = [7/8,  3/8,  1/8]

Note that the unitcell in this implementation is not bipartite (although the lattice itself is bipartite)!

For version `1`, all couplings have strength 1.0.

For version `4`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`.







# Examples

```julia-repl
julia> unitcell = getUnitcell_10_3_a()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcell_10_3_a(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = [1, 0, 0]
        a2 = [0.5, 0.5, -0.5]
        a3 = [0.5, 0.5, 0.5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [1/8., 1/8.,  1/8.],
            [5/8., 3/8., -1/8.],
            [3/8., 1/8., -1/8.],
            [7/8., 3/8.,  1/8.]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 3; 1.0; (0, 0, 0)],
            [3; 2; 1.0; (0, 0, 0)],  # zz
            [2; 4; 1.0; (0, 0, 0)],

            [3; 1; 1.0; (0, 0, 0)],
            [2; 3; 1.0; (0, 0, 0)],  # zz
            [4; 2; 1.0; (0, 0, 0)],

            [4; 1; 1.0; (1, 0, 0)],  # zz
            [1; 4; 1.0; (-1, 0, 0)], # zz

            [2; 1; 1.0; (0, 1, 0)],
            [1; 2; 1.0; (0, -1, 0)],

            [3; 4; 1.0; (0, 0, -1)],
            [4; 3; 1.0; (0, 0, 1)],
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_hyperoctagon_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [1, 0, 0]
        a2 = [0.5, 0.5, -0.5]
        a3 = [0.5, 0.5, 0.5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [1/8., 1/8.,  1/8.],
            [5/8., 3/8., -1/8.],
            [3/8., 1/8., -1/8.],
            [7/8., 3/8.,  1/8.]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 3; "tx"; (0, 0, 0)],
            [3; 2; "tz"; (0, 0, 0)],  # zz
            [2; 4; "tx"; (0, 0, 0)],

            [3; 1; "tx"; (0, 0, 0)],
            [2; 3; "tz"; (0, 0, 0)],  # zz
            [4; 2; "tx"; (0, 0, 0)],

            [4; 1; "tz"; (1, 0, 0)],  # zz
            [1; 4; "tz"; (-1, 0, 0)], # zz

            [2; 1; "ty"; (0, 1, 0)],
            [1; 2; "ty"; (0, -1, 0)],

            [3; 4; "ty"; (0, 0, -1)],
            [4; 3; "ty"; (0, 0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_hyperoctagon_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
getUnitcellHyperoctagon = getUnitcell_10_3_a

export getUnitcell_10_3_a
export getUnitcellHyperoctagon






"""
    getUnitcell_10_3_b([version::Int64=1; save::Bool=false])
    getUnitcellHyperhoneycomb([version::Int64=1; save::Bool=false])

get the implementation of the *3D (10,3)b / hyperhoneycomb lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 4 - simple (DEFAULT) & simple Kitaev

Bravais lattice vectors are

    a1 = [-1.0,  1.0, -2.0]
    a2 = [-1.0,  1.0,  2.0]
    a3 = [ 2.0,  4.0,  0.0]

4 sites per unitcell, located at

    r1 = [0.0,  0.0,  0.0]
    r2 = [1.0,  1.0,  0.0]
    r3 = [1.0,  2.0,  1.0]
    r4 = [0.0, -1.0,  1.0]

For version `1`, all couplings have strength 1.0.

For version `4`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`.




#### 2 & 5 - shifted & shifted Kitaev

Bravais lattice vectors are

    a1 = [-1.0,  1.0, -2.0]
    a2 = [-1.0,  1.0,  2.0]
    a3 = [ 2.0,  4.0,  0.0]

4 sites per unitcell, located at

    r1 = [0.0,  0.0,  0.0]
    r2 = [1.0,  1.0,  0.0]
    r3 = [1.0,  2.0,  1.0]
    r4 = [2.0,  3.0,  1.0]

which is the same definition as 1 & 4 but with site 4 shifted through the unitcell along a3 direction.

For version `2`, all couplings have strength 1.0.

For version `5`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`.





# Examples

```julia-repl
julia> unitcell = getUnitcell_10_3_b()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcell_10_3_b(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = [-1.0,  1.0, -2.0]
        a2 = [-1.0,  1.0,  2.0]
        a3 = [ 2.0,  4.0,  0.0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 2.0, 1.0],
            [0.0, -1.0, 1.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0, 0, 0)],
            [1; 4; 1.0; (0, 0, 0)],
            [1; 4; 1.0; (1, 0, 0)],

            [2; 1; 1.0; (0, 0, 0)],
            [2; 3; 1.0; (0, 0, 0)],
            [2; 3; 1.0; (0, -1, 0)],

            [3; 2; 1.0; (0, 0, 0)],
            [3; 2; 1.0; (0, 1, 0)],
            [3; 4; 1.0; (0, 0, 1)],

            [4; 1; 1.0; (0, 0, 0)],
            [4; 1; 1.0; (-1, 0, 0)],
            [4; 3; 1.0; (0, 0, -1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_hyperhoneycomb_v1_unitcell.jld"
    elseif version == 2
        # the lattice vectors
        a1 = [-1.0,  1.0, -2.0]
        a2 = [-1.0,  1.0,  2.0]
        a3 = [ 2.0,  4.0,  0.0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 2.0, 1.0],
            [2.0, 3.0, 1.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0, 0, 0)],
            [1; 4; 1.0; (0, 0, -1)],
            [1; 4; 1.0; (1, 0, -1)],

            [2; 1; 1.0; (0, 0, 0)],
            [2; 3; 1.0; (0, 0, 0)],
            [2; 3; 1.0; (0, -1, 0)],

            [3; 2; 1.0; (0, 0, 0)],
            [3; 2; 1.0; (0, 1, 0)],
            [3; 4; 1.0; (0, 0, 0)],

            [4; 1; 1.0; (0, 0, 1)],
            [4; 1; 1.0; (-1, 0, 1)],
            [4; 3; 1.0; (0, 0, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_hyperhoneycomb_v2_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [-1.0,  1.0, -2.0]
        a2 = [-1.0,  1.0,  2.0]
        a3 = [ 2.0,  4.0,  0.0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 2.0, 1.0],
            [0.0, -1.0, 1.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; "tz"; (0, 0, 0)],
            [1; 4; "tx"; (0, 0, 0)],
            [1; 4; "ty"; (1, 0, 0)],

            [2; 1; "tz"; (0, 0, 0)],
            [2; 3; "tx"; (0, 0, 0)],
            [2; 3; "ty"; (0, -1, 0)],

            [3; 2; "tx"; (0, 0, 0)],
            [3; 2; "ty"; (0, 1, 0)],
            [3; 4; "tz"; (0, 0, 1)],

            [4; 1; "tx"; (0, 0, 0)],
            [4; 1; "ty"; (-1, 0, 0)],
            [4; 3; "tz"; (0, 0, -1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_hyperhoneycomb_kitaev_unitcell.jld"
    elseif version == 5
        # the lattice vectors
        a1 = [-1.0,  1.0, -2.0]
        a2 = [-1.0,  1.0,  2.0]
        a3 = [ 2.0,  4.0,  0.0]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 2.0, 1.0],
            [2.0, 3.0, 1.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; "tz"; (0, 0, 0)],
            [1; 4; "tx"; (0, 0, -1)],
            [1; 4; "ty"; (1, 0, -1)],

            [2; 1; "tz"; (0, 0, 0)],
            [2; 3; "tx"; (0, 0, 0)],
            [2; 3; "ty"; (0, -1, 0)],

            [3; 2; "tx"; (0, 0, 0)],
            [3; 2; "ty"; (0, 1, 0)],
            [3; 4; "tz"; (0, 0, 0)],

            [4; 1; "tx"; (0, 0, 1)],
            [4; 1; "ty"; (-1, 0, 1)],
            [4; 3; "tz"; (0, 0, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_hyperhoneycomb_v2_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
getUnitcellHyperhoneycomb = getUnitcell_10_3_b

export getUnitcellHyperhoneycomb
export getUnitcell_10_3_b



"""
    getUnitcell_10_3_c([version::Int64=1; save::Bool=false])

get the implementation of the *3D (10,3)c lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 4 - simple (DEFAULT) & simple Kitaev

Bravais lattice vectors are

    a1 = [ 1.0,       0.0,         0.0]
    a2 = [-0.5, sqrt(3)/2,         0.0]
    a3 = [ 0.0,       0.0, 3*sqrt(3)/2]

6 sites per unitcell, located at

    r1 = [0.25, 1/(4*sqrt(3)), 1/(2*sqrt(3))],
    r2 = [0.75, 1/(4*sqrt(3)), 2/sqrt(3)],
    r3 = [0.5,  1/sqrt(3),     7/(2*sqrt(3))],
    r4 = [0.75, 1/(4*sqrt(3)), 1/sqrt(3)],
    r5 = [0.5,  1/sqrt(3),     5/(2*sqrt(3))],
    r6 = [0.25, 1/(4*sqrt(3)), 4/sqrt(3)]

For version `1`, all couplings have strength 1.0.

For version `4`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`.







# Examples

```julia-repl
julia> unitcell = getUnitcell_10_3_c()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcell_10_3_c(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = [1, 0, 0]
        a2 = [-0.5, sqrt(3)/2., 0.0]
        a3 = [0.0, 0.0, (3*sqrt(3))/2.]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.25, 1/(4*sqrt(3)), 1/(2*sqrt(3))],
            [0.75, 1/(4*sqrt(3)), 2/sqrt(3)],
            [0.5, 1/sqrt(3), 7/(2*sqrt(3))],
            [0.75, 1/(4*sqrt(3)), 1/sqrt(3)],
            [0.5, 1/sqrt(3), 5/(2*sqrt(3))],
            [0.25, 1/(4*sqrt(3)), 4/sqrt(3)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 4; 1.0; (0, 0, 0)],
            [4; 2; 1.0; (0, 0, 0)],  # zz
            [2; 5; 1.0; (0, 0, 0)],
            [5; 3; 1.0; (0, 0, 0)],  # zz
            [3; 6; 1.0; (0, 0, 0)],

            [4; 1; 1.0; (0, 0, 0)],
            [2; 4; 1.0; (0, 0, 0)],  # zz
            [5; 2; 1.0; (0, 0, 0)],
            [3; 5; 1.0; (0, 0, 0)],  # zz
            [6; 3; 1.0; (0, 0, 0)],

            [4; 1; 1.0; (1, 0, 0)],
            [1; 4; 1.0; (-1, 0, 0)],

            [5; 2; 1.0; (0, 1, 0)],
            [2; 5; 1.0; (0, -1, 0)],

            [3; 6; 1.0; (1, 1, 0)],
            [6; 3; 1.0; (-1, -1, 0)],

            [1; 6; 1.0; (0, 0, -1)],  # zz
            [6; 1; 1.0; (0, 0, 1)], # zz
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_10_3_c_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [1, 0, 0]
        a2 = [-0.5, sqrt(3)/2., 0.0]
        a3 = [0.0, 0.0, (3*sqrt(3))/2.]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.25, 1/(4*sqrt(3)), 1/(2*sqrt(3))],
            [0.75, 1/(4*sqrt(3)), 2/sqrt(3)],
            [0.5, 1/sqrt(3), 7/(2*sqrt(3))],
            [0.75, 1/(4*sqrt(3)), 1/sqrt(3)],
            [0.5, 1/sqrt(3), 5/(2*sqrt(3))],
            [0.25, 1/(4*sqrt(3)), 4/sqrt(3)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 4; "tx"; (0, 0, 0)],
            [4; 2; "tz"; (0, 0, 0)],  # zz
            [2; 5; "tx"; (0, 0, 0)],
            [5; 3; "tz"; (0, 0, 0)],  # zz
            [3; 6; "tx"; (0, 0, 0)],

            [4; 1; "tx"; (0, 0, 0)],
            [2; 4; "tz"; (0, 0, 0)],  # zz
            [5; 2; "tx"; (0, 0, 0)],
            [3; 5; "tz"; (0, 0, 0)],  # zz
            [6; 3; "tx"; (0, 0, 0)],

            [4; 1; "ty"; (1, 0, 0)],
            [1; 4; "ty"; (-1, 0, 0)],

            [5; 2; "ty"; (0, 1, 0)],
            [2; 5; "ty"; (0, -1, 0)],

            [3; 6; "ty"; (1, 1, 0)],
            [6; 3; "ty"; (-1, -1, 0)],

            [1; 6; "tz"; (0, 0, -1)],  # zz
            [6; 1; "tz"; (0, 0, 1)] # zz
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_10_3_c_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcell_10_3_c



"""
    getUnitcell_10_3_d([version::Int64=1; save::Bool=false])

get the implementation of the *3D (10,3)d lattice* unitcell. The `version` integer
specifies the exact implementation convention that is used and the boolean `save`
can be passed if one wants to save the unitcell after creation.



# Versions

#### 1 & 4 - simple (DEFAULT) & simple Kitaev

Bravais lattice vectors are

    a1 = [ 1.0,       0.0,         0.0]
    a2 = [-0.5, sqrt(3)/2,         0.0]
    a3 = [ 0.0,       0.0, 3*sqrt(3)/2]

8 sites per unitcell, located at

    r1 = [0.0,     -a, 0.75*c]
    r2 = [ -a,    0.0,  0.5*c]
    r3 = [0.0,      a, 0.25*c]
    r4 = [  a,    0.0,    0.0]
    r5 = [ -a,   -0.5, 0.25*c]
    r6 = [0.0,  a-0.5,  0.5*c]
    r7 = [  a,   -0.5, 0.75*c]
    r8 = [0.0, -a-0.5,    0.0]

with the additional definition of

    a = 0.25*(2 - sqrt(2))
    c = 0.5

For version `1`, all couplings have strength 1.0.

For version `4`, the couplings follow the Kitaev scheme and are labeled `"tx"`, `"ty"` and `"tz"`.







# Examples

```julia-repl
julia> unitcell = getUnitcell_10_3_d()
LatticePhysics.Unitcell(...)
```
"""
function getUnitcell_10_3_d(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a = 0.25*(2 - sqrt(2))
        c = 0.5
        a1 = [0.5, -0.5,  0.0]
        a2 = [0.5,  0.5,  0.0]
        a3 = [0.0,  0.0,  0.5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0,     -a, 0.75*c],
            [ -a,    0.0,  0.5*c],
            [0.0,      a, 0.25*c],
            [  a,    0.0,    0.0],
            [ -a,   -0.5, 0.25*c],
            [0.0,  a-0.5,  0.5*c],
            [  a,   -0.5, 0.75*c],
            [0.0, -a-0.5,    0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; 1.0; (0, 0, 0)],
            [3; 4; 1.0; (0, 0, 0)],
            [5; 6; 1.0; (0, 0, 0)],

            [2; 1; 1.0; (0, 0, 0)],
            [4; 3; 1.0; (0, 0, 0)],
            [6; 5; 1.0; (0, 0, 0)],

            [2; 3; 1.0; (0, 0, 0)],
            [3; 2; 1.0; (0, 0, 0)],

            [6; 7; 1.0; (0, 0, 0)],
            [7; 6; 1.0; (0, 0, 0)],

            [1; 6; 1.0; (0, 0, 0)], #zz
            [6; 1; 1.0; (0, 0, 0)], #zz

            [4; 5; 1.0; (0, 1, 0)], #zz
            [5; 4; 1.0; (0, -1, 0)], #zz

            [8; 5; 1.0; (0, 0, 0)],
            [5; 8; 1.0; (0, 0, 0)],

            [8; 7; 1.0; (0, 0, -1)],
            [7; 8; 1.0; (0, 0, 1)],

            [2; 7; 1.0; (-1, 0, 0)], # zz
            [7; 2; 1.0; (1, 0, 0)],  # zz

            [3; 8; 1.0; (-1, 1, 0)], # zz
            [8; 3; 1.0; (1, -1, 0)],  # zz

            [1; 4; 1.0; (0, 0, 1)],
            [4; 1; 1.0; (0, 0, -1)],

        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_10_3_d_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a = 0.25*(2 - sqrt(2))
        c = 0.5
        a1 = [0.5, -0.5,  0.0]
        a2 = [0.5,  0.5,  0.0]
        a3 = [0.0,  0.0,  0.5]
        lattice_vectors = Array{Float64, 1}[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array{Float64, 1}[
            [0.0,     -a, 0.75*c],
            [ -a,    0.0,  0.5*c],
            [0.0,      a, 0.25*c],
            [  a,    0.0,    0.0],
            [ -a,   -0.5, 0.25*c],
            [0.0,  a-0.5,  0.5*c],
            [  a,   -0.5, 0.75*c],
            [0.0, -a-0.5,    0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array{Any, 1}[
            [1; 2; "tx"; (0, 0, 0)],
            [3; 4; "tx"; (0, 0, 0)],
            [5; 6; "ty"; (0, 0, 0)],

            [2; 1; "tx"; (0, 0, 0)],
            [4; 3; "tx"; (0, 0, 0)],
            [6; 5; "ty"; (0, 0, 0)],

            [2; 3; "ty"; (0, 0, 0)],
            [3; 2; "ty"; (0, 0, 0)],

            [6; 7; "tx"; (0, 0, 0)],
            [7; 6; "tx"; (0, 0, 0)],

            [1; 6; "tz"; (0, 0, 0)], #zz
            [6; 1; "tz"; (0, 0, 0)], #zz

            [4; 5; "tz"; (0, 1, 0)], #zz
            [5; 4; "tz"; (0, -1, 0)], #zz

            [8; 5; "tx"; (0, 0, 0)],
            [5; 8; "tx"; (0, 0, 0)],

            [8; 7; "ty"; (0, 0, -1)],
            [7; 8; "ty"; (0, 0, 1)],

            [2; 7; "tz"; (-1, 0, 0)], # zz
            [7; 2; "tz"; (1, 0, 0)],  # zz

            [3; 8; "tz"; (-1, 1, 0)], # zz
            [8; 3; "tz"; (1, -1, 0)],  # zz

            [1; 4; "ty"; (0, 0, 1)],
            [4; 1; "ty"; (0, 0, -1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_10_3_d_kitaev_unitcell.jld"
    end
    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
end
export getUnitcell_10_3_d
