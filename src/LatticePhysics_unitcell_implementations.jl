################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT UNITCELLS in 2D and 3D
#
#   STRUCTURE OF THE FILE
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [0.5, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [0.5, 0.0],
            [0.5, 0.5],
            [0.0, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [0.5, 0.0],
            [0.0, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [0.5, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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

#### 1 - simple (DEFAULT)

Bravais lattice vectors are

    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]

4 sites per unitcell, located at

    r1 = [0.0, 0.0]
    r2 = [0.5, 0.0]
    r3 = [0.5, 0.5]
    r4 = [0.0, 0.5]

The coupling values are `["t1", ..., "t5"]` for the distinct bonds.



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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [0.5, 0.0],
            [0.5, 0.5],
            [0.0, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        filename = "$(FOLDER_UNITCELLS)2d_shastry_sutherland_unitcell.jld"
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [0.0, 1/sqrt(3.0)],
			[-1/sqrt(3.0), 0.0],
			[-1/sqrt(3.0), 1/sqrt(3.0)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [ 0.5, -0.5],
            [ 0.5,  0.5],
            [-0.5, -0.5],
            [-0.5,  0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [0.0, 1/sqrt(3.0)],
			[-1/sqrt(3.0), 0.0],
			[-1/sqrt(3.0), 1/sqrt(3.0)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [ 0.5, -0.5],
            [ 0.5,  0.5],
            [-0.5, -0.5],
            [-0.5,  0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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

    # generate unitcell
    uc = Unitcell(lattice_vectors, basis, connections, filename)
    if save
        saveUnitcell(uc)
    end
    # return the unitcell
    return uc
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
function getUnitcellFullyConnectedSquare(version::Int64=1; save::Bool=true, J1=1.0, J1X=0.5)
    if version == 1
        if J1==1.0 && J1X==0.5
            # the lattice vectors
            a1 = [1.0, 0.0]
            a2 = [0.0, 1.0]
            lattice_vectors = Array[]
            push!(lattice_vectors, a1)
            push!(lattice_vectors, a2)
            # Basis Definition
            basis = Array[
                [0.0, 0.0]
            ]
            # Connection Definition
            # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            connections = Array[
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
            lattice_vectors = Array[]
            push!(lattice_vectors, a1)
            push!(lattice_vectors, a2)
            # Basis Definition
            basis = Array[
                [0.0, 0.0]
            ]
            # Connection Definition
            # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
            connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [1/sqrt(3.0), 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
            [0.0, 0.0],
            [sqrt(3.0)/4, -0.25],
            [sqrt(3.0)/4, +0.25]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        b1 = [0.0, 0.0]
        b2 = [1/sqrt(3.0), 0.0]
        basis = Array[
            b1 .+ 1/4 * (b2 .- b1),
            b1 .+ 1/4 * (b2 .- a2),
            b1 .+ 1/4 * (b2 .- a1),
            b2 .+ 1/4 * (b1 .- b2),
            b2 .+ 1/4 * (a2 .- b2),
            b2 .+ 1/4 * (a1 .- b2)
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        b1 = [0.0, 0.0]
        b2 = [1/sqrt(3.0), 0.0]
        basis = Array[
            b1 .+ 1/4 * (b2 .- b1),
            b1 .+ 1/4 * (b2 .- a2),
            b1 .+ 1/4 * (b2 .- a1),
            b2 .+ 1/4 * (b1 .- b2),
            b2 .+ 1/4 * (a2 .- b2),
            b2 .+ 1/4 * (a1 .- b2)
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
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
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        # Basis Definition
        basis = Array[
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
        connections = Array[
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [0.25, 0.25, 0.25]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_diamond_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_diamond_$(J1)_unitcell.jld"
        end
    elseif version == 2
        # the lattice vectors
        a1 = 0.5 .* [0, 1, 1]
        a2 = 0.5 .* [1, 0, 1]
        a3 = 0.5 .* [1, 1, 0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [0.25, 0.25, 0.25]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        if J1==1.0 && J2==0.5
            filename = "$(FOLDER_UNITCELLS)3d_diamond_NN_NNN_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_diamond_NN_NNN_$(J1)_$(J2)_unitcell.jld"
        end
    elseif version == 3
        # the lattice vectors
        a1 = 0.5 .* [0, 1, 1]
        a2 = 0.5 .* [1, 0, 1]
        a3 = 0.5 .* [1, 1, 0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [0.25, 0.25, 0.25]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        if J1==1.0 && J21=="J21" && J22=="J22"
            filename = "$(FOLDER_UNITCELLS)3d_diamond_aniso_NNN_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_diamond_aniso_NNN_$(J1)_$(J21)_$(J22)_unitcell.jld"
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
export getUnitcellDiamond


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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_bcc_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_bcc_$(J1)_unitcell.jld"
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0., 0., 0.],
            [0., 0.25, 0.25],
            [0.25, 0., 0.25],
            [0.25, 0.25, 0.]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_pyrochlore_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_pyrochlore_$(J1)_unitcell.jld"
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.5, sqrt(3)/10., 0.0],
            [3/5., sqrt(3)/5., (2*sqrt(2))/5.],
            [0.1, (3*sqrt(3))/10., sqrt(2)/5.],
            [0.4, sqrt(3)/5., sqrt(2)/5.],
            [0.0, (2*sqrt(3))/5., 0.0],
            [-0.1, (3*sqrt(3))/10., (2*sqrt(2))/5.]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_8_3_a_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_8_3_a_$(J1)_unitcell.jld"
        end
    elseif version == 4
        # the lattice vectors
        a1 = [1.0, 0.0, 0.0]
        a2 = [-0.5, sqrt(3)/2., 0.0]
        a3 = [0.0, 0.0, (3*sqrt(2))/5.]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.5, sqrt(3)/10., 0.0],
            [3/5., sqrt(3)/5., (2*sqrt(2))/5.],
            [0.1, (3*sqrt(3))/10., sqrt(2)/5.],
            [0.4, sqrt(3)/5., sqrt(2)/5.],
            [0.0, (2*sqrt(3))/5., 0.0],
            [-0.1, (3*sqrt(3))/10., (2*sqrt(2))/5.]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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

#-----------------------------------------------------------------------------------------------------------------------------
# LATTICE (8,3)b
# 1 - simple, 6 sites per UC, all connections have interaction strength J1
# 4 - Kitaev
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcell_8_3_b(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = [1/2., 1/(2*sqrt(3)), sqrt(2)/(5*sqrt(3))]
        a2 = [0, 1/sqrt(3), (2*sqrt(2))/(5*sqrt(3))]
        a3 = [0, 0, sqrt(6)/5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [1/10., 1/(2*sqrt(3)), sqrt(2)/(5*sqrt(3))],
            [1/5., sqrt(3)/5, sqrt(6)/5],
            [3/10., 11/(10*sqrt(3)), (4*sqrt(2))/(5*sqrt(3))],
            [1/5., 2/(5*sqrt(3)), (2*sqrt(2))/(5*sqrt(3))],
            [3/10., (3*sqrt(3))/10., sqrt(6)/5],
            [2/5., 1/sqrt(3), sqrt(2)/sqrt(3)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 4; J1; (0, 0, 0)], # zz
            [4; 2; J1; (0, 0, 0)],
            [2; 5; J1; (0, 0, 0)], # zz
            [5; 3; J1; (0, 0, 0)],
            [3; 6; J1; (0, 0, 0)], # zz

            [4; 1; J1; (0, 0, 0)], # zz
            [2; 4; J1; (0, 0, 0)],
            [5; 2; J1; (0, 0, 0)], # zz
            [3; 5; J1; (0, 0, 0)],
            [6; 3; J1; (0, 0, 0)], # zz

            [6; 1; J1; (1, 0, 1)],
            [1; 6; J1; (-1, 0, -1)],

            [4; 3; J1; (0, -1, 0)],
            [3; 4; J1; (0, 1, 0)],

            [1; 2; J1; (0, 0, -1)],
            [2; 1; J1; (0, 0, 1)],

            [5; 6; J1; (0, 0, -1)],
            [6; 5; J1; (0, 0, 1)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_8_3_b_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a1 = [1/2., 1/(2*sqrt(3)), sqrt(2)/(5*sqrt(3))]
        a2 = [0, 1/sqrt(3), (2*sqrt(2))/(5*sqrt(3))]
        a3 = [0, 0, sqrt(6)/5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [1/10., 1/(2*sqrt(3)), sqrt(2)/(5*sqrt(3))],
            [1/5., sqrt(3)/5, sqrt(6)/5],
            [3/10., 11/(10*sqrt(3)), (4*sqrt(2))/(5*sqrt(3))],
            [1/5., 2/(5*sqrt(3)), (2*sqrt(2))/(5*sqrt(3))],
            [3/10., (3*sqrt(3))/10., sqrt(6)/5],
            [2/5., 1/sqrt(3), sqrt(2)/sqrt(3)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
    saveUnitcell(uc)
    # return the unitcell
    return uc
end
export getUnitcell_8_3_b

#-----------------------------------------------------------------------------------------------------------------------------
# LATTICE (8,3)c
# 1 - simple, 8 sites per UC, all connections have interaction strength J1
# 4 - Kitaev
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcell_8_3_c(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = [1., 0., 0.]
        a2 = [-1/2., sqrt(3)/2., 0.]
        a3 = [0., 0., 2/5.]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
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
        connections = Array[
            [1; 5; J1; (0, 0, 0)],
            [2; 5; J1; (0, 0, 0)], # zz
            [5; 3; J1; (0, 0, 0)],
            [3; 6; J1; (0, 0, 0)], # zz
            [6; 4; J1; (0, 0, 0)],
            [4; 7; J1; (0, 0, 0)], # zz
            [4; 8; J1; (0, 0, 0)],

            [5; 1; J1; (0, 0, 0)],
            [5; 2; J1; (0, 0, 0)], # zz
            [3; 5; J1; (0, 0, 0)],
            [6; 3; J1; (0, 0, 0)], # zz
            [4; 6; J1; (0, 0, 0)],
            [7; 4; J1; (0, 0, 0)], # zz
            [8; 4; J1; (0, 0, 0)],

            [8; 1; J1; (1, 0, 0)],
            [1; 8; J1; (-1, 0, 0)],

            [8; 1; J1; (1, 0, 1)],   # zz
            [1; 8; J1; (-1, 0, -1)], # zz

            [2; 7; J1; (0, 1, 0)],
            [7; 2; J1; (0, -1, 0)],

            [2; 7; J1; (0, 1, -1)],
            [7; 2; J1; (0, -1, 1)],

            [3; 6; J1; (0, 0, -1)],
            [6; 3; J1; (0, 0, 1)]
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_8_3_c_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_8_3_c_$(J1)_unitcell.jld"
        end
    elseif version == 4
        # the lattice vectors
        a1 = [1., 0., 0.]
        a2 = [-1/2., sqrt(3)/2., 0.]
        a3 = [0., 0., 2/5.]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
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
        connections = Array[
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

#-----------------------------------------------------------------------------------------------------------------------------
# LATTICE (8,3)n
# 1 - simple, 16 sites per UC, all connections have interaction strength J1
# 4 - Kitaev
#-----------------------------------------------------------------------------------------------------------------------------
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            x*a + (0.5 - x)*b + c/4.,
            (1-x)*a + (0.5 - x)*b + c/4.,
            (0.5 + x)*a + b/2. + (0.5 - z)*c,
            (1-x)*a + (0.5 + x)*b + c/4.,
            x*a + (0.5 + x)*b + c/4.,
            (0.5 - x)*a + b/2. + (0.5 - z)*c,
            (1-x)*b + z*c,
            x*b + z*c,
            (0.5 - x)*a + x*b + c/4.,
            a/2. + (0.5 - x)*b + (0.5 - z)*c,
            (0.5 + x)*a + x*b + c/4.,
            (0.5 + x)*a + (1 - x)*b + c/4.,
            a/2. + (0.5 + x)*b + (0.5 - z)*c,
            (0.5 - x)*a + (1 - x)*b + c/4.,
            x*a + z*c,
            (1-x)*a + z*c
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 10; J1; (0, 0, 0)],
            [10; 2; J1; (0, 0, 0)],
            [2; 11; J1; (0, 0, 0)], # zz
            [11; 3; J1; (0, 0, 0)],
            [3; 12; J1; (0, 0, 0)],
            [12; 4; J1; (0, 0, 0)], # zz
            [4; 13; J1; (0, 0, 0)],
            [13; 5; J1; (0, 0, 0)],
            [5; 14; J1; (0, 0, 0)], # zz
            [14; 6; J1; (0, 0, 0)],
            [6; 9; J1; (0, 0, 0)],
            [9; 8; J1; (0, 0, 0)],
            [9; 1; J1; (0, 0, 0)], # zz
            [1; 15; J1; (0, 0, 0)],
            [2; 16; J1; (0, 0, 0)],
            [14; 7; J1; (0, 0, 0)],

            [10; 1; J1; (0, 0, 0)],
            [2; 10; J1; (0, 0, 0)],
            [11; 2; J1; (0, 0, 0)], # zz
            [3; 11; J1; (0, 0, 0)],
            [12; 3; J1; (0, 0, 0)],
            [4; 12; J1; (0, 0, 0)], # zz
            [13; 4; J1; (0, 0, 0)],
            [5; 13; J1; (0, 0, 0)],
            [14; 5; J1; (0, 0, 0)], # zz
            [6; 14; J1; (0, 0, 0)],
            [9; 6; J1; (0, 0, 0)],
            [8; 9; J1; (0, 0, 0)],
            [1; 9; J1; (0, 0, 0)], # zz
            [15; 1; J1; (0, 0, 0)],
            [16; 2; J1; (0, 0, 0)],
            [7; 14; J1; (0, 0, 0)],

            [11; 8; J1; (1, 0, 0)],
            [8; 11; J1; (-1, 0, 0)],

            [12; 7; J1; (1, 0, 0)],
            [7; 12; J1; (-1, 0, 0)],

            [7; 10; J1; (0, 1, -1)], # zz
            [10; 7; J1; (0, -1, 1)], # zz

            [8; 13; J1; (0, 0, -1)], # zz
            [13; 8; J1; (0, 0, 1)],  # zz

            [16; 4; J1; (0, -1, 0)],
            [4; 16; J1; (0, 1, 0)],

            [15; 5; J1; (0, -1, 0)],
            [5; 15; J1; (0, 1, 0)],

            [15; 3; J1; (0, 0, -1)], # zz
            [3; 15; J1; (0, 0, 1)],  # zz

            [16; 6; J1; (1, 0, -1)], # zz
            [6; 16; J1; (-1, 0, 1)]  # zz
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_8_3_n_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_8_3_n_$(J1)_unitcell.jld"
        end
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
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            x*a + (0.5 - x)*b + c/4.,
            (1-x)*a + (0.5 - x)*b + c/4.,
            (0.5 + x)*a + b/2. + (0.5 - z)*c,
            (1-x)*a + (0.5 + x)*b + c/4.,
            x*a + (0.5 + x)*b + c/4.,
            (0.5 - x)*a + b/2. + (0.5 - z)*c,
            (1-x)*b + z*c,
            x*b + z*c,
            (0.5 - x)*a + x*b + c/4.,
            a/2. + (0.5 - x)*b + (0.5 - z)*c,
            (0.5 + x)*a + x*b + c/4.,
            (0.5 + x)*a + (1 - x)*b + c/4.,
            a/2. + (0.5 + x)*b + (0.5 - z)*c,
            (0.5 - x)*a + (1 - x)*b + c/4.,
            x*a + z*c,
            (1-x)*a + z*c

        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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


#-----------------------------------------------------------------------------------------------------------------------------
# LATTICE (9,3)a
# 1 - simple, 12 sites per UC, all connections have interaction strength J1
# 2 - simple shifted
# 4 - Kitaev
# 5 - Kitaev shifted
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcell_9_3_a(version::Int64=1; save::Bool=false)
    if version==1
        # the lattice vectors
        a = [1.0, 0.0, 0.0]
        b = [-0.5, sqrt(3)/2., 0.0]
        c = [0.0, 0.0, sqrt(6*(4 + sqrt(3)))/(1 + 2*sqrt(3))]
        a1 = -a/3. + b/3. + c/3.
        a2 = -a/3. -2*b/3. + c/3.
        a3 = 2*a/3. + b/3. + c/3.
        d_f = sqrt(3)/(1+2*sqrt(3))
        d_h = (29 - 3*sqrt(3))/132.
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            d_f * a,
            2*d_h * a + d_h * b + c/12.,
            d_f * (a + b),
            d_h * a + 2*d_h * b - c/12.,
            d_f * b,
            -d_h *a + d_h * b + c/12.,
            -d_f * a,
            -2*d_h * a - d_h * b - c/12.,
            -d_f * (a + b),
            -d_h * a - 2*d_h*b + c/12.,
            -d_f * b,
            d_h * a - d_h * b - c/12.

        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)],
            [2; 3; J1; (0, 0, 0)],
            [3; 4; J1; (0, 0, 0)],
            [4; 5; J1; (0, 0, 0)],
            [5; 6; J1; (0, 0, 0)],
            [6; 7; J1; (0, 0, 0)],
            [7; 8; J1; (0, 0, 0)],
            [8; 9; J1; (0, 0, 0)],
            [9; 10; J1; (0, 0, 0)],
            [10; 11; J1; (0, 0, 0)],
            [11; 12; J1; (0, 0, 0)],
            [12; 1; J1; (0, 0, 0)],

            [2; 1; J1; (0, 0, 0)],
            [3; 2; J1; (0, 0, 0)],
            [4; 3; J1; (0, 0, 0)],
            [5; 4; J1; (0, 0, 0)],
            [6; 5; J1; (0, 0, 0)],
            [7; 6; J1; (0, 0, 0)],
            [8; 7; J1; (0, 0, 0)],
            [9; 8; J1; (0, 0, 0)],
            [10; 9; J1; (0, 0, 0)],
            [11; 10; J1; (0, 0, 0)],
            [12; 11; J1; (0, 0, 0)],
            [1; 12; J1; (0, 0, 0)],


            [3; 9; J1; (0, -1, 1)], # zz
            [9; 3; J1; (0, 1, -1)], # zz

            [1; 7; J1; (-1, 0, 1)], # zz
            [7; 1; J1; (1, 0, -1)], # zz

            [5; 11; J1; (1, -1, 0)], # zz
            [11; 5; J1; (-1, 1, 0)], # zz

            [12; 6; J1; (-1, 0, 0)], # zz
            [6; 12; J1; (1, 0, 0)], # zz

            [8; 2; J1; (0, 0, -1)], # zz
            [2; 8; J1; (0, 0, 1)], # zz

            [4; 10; J1; (0, -1, 0)], # zz
            [10; 4; J1; (0, 1, 0)]   # zz

        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_9_3_a_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_9_3_a_$(J1)_unitcell.jld"
        end
    elseif version == 2
        # the lattice vectors
        a1 = [-sqrt(3)/2., 1/2., 1/sqrt(3)]
        a2 = [0, -1, 1/sqrt(3)]
        a3 = [sqrt(3)/2, 1/2, 1/sqrt(3)]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0, 0, 0],
            a1/6. - a2/6.,
            a1/3. - a2/3.,
            a1/2. - a2/3. - a3/6.,
            2/3.*a1 - a2/3. - a3/3.,
            2/3.*a1 - a2/6. - a3/2.,
            2/3.*a1 - 2/3.*a3,
            a1/2. + a2/6. - 2/3.*a3,
            a1/3. + a2/3. - 2/3.*a3,
            a1/6. + a2/3. - a3/2.,
            a2/3. - a3/3.,
            a2/6. - a3/6.
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)],
            [2; 3; J1; (0, 0, 0)],
            [3; 4; J1; (0, 0, 0)],
            [4; 5; J1; (0, 0, 0)],
            [5; 6; J1; (0, 0, 0)],
            [6; 7; J1; (0, 0, 0)],
            [7; 8; J1; (0, 0, 0)],
            [8; 9; J1; (0, 0, 0)],
            [9; 10; J1; (0, 0, 0)],
            [10; 11; J1; (0, 0, 0)],
            [11; 12; J1; (0, 0, 0)],
            [12; 1; J1; (0, 0, 0)],

            [2; 1; J1; (0, 0, 0)],
            [3; 2; J1; (0, 0, 0)],
            [4; 3; J1; (0, 0, 0)],
            [5; 4; J1; (0, 0, 0)],
            [6; 5; J1; (0, 0, 0)],
            [7; 6; J1; (0, 0, 0)],
            [8; 7; J1; (0, 0, 0)],
            [9; 8; J1; (0, 0, 0)],
            [10; 9; J1; (0, 0, 0)],
            [11; 10; J1; (0, 0, 0)],
            [12; 11; J1; (0, 0, 0)],
            [1; 12; J1; (0, 0, 0)],

            [1; 7; J1; (-1, 0, 1)],
            [7; 1; J1; (1, 0, -1)],
            [2; 8; J1; (0, 0, 1)],
            [8; 2; J1; (0, 0, -1)],
            [3; 9; J1; (0, -1, 1)],
            [9; 3; J1; (0, 1, -1)],
            [4; 10; J1; (0, -1, 0)],
            [10; 4; J1; (0, 1, 0)],
            [5; 11; J1; (1, -1, 0)],
            [11; 5; J1; (-1, 1, 0)],
            [6; 12; J1; (1, 0, 0)],
            [12; 6; J1; (-1, 0, 0)]
        ]
        # filename
        filename = "$(FOLDER_UNITCELLS)3d_9_3_a_V2_unitcell.jld"
    elseif version == 4
        # the lattice vectors
        a = [1.0, 0.0, 0.0]
        b = [-0.5, sqrt(3)/2., 0.0]
        c = [0.0, 0.0, sqrt(6*(4 + sqrt(3)))/(1 + 2*sqrt(3))]
        a1 = -a/3. + b/3. + c/3.
        a2 = -a/3. -2*b/3. + c/3.
        a3 = 2*a/3. + b/3. + c/3.
        d_f = sqrt(3)/(1+2*sqrt(3))
        d_h = (29 - 3*sqrt(3))/132.
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            d_f * a,
            2*d_h * a + d_h * b + c/12.,
            d_f * (a + b),
            d_h * a + 2*d_h * b - c/12.,
            d_f * b,
            -d_h *a + d_h * b + c/12.,
            -d_f * a,
            -2*d_h * a - d_h * b - c/12.,
            -d_f * (a + b),
            -d_h * a - 2*d_h*b + c/12.,
            -d_f * b,
            d_h * a - d_h * b - c/12.

        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        a2 = [0, -1, 1/sqrt(3)]
        a3 = [sqrt(3)/2, 1/2, 1/sqrt(3)]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0, 0, 0],
            a1/6. - a2/6.,
            a1/3. - a2/3.,
            a1/2. - a2/3. - a3/6.,
            2/3.*a1 - a2/3. - a3/3.,
            2/3.*a1 - a2/6. - a3/2.,
            2/3.*a1 - 2/3.*a3,
            a1/2. + a2/6. - 2/3.*a3,
            a1/3. + a2/3. - 2/3.*a3,
            a1/6. + a2/3. - a3/2.,
            a2/3. - a3/3.,
            a2/6. - a3/6.
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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


#-----------------------------------------------------------------------------------------------------------------------------
# HYPEROCTAGON LATTICE (10,3)a
# 1 - simple, 4 sites per UC, all connections have interaction strength J1
# 4 - Kitaev
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellHyperoctagon(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = [1, 0, 0]
        a2 = [0.5, 0.5, -0.5]
        a3 = [0.5, 0.5, 0.5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.125, 0.125, 0.125],
            [5/8., 3/8., -1/8.],
            [3/8., 1/8., -1/8.],
            [7/8., 3/8., 1/8.]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 3; J1; (0, 0, 0)],
            [3; 2; J1; (0, 0, 0)],  # zz
            [2; 4; J1; (0, 0, 0)],

            [3; 1; J1; (0, 0, 0)],
            [2; 3; J1; (0, 0, 0)],  # zz
            [4; 2; J1; (0, 0, 0)],

            [4; 1; J1; (1, 0, 0)],  # zz
            [1; 4; J1; (-1, 0, 0)], # zz

            [2; 1; J1; (0, 1, 0)],
            [1; 2; J1; (0, -1, 0)],

            [3; 4; J1; (0, 0, -1)],
            [4; 3; J1; (0, 0, 1)],
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_hyperoctagon_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_hyperoctagon_$(J1)_unitcell.jld"
        end
    elseif version == 4
        # the lattice vectors
        a1 = [1, 0, 0]
        a2 = [0.5, 0.5, -0.5]
        a3 = [0.5, 0.5, 0.5]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.125, 0.125, 0.125],
            [5/8., 3/8., -1/8.],
            [3/8., 1/8., -1/8.],
            [7/8., 3/8., 1/8.]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
function getUnitcell_10_3_a(version::Int64=1; save::Bool=false)
    return getUnitcellHyperoctagon(version, save=save, J1=J1)
end
export getUnitcellHyperoctagon
export getUnitcell_10_3_a

#-----------------------------------------------------------------------------------------------------------------------------
# HYPERHONEYCOMB LATTICE (10,3)b
# 1 - simple, 4 sites per UC, all connections have interaction strength J1
# 2 - simple, 4 sites per UC, all connections have interaction strength J1, site 4 shifted through the unitcell
# 4 - Kitaev
# 5 - Kitaev shifted unitcell
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcellHyperhoneycomb(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = [-1, 1, -2]
        a2 = [-1, 1, 2]
        a3 = [2, 4, 0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 2.0, 1.0],
            [0.0, -1.0, 1.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)],
            [1; 4; J1; (0, 0, 0)],
            [1; 4; J1; (1, 0, 0)],

            [2; 1; J1; (0, 0, 0)],
            [2; 3; J1; (0, 0, 0)],
            [2; 3; J1; (0, -1, 0)],

            [3; 2; J1; (0, 0, 0)],
            [3; 2; J1; (0, 1, 0)],
            [3; 4; J1; (0, 0, 1)],

            [4; 1; J1; (0, 0, 0)],
            [4; 1; J1; (-1, 0, 0)],
            [4; 3; J1; (0, 0, -1)]
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_hyperhoneycomb_v1_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_hyperhoneycomb_v1_$(J1)_unitcell.jld"
        end
    elseif version == 2
        # the lattice vectors
        a1 = [-1, 1, -2]
        a2 = [-1, 1, 2]
        a3 = [2, 4, 0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 2.0, 1.0],
            [2.0, 3.0, 1.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)],
            [1; 4; J1; (0, 0, -1)],
            [1; 4; J1; (1, 0, -1)],

            [2; 1; J1; (0, 0, 0)],
            [2; 3; J1; (0, 0, 0)],
            [2; 3; J1; (0, -1, 0)],

            [3; 2; J1; (0, 0, 0)],
            [3; 2; J1; (0, 1, 0)],
            [3; 4; J1; (0, 0, 0)],

            [4; 1; J1; (0, 0, 1)],
            [4; 1; J1; (-1, 0, 1)],
            [4; 3; J1; (0, 0, 0)]
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_hyperhoneycomb_v2_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_hyperhoneycomb_v2_$(J1)_unitcell.jld"
        end
    elseif version == 4
        # the lattice vectors
        a1 = [-1, 1, -2]
        a2 = [-1, 1, 2]
        a3 = [2, 4, 0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 2.0, 1.0],
            [0.0, -1.0, 1.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
        a1 = [-1, 1, -2]
        a2 = [-1, 1, 2]
        a3 = [2, 4, 0]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 2.0, 1.0],
            [2.0, 3.0, 1.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
function getUnitcell_10_3_b(version::Int64=1; save::Bool=false)
    return getUnitcellHyperhoneycomb(version, save=save, J1=J1)
end
export getUnitcellHyperhoneycomb
export getUnitcell_10_3_b

#-----------------------------------------------------------------------------------------------------------------------------
# LATTICE (10,3)c
# 1 - simple, 6 sites per UC, all connections have interaction strength J1
# 4 - Kitaev
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcell_10_3_c(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a1 = [1, 0, 0]
        a2 = [-0.5, sqrt(3)/2., 0.0]
        a3 = [0.0, 0.0, (3*sqrt(3))/2.]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.25, 1/(4*sqrt(3)), 1/(2*sqrt(3))],
            [0.75, 1/(4*sqrt(3)), 2/sqrt(3)],
            [0.5, 1/sqrt(3), 7/(2*sqrt(3))],
            [0.75, 1/(4*sqrt(3)), 1/sqrt(3)],
            [0.5, 1/sqrt(3), 5/(2*sqrt(3))],
            [0.25, 1/(4*sqrt(3)), 4/sqrt(3)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 4; J1; (0, 0, 0)],
            [4; 2; J1; (0, 0, 0)],  # zz
            [2; 5; J1; (0, 0, 0)],
            [5; 3; J1; (0, 0, 0)],  # zz
            [3; 6; J1; (0, 0, 0)],

            [4; 1; J1; (0, 0, 0)],
            [2; 4; J1; (0, 0, 0)],  # zz
            [5; 2; J1; (0, 0, 0)],
            [3; 5; J1; (0, 0, 0)],  # zz
            [6; 3; J1; (0, 0, 0)],

            [4; 1; J1; (1, 0, 0)],
            [1; 4; J1; (-1, 0, 0)],

            [5; 2; J1; (0, 1, 0)],
            [2; 5; J1; (0, -1, 0)],

            [3; 6; J1; (1, 1, 0)],
            [6; 3; J1; (-1, -1, 0)],

            [1; 6; J1; (0, 0, -1)],  # zz
            [6; 1; J1; (0, 0, 1)], # zz
        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_10_3_c_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_10_3_c_$(J1)_unitcell.jld"
        end
    elseif version == 4
        # the lattice vectors
        a1 = [1, 0, 0]
        a2 = [-0.5, sqrt(3)/2., 0.0]
        a3 = [0.0, 0.0, (3*sqrt(3))/2.]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.25, 1/(4*sqrt(3)), 1/(2*sqrt(3))],
            [0.75, 1/(4*sqrt(3)), 2/sqrt(3)],
            [0.5, 1/sqrt(3), 7/(2*sqrt(3))],
            [0.75, 1/(4*sqrt(3)), 1/sqrt(3)],
            [0.5, 1/sqrt(3), 5/(2*sqrt(3))],
            [0.25, 1/(4*sqrt(3)), 4/sqrt(3)]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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

#-----------------------------------------------------------------------------------------------------------------------------
# LATTICE (10,3)d
# 1 - simple, 8 sites per UC, all connections have interaction strength J1
# 4 - Kitaev
#-----------------------------------------------------------------------------------------------------------------------------
function getUnitcell_10_3_d(version::Int64=1; save::Bool=false)
    if version == 1
        # the lattice vectors
        a = 0.25*(2 - sqrt(2))
        c = 0.5
        a1 = [1, 0, 0]
        a2 = [0.5, 0.5, 0.0]
        a3 = [0.0, 0.0, c]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1 - a2)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, -a, 0.75*c],
            [-a, 0.0, 0.5*c],
            [0, a, 0.25*c],
            [a, 0.0, 0.0],
            [-a, -0.5, 0.25*c],
            [0, a - 0.5, 0.5*c],
            [a, -0.5, 0.75*c],
            [0, -a - 0.5, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
            [1; 2; J1; (0, 0, 0)],
            [3; 4; J1; (0, 0, 0)],
            [5; 6; J1; (0, 0, 0)],

            [2; 1; J1; (0, 0, 0)],
            [4; 3; J1; (0, 0, 0)],
            [6; 5; J1; (0, 0, 0)],

            [2; 3; J1; (0, 0, 0)],
            [3; 2; J1; (0, 0, 0)],

            [6; 7; J1; (0, 0, 0)],
            [7; 6; J1; (0, 0, 0)],

            [1; 6; J1; (0, 0, 0)], #zz
            [6; 1; J1; (0, 0, 0)], #zz

            [4; 5; J1; (0, 1, 0)], #zz
            [5; 4; J1; (0, -1, 0)], #zz

            [8; 5; J1; (0, 0, 0)],
            [5; 8; J1; (0, 0, 0)],

            [8; 7; J1; (0, 0, -1)],
            [7; 8; J1; (0, 0, 1)],

            [2; 7; J1; (-1, 0, 0)], # zz
            [7; 2; J1; (1, 0, 0)],  # zz

            [3; 8; J1; (-1, 1, 0)], # zz
            [8; 3; J1; (1, -1, 0)],  # zz

            [1; 4; J1; (0, 0, 1)],
            [4; 1; J1; (0, 0, -1)],

        ]
        # filename
        if J1==1.0
            filename = "$(FOLDER_UNITCELLS)3d_10_3_d_unitcell.jld"
        else
            filename = "$(FOLDER_UNITCELLS)3d_10_3_d_$(J1)_unitcell.jld"
        end
    elseif version == 4
        # the lattice vectors
        a = 0.25*(2 - sqrt(2))
        c = 0.5
        a1 = [1, 0, 0]
        a2 = [0.5, 0.5, 0.0]
        a3 = [0.0, 0.0, c]
        lattice_vectors = Array[]
        push!(lattice_vectors, a1 - a2)
        push!(lattice_vectors, a2)
        push!(lattice_vectors, a3)
        # Basis Definition
        basis = Array[
            [0.0, -a, 0.75*c],
            [-a, 0.0, 0.5*c],
            [0, a, 0.25*c],
            [a, 0.0, 0.0],
            [-a, -0.5, 0.25*c],
            [0, a - 0.5, 0.5*c],
            [a, -0.5, 0.75*c],
            [0, -a - 0.5, 0.0]
        ]
        # Connection Definition
        # [<from index>; <to index>; <strength>; (<lattice displaced by lattice vector j>)]
        connections = Array[
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
