# LatticePhysics [![pipeline status](http://gitsrv.thp.uni-koeln.de/attig/LatticePhysics.jl/badges/dev-julia-1.0/pipeline.svg)](http://gitsrv.thp.uni-koeln.de/attig/LatticePhysics.jl/commits/dev-julia-1.0) [![coverage report](http://gitsrv.thp.uni-koeln.de/attig/LatticePhysics.jl/badges/dev-julia-1.0/coverage.svg)](http://gitsrv.thp.uni-koeln.de/attig/LatticePhysics.jl/commits/dev-julia-1.0)


Lattice based calculations and plotting for Julia v.1.0.*


## Usage

See the [wiki pages](http://gitsrv.thp.uni-koeln.de/attig/LatticePhysics.jl/wikis/home).


## Installation:

You can use the package via the package mode in Julia (Pkg). However, since the package
is not listed in the Julia package repositories, you have to use
```julia
(v1.0) pkg> add "git@gitsrv.thp.uni-koeln.de:attig/LatticePhysics.jl.git"
```
Note: this can lead to Errors under Windows 10 due to incorrect SSH access. Use the following command instead:
```julia
(v1.0) pkg> add "http://gitsrv.thp.uni-koeln.de/attig/LatticePhysics.jl.git"
```
You will be prompted a username and password validation but it should work the same way.


## Sub-modules

The main code is divided into various sub-modules, which can be found here:
1.  [LatPhysBase](http://gitsrv.thp.uni-koeln.de/attig/LatPhysBase.jl.git)
2.  [LatPhysUnitcellLibrary](http://gitsrv.thp.uni-koeln.de/attig/LatPhysUnitcellLibrary.jl.git)
3.  [LatPhysLatticeConstruction](http://gitsrv.thp.uni-koeln.de/attig/LatPhysLatticeConstruction.jl.git)
4.  [LatPhysReciprocal](http://gitsrv.thp.uni-koeln.de/attig/LatPhysReciprocal.jl.git)
