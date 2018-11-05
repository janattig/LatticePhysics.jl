# LatticePhysics

Lattice based calculations and plotting for Julia v.1.0.*


## Usage

See the [wiki pages](http://gitsrv.thp.uni-koeln.de/attig/LatticePhysics.jl/wikis/home).


## Installation:

You can use the package via the package mode in Julia (Pkg). However, since the package
is not listed in the Julia package repositories, you have to use
```julia
pkg(1.0)> add "git@gitsrv.thp.uni-koeln.de:attig/LatticePhysics.jl.git"
```
Note: this can lead to Errors under Windows 10 due to incorrect SSH access. Use the following command instead:
```julia
pkg(1.0)> add "http://gitsrv.thp.uni-koeln.de/attig/LatticePhysics.jl.git"
```
You will be prompted a username and password validation but it should work the same way.
