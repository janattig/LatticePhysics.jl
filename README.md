![LatticePhysics.jl logo](https://github.com/janattig/LatticePhysics.jl/blob/master/logos/LatticePhysics_banner.png)


Lattice based calculations and plotting for Julia v.1.0.*


## Usage

To use the lattice library in your code, install it as specified below.
To use it in your code, start by importing the library
```julia-REPL
julia> using LatticePhysics
```




## Installation

You can install the package via the package mode in Julia (Pkg). However, since the package
is not listed in the Julia package repositories, you have to first install the unregistered
dependencies manually with
```julia
(v1.0) pkg> add "https://github.com/janattig/LatPhysBase.jl"
(v1.0) pkg> add "https://github.com/janattig/LatPhysUnitcellLibrary.jl"
(v1.0) pkg> add "https://github.com/janattig/LatPhysLatticeConstruction.jl"
(v1.0) pkg> add "https://github.com/janattig/LatPhysLatticeModification.jl"
(v1.0) pkg> add "https://github.com/janattig/LatPhysReciprocal.jl"
(v1.0) pkg> add "https://github.com/janattig/LatPhysLuttingerTisza.jl"
```
to finally install the main package with
```julia
(v1.0) pkg> add "https://github.com/janattig/LatticePhysics.jl"
```



## Sub-modules

The main code is divided into various sub-modules, which can be found here:
1.  [LatPhysBase](https://github.com/janattig/LatPhysBase.jl.git)
2.  [LatPhysUnitcellLibrary](https://github.com/janattig/LatPhysUnitcellLibrary.jl.git)
3.  [LatPhysLatticeConstruction](https://github.com/janattig/LatPhysLatticeConstruction.jl.git)
4.  [LatPhysLatticeModification](https://github.com/janattig/LatPhysLatticeModification.jl.git)
5.  [LatPhysReciprocal](https://github.com/janattig/LatPhysReciprocal.jl.git)
6.  [LatPhysLuttingerTisza](https://github.com/janattig/LatPhysLuttingerTisza.jl.git)

Additionally, there are modules to provide plotting functions, which can be found here:
1.  [LatPhysPlottingPyPlot](https://github.com/janattig/LatPhysPlottingPyPlot.jl.git)
2.  [LatPhysPlottingSVG](https://github.com/janattig/LatPhysPlottingSVG.jl.git)
3.  [LatPhysReciprocalPlottingPyPlot](https://github.com/janattig/LatPhysReciprocalPlottingPyPlot.jl.git)
4.  [LatPhysLuttingerTiszaPlottingPyPlot](https://github.com/janattig/LatPhysLuttingerTiszaPlottingPyPlot.jl.git)

The plotting modules have to be used explicitly inside user code as they are not a direct
part of `LatticePhysics.jl`.
