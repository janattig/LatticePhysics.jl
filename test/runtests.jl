using LatticePhysics
using Test
using Pkg

@testset "LatticePhysics.jl" begin

end


# test all submodules
Pkg.test("LatPhysBase")
Pkg.test("LatPhysUnitcellLibrary")
Pkg.test("LatPhysLatticeConstruction")
Pkg.test("LatPhysLatticeModification")
Pkg.test("LatPhysReciprocal")
Pkg.test("LatPhysBandstructures")
Pkg.test("LatPhysLuttingerTisza")
