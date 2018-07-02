# include the package
using LatticePhysics

# obtain the unitcell of the honeycomb lattice
unitcell = getUnitcellHoneycomb()

# build a lattice for visual confirmation that this is the correct unitcell
lattice = getLatticeInBox(unitcell, [10.0, 8.0])
# plot the lattice for visual confirmation
plotLattice2D(lattice, openfile=true)
