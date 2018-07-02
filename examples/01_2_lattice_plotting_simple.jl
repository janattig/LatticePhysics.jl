# include the package
using LatticePhysics

# obtain the unitcell of the square octagon lattice (kitaev bonds)
unitcell = getUnitcellSquareOctagon(5)

# build a lattice for visual confirmation that this is the correct unitcell
lattice = getLatticeByBondDistance(unitcell, 7)
# plot the lattice for visual confirmation (including bond coloring)
plotLattice2D(lattice, openfile=true, colorcode_bonds_automation="COLOR")
