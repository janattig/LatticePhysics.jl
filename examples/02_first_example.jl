# include the package
using LatticePhysics

# obtain the unitcell of the honeycomb lattice
unitcell = getUnitcellHoneycomb()

# build a lattice for visual confirmation that this is the correct unitcell
lattice = getLatticeInBox(unitcell, [10.0, 8.0])
# plot the lattice for visual confirmation
plotLattice(lattice, openfile=true)

# get a path for band structure calculations, use the default path
path = DEFAULT_PATH_TRIANGULAR
# calculate the band structure (using the unitcell)
calculateBandStructureAlongPath(unitcell, path, showPlot=true)