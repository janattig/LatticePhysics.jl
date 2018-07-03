# include the package
using LatticePhysics

# obtain the unitcell of the square octagon lattice (kitaev bonds)
unitcell = getUnitcellSquareOctagon(5)

# build a lattice for visual confirmation that this is the correct unitcell
lattice = getLatticeInBox(unitcell, [12.0, 6.0])
# plot the lattice for visual confirmation (including bond coloring)
plotLattice2D(
    lattice,
    openfile=true,
    colorcode_bonds=Dict("tx"=>[255,0,0], "ty"=>[175, 215, 175], "tz"=>[175, 185, 215]),
    site_radius=-0.25,
    bond_thickness=-0.25
)
