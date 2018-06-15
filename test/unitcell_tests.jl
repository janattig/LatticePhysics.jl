# ALL TESTS THAT CONCERN UNITCELLS

# include
using LatticePhysics
using Base.Test







################################################################################
#
#   IMPLEMENTATIONS OF DIFFERENT UNITCELLS in 2D and 3D
#
################################################################################



#   2D UNITCELLS
#   - SQUARE / RECTANGLE
#       - getUnitcellSquare
@test testUnitcell(getUnitcellSquare(), 2, 2)
@test testUnitcell(getUnitcellSquare(1), 2, 2)
@test testUnitcell(getUnitcellSquare(2), 2, 2)
@test testUnitcell(getUnitcellSquare(3), 2, 2)
#       - getUnitcellExtendedSquare
@test testUnitcell(getUnitcellExtendedSquare(), 2, 2)
@test testUnitcell(getUnitcellExtendedSquare(1), 2, 2)
#       - getUnitcellCheckerboard
@test testUnitcell(getUnitcellCheckerboard(), 2, 2)
@test testUnitcell(getUnitcellCheckerboard(1), 2, 2)
#       - getUnitcellShastrySutherland
@test testUnitcell(getUnitcellShastrySutherland(), 2, 2)
@test testUnitcell(getUnitcellShastrySutherland(1), 2, 2)
@test testUnitcell(getUnitcellShastrySutherland(4), 2, 2)
#       - getUnitcellAdvancedSquare
@test testUnitcell(getUnitcellAdvancedSquare(), 2, 2)
@test testUnitcell(getUnitcellAdvancedSquare(1), 2, 2)
#       - getUnitcellSquareOctagon
@test testUnitcell(getUnitcellSquareOctagon(), 2, 2)
@test testUnitcell(getUnitcellSquareOctagon(1), 2, 2)
@test testUnitcell(getUnitcellSquareOctagon(4), 2, 2)
@test testUnitcell(getUnitcellSquareOctagon(2), 2, 2)
@test testUnitcell(getUnitcellSquareOctagon(5), 2, 2)
#       - getUnitcellBCC2D
@test testUnitcell(getUnitcellBCC2D(), 2, 2)
@test testUnitcell(getUnitcellBCC2D(1), 2, 2)
#       - getUnitcellFullyConnectedSquare
@test testUnitcell(getUnitcellFullyConnectedSquare(), 2, 2)
@test testUnitcell(getUnitcellFullyConnectedSquare(1), 2, 2)


#   - TRIANGULAR
#       - getUnitcellTriangular
@test testUnitcell(getUnitcellTriangular(), 2, 2)
@test testUnitcell(getUnitcellTriangular(1), 2, 2)
@test testUnitcell(getUnitcellTriangular(3), 2, 2)
#       - getUnitcellHoneycomb
@test testUnitcell(getUnitcellHoneycomb(), 2, 2)
@test testUnitcell(getUnitcellHoneycomb(1), 2, 2)
@test testUnitcell(getUnitcellHoneycomb(4), 2, 2)
@test testUnitcell(getUnitcellHoneycomb(2), 2, 2)
@test testUnitcell(getUnitcellHoneycomb(5), 2, 2)
@test testUnitcell(getUnitcellHoneycomb(3), 2, 2)
#       - getUnitcellKagome
@test testUnitcell(getUnitcellKagome(), 2, 2)
@test testUnitcell(getUnitcellKagome(1), 2, 2)
#       - getUnitcellKagomeMinus
@test testUnitcell(getUnitcellKagomeMinus(), 2, 2)
@test testUnitcell(getUnitcellKagomeMinus(1), 2, 2)
@test testUnitcell(getUnitcellKagomeMinus(4), 2, 2)
#       - getUnitcellHoneycombXXX
@test testUnitcell(getUnitcellHoneycombXXX(), 2, 2)
@test testUnitcell(getUnitcellHoneycombXXX(1), 2, 2)
@test testUnitcell(getUnitcellHoneycombXXX(2), 2, 2)





#   3D UNITCELLS
#   - CUBIC / FCC
#       - getUnitcellDiamond
@test testUnitcell(getUnitcellDiamond(), 3, 3)
@test testUnitcell(getUnitcellDiamond(1), 3, 3)
@test testUnitcell(getUnitcellDiamond(2), 3, 3)
@test testUnitcell(getUnitcellDiamond(3), 3, 3)
#       - getUnitcellBCC
@test testUnitcell(getUnitcellBCC(), 3, 3)
@test testUnitcell(getUnitcellBCC(1), 3, 3)
#       - getUnitcellPyrochlore
@test testUnitcell(getUnitcellPyrochlore(), 3, 3)
@test testUnitcell(getUnitcellPyrochlore(1), 3, 3)

#   - (X,3)y FAMILY
#       - getUnitcell_8_3_a
@test testUnitcell(getUnitcell_8_3_a(), 3, 3)
@test testUnitcell(getUnitcell_8_3_a(1), 3, 3)
@test testUnitcell(getUnitcell_8_3_a(4), 3, 3)
#       - getUnitcell_8_3_b
@test testUnitcell(getUnitcell_8_3_b(), 3, 3)
@test testUnitcell(getUnitcell_8_3_b(1), 3, 3)
@test testUnitcell(getUnitcell_8_3_b(4), 3, 3)
#       - getUnitcell_8_3_c
@test testUnitcell(getUnitcell_8_3_c(), 3, 3)
@test testUnitcell(getUnitcell_8_3_c(1), 3, 3)
@test testUnitcell(getUnitcell_8_3_c(4), 3, 3)
#       - getUnitcell_8_3_n
@test testUnitcell(getUnitcell_8_3_n(), 3, 3)
@test testUnitcell(getUnitcell_8_3_n(1), 3, 3)
@test testUnitcell(getUnitcell_8_3_n(4), 3, 3)
#       - getUnitcell_9_3_a
@test testUnitcell(getUnitcell_9_3_a(), 3, 3)
@test testUnitcell(getUnitcell_9_3_a(1), 3, 3)
@test testUnitcell(getUnitcell_9_3_a(4), 3, 3)
@test testUnitcell(getUnitcell_9_3_a(2), 3, 3)
@test testUnitcell(getUnitcell_9_3_a(5), 3, 3)
#       - getUnitcell_10_3_a / getUnitcellHyperoctagon
@test testUnitcell(getUnitcell_10_3_a(), 3, 3)
@test testUnitcell(getUnitcell_10_3_a(1), 3, 3)
@test testUnitcell(getUnitcell_10_3_a(4), 3, 3)
@test testUnitcell(getUnitcellHyperoctagon(), 3, 3)
@test testUnitcell(getUnitcellHyperoctagon(1), 3, 3)
@test testUnitcell(getUnitcellHyperoctagon(4), 3, 3)
#       - getUnitcell_10_3_b / getUnitcellHyperhoneycomb
@test testUnitcell(getUnitcell_10_3_b(), 3, 3)
@test testUnitcell(getUnitcell_10_3_b(1), 3, 3)
@test testUnitcell(getUnitcell_10_3_b(4), 3, 3)
@test testUnitcell(getUnitcell_10_3_b(2), 3, 3)
@test testUnitcell(getUnitcell_10_3_b(5), 3, 3)
@test testUnitcell(getUnitcellHyperhoneycomb(), 3, 3)
@test testUnitcell(getUnitcellHyperhoneycomb(1), 3, 3)
@test testUnitcell(getUnitcellHyperhoneycomb(4), 3, 3)
@test testUnitcell(getUnitcellHyperhoneycomb(2), 3, 3)
@test testUnitcell(getUnitcellHyperhoneycomb(5), 3, 3)
#       - getUnitcell_10_3_c
@test testUnitcell(getUnitcell_10_3_c(), 3, 3)
@test testUnitcell(getUnitcell_10_3_c(1), 3, 3)
@test testUnitcell(getUnitcell_10_3_c(4), 3, 3)
#       - getUnitcell_10_3_d
@test testUnitcell(getUnitcell_10_3_d(), 3, 3)
@test testUnitcell(getUnitcell_10_3_d(1), 3, 3)
@test testUnitcell(getUnitcell_10_3_d(4), 3, 3)
