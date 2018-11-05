# CREATE LATTICES FROM UNITCELLS
# PUT TOGEHTER PATTERNS OF UNITCELLS PERIODICALLY


# 2d
function getLatticePeriodic(
        unitcell    :: U,
        extent      :: NTuple{2,Int64},
        lattice_type:: Type{LA}
    ) :: LA where {D,L,S<:AbstractSite{L,D},B<:AbstractBond{L,2},U<:AbstractUnitcell{D,2,L,S,B},LA<:AbstractLattice{D,2,L,S,B}}

end
