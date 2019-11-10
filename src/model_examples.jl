# examples UnitCells
chainunitcell = UnitCell{1}(1, [0.], +1.)
squareunitcell = UnitCell{2}(1, [(0, 0)], [(0, 1), (1, 0)])
triangularunitcell = UnitCell{2}(1, [(0, 0)], [(1, 0), (0.5, sin(pi/3))])

# example types
spinhalf = SpinType(2)

# example models
function j1j2model(j1::Float64, j2::Float64, lx::Int)
    heisj1 = QModelTerm{1, 2, Float64}(
        j1,
        (1, 1),
        ((0,), (1,)),
        (1,),
        [(sz, sz), (sp, 0.5*sm), (sm, 0.5*sp)])
    heisj2 = QModelTerm{1, 2, Float64}(
        j2,
        (1, 1),
        ((0,), (2,)),
        (1,),
        [(sz, sz), (sp, 0.5*sm), (sm, 0.5*sp)])

    QModel(spinhalf,
           QLattice(chainunitcell, lx),
           [heisj1, heisj2])
end
