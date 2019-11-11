# examples UnitCells
chainunitcell = UnitCell{1}(1, [0.], +1.)
squareunitcell = UnitCell{2}(1, [(0, 0)], [(0, 1), (1, 0)])
triangularunitcell = UnitCell{2}(1, [(0, 0)], [(1, 0), (0.5, sin(pi/3))])

# example types
spinhalf = SpinType(2)

sz, sp, sm = Matrix(sz_half), Matrix(sp_half), Matrix(sm_half)

# example models
function j1j2model(lx::Int, j1::Float64, j2::Float64)
    heis = nbodyopexpansion(2, 0.25 * (ringexchangeoperator(2) - I(4)))
    heis1 = QModelInteraction{1, 2, Float64}(
        j1,
        (1, 1),
        ((0,), (1,)), heis)
    heis2 = QModelInteraction{1, 2, Float64}(
        j2,
        (1, 1),
        ((0,), (2,)), heis)

    UnitCellQModel{SpinType, 1}(spinhalf,
                                QLattice(chainunitcell, lx, :OBC),
                                [heis1, heis2])
end

# triangular
# ls is (ly, lx)
function triangularmodel(ls::Tuple{Int, Int},
                         j1::Float64,
                         j2::Float64,
                         j3::Float64,
                         k1::Float64,
                         k2::Float64,
                         k3::Float64)

    heis = nbodyopexpansion(0.25 * (ringexchangeoperator(2) - I(4)))
    ly, lx = ls
    heis1 = QModelInteraction{1, 2, Float64}(
        j1,
        (1, 1),
        ((0, 0), (0, 1)),
        [(sz, sz), (sp, 0.5*sm), (sm, 0.5*sp)])

    heis2 = QModelInteraction{1, 2, Float64}(
        j2,
        (1, 1),
        ((0, 0), (1, 0)),
        [(sz, sz), (sp, 0.5*sm), (sm, 0.5*sp)])

    heis3 = QModelInteraction{1, 2, Float64}(
        j3,
        (1, 1),
        ((0, 0), (1, -1)),
        [(sz, sz), (sp, 0.5*sm), (sm, 0.5*sp)])

    R4 = nbodyopexpansion(ringexchangeoperator(4) - 0.25 * I(16))

    ring1 = QModelInteraction{1, 2, Float64}(
        k1,
        (1, 1, 1, 1),
        ((0, 0), (0, 1), (1, 1), (1, 0)),
        R4)

    UnitCellQModel{SpinType, 2}(spinhalf,
                                Qlattice(triangularunitcell, (ly, lx), :OBC),
                                [heis1, heis2, heis3, ring1])
end
