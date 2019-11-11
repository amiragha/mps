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

    terms = []
    j1 != 0 && push!(terms, heis1)
    j2 != 0 && push!(terms, heis2)

    UnitCellQModel{SpinType, 1}(spinhalf,
                                QLattice(chainunitcell, lx, :OBC),
                                [heis1, heis2])
end

# triangular
# ls is (ly, lx)
function triangularspinmodel(ls::Tuple{Int, Int},
                             j1::Float64,
                             j2::Float64,
                             j3::Float64,
                             k1::Float64,
                             k2::Float64,
                             k3::Float64)

    heis = nbodyopexpansion(2, 0.25 * (ringexchangeoperator(2) - I(4)))
    ly, lx = ls

    heis1 = QModelInteraction{2, 2, Float64}(
        j1,
        (1, 1),
        ((0, 0), (0, 1)), heis)

    heis2 = QModelInteraction{2, 2, Float64}(
        j2,
        (1, 1),
        ((0, 0), (1, 0)), heis)

    heis3 = QModelInteraction{2, 2, Float64}(
        j3,
        (1, 1),
        ((0, 0), (1, -1)), heis)

    R4 = nbodyopexpansion(4, ringexchangeoperator(4) - 0.25 * I(16))

    ring1 = QModelInteraction{2, 4, Float64}(
        k1,
        (1, 1, 1, 1),
        ((0, 0), (0, 1), (1, 1), (1, 0)),
        R4)

    ring2 = QModelInteraction{2, 4, Float64}(
        k1,
        (1, 1, 1, 1),
        ((0, 0), (0, 1), (1, 0), (1, -1)),
        R4)

    ring3 = QModelInteraction{2, 4, Float64}(
        k1,
        (1, 1, 1, 1),
        ((0, 0), (1, 0), (0, 1), (-1, 1)),
        R4)

    terms = []
    j1 != 0 && push!(terms, heis1)
    j2 != 0 && push!(terms, heis2)
    j3 != 0 && push!(terms, heis3)
    k1 != 0 && push!(terms, ring1)
    k2 != 0 && push!(terms, ring2)
    k3 != 0 && push!(terms, ring3)

    UnitCellQModel{SpinType, 2}(spinhalf,
                                QLattice(triangularunitcell, (ly, lx), :OBC),
                                terms)
end
