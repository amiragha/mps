# examples UnitCells
chainunitcell = UnitCell{1}(1, [0.], +1.)
squareunitcell = UnitCell{2}(1, [(0, 0)], [(0, 1), (1, 0)])
triangularunitcell = UnitCell{2}(1, [(0, 0)], [(1, 0), (0.5, sin(pi/3))])

# example types
spinhalf = SpinType(2)

sz, sp, sm = spinoperators(1/2)

# example models
function j1j2model(lx::Int, j1::Float64, j2::Float64;
                   boundary::Symbol=:OBC,
                   symmetry::Symbol=:NONE)

    heis = nbodyopexpansion(2,
                            0.25 * (ringexchangeoperator(2) - I(4)),
                            symmetry=symmetry)

    if symmetry == :NONE
        qitype = QModelInteraction
    elseif symmetry == :U1
        qitype = SymQModelInteraction
    else
        error()
    end
    heis1 = qitype{1, 2, Float64}(
        j1,
        (1, 1),
        ((0,), (1,)), heis)
    heis2 = qitype{1, 2, Float64}(
        j2,
        (1, 1),
        ((0,), (2,)), heis)
    terms = []
    j1 != 0 && push!(terms, heis1)
    j2 != 0 && push!(terms, heis2)

    UnitCellQModel{SpinType, 1}(spinhalf,
                                QLattice(chainunitcell, lx, boundary),
                                terms)
end

# triangular
# ls is (ly, lx)
function triangularspinmodel(ls::Tuple{Int, Int},
                             j1::Float64,
                             j2::Float64,
                             j3::Float64,
                             k1::Float64,
                             k2::Float64,
                             k3::Float64;
                             boundary::Tuple{Symbol,Symbol}=(:PBC, :OBC),
                             symmetry::Symbol=:NONE)

    heis = nbodyopexpansion(2,
                            0.25 * (ringexchangeoperator(2) - I(4)),
                            symmetry=symmetry)
    ly, lx = ls

    if symmetry == :NONE
        qitype = QModelInteraction
    elseif symmetry == :U1
        qitype = SymQModelInteraction
    else
        error()
    end

    heis1 = qitype{2, 2, Float64}(
        j1,
        (1, 1),
        ((0, 0), (0, 1)), heis)

    heis2 = qitype{2, 2, Float64}(
        j2,
        (1, 1),
        ((0, 0), (1, 0)), heis)

    heis3 = qitype{2, 2, Float64}(
        j3,
        (1, 1),
        ((0, 0), (1, -1)), heis)

    R4 = nbodyopexpansion(4,
                          ringexchangeoperator(4) - 0.25 * I(16),
                          symmetry=symmetry)

    ring1 = qitype{2, 4, Float64}(
        k1,
        (1, 1, 1, 1),
        ((0, 0), (0, 1), (1, 1), (1, 0)),
        R4)

    ring2 = qitype{2, 4, Float64}(
        k2,
        (1, 1, 1, 1),
        ((0, 0), (0, 1), (1, 0), (1, -1)),
        R4)

    ring3 = qitype{2, 4, Float64}(
        k3,
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
                                QLattice(triangularunitcell, (ly, lx), boundary),
                                terms)
end
