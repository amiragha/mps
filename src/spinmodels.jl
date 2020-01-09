# examples UnitCells
chainunitcell = UnitCell{1}(1, [0.], +1.)
squareunitcell = UnitCell{2}(1, [(0, 0)], [(0, 1), (1, 0)])
triangularunitcell = UnitCell{2}(1, [(0, 0)], [(1, 0), (0.5, sin(pi/3))])
kagomestripunitcell = UnitCell{1}(3, [0., 0., 0.], +1.)
honeycombunitcell = UnitCell{2}(2, [(0, 0), (0, 1)], [(sqrt(3), 0), (sqrt(3)/2, 3/2)])

# example types
spinhalf = SpinType(2)

sz, sp, sm = spinoperators(1/2)

# example models
function j1j2model(lx::Int, j1::Float64, j2::Float64;
                   boundary::Symbol=:OBC,
                   symmetry::Type{<:AbstractCharge}=Trivial)

    heis = nbodyopexpansion(2,
                            0.25 * (ringexchangeoperator(2) - I(4)),
                            mode=:SYMBOL)

    heis1 = QModelInteraction{1, 2, Float64}(
        j1,
        (1, 1),
        ((0,), (1,)), heis)
    heis2 = QModelInteraction{1, 2, Float64}(
        j2,
        (1, 1),
        ((0,), (2,)), heis)
    terms = QModelInteraction[]
    j1 != 0 && push!(terms, heis1)
    j2 != 0 && push!(terms, heis2)

    UnitCellQModel(spinhalf,
                   QLattice(chainunitcell, lx, boundary),
                   symmetry,
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
                             symmetry::Type{<:AbstractCharge}=Trivial)

    heis = nbodyopexpansion(2,
                            0.25 * (ringexchangeoperator(2) - I(4)),
                            mode=:SYMBOL)
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

    R4 = nbodyopexpansion(4,
                          ringexchangeoperator(4) - 0.25 * I(16),
                          mode=:SYMBOL)

    ring1 = QModelInteraction{2, 4, Float64}(
        k1,
        (1, 1, 1, 1),
        ((0, 0), (0, 1), (1, 1), (1, 0)),
        R4)

    ring2 = QModelInteraction{2, 4, Float64}(
        k2,
        (1, 1, 1, 1),
        ((0, 0), (0, 1), (1, 0), (1, -1)),
        R4)

    ring3 = QModelInteraction{2, 4, Float64}(
        k3,
        (1, 1, 1, 1),
        ((0, 0), (1, 0), (0, 1), (-1, 1)),
        R4)

    terms = QModelInteraction[]
    j1 != 0 && push!(terms, heis1)
    j2 != 0 && push!(terms, heis2)
    j3 != 0 && push!(terms, heis3)
    k1 != 0 && push!(terms, ring1)
    k2 != 0 && push!(terms, ring2)
    k3 != 0 && push!(terms, ring3)

    UnitCellQModel(spinhalf,
                   QLattice(triangularunitcell, (ly, lx), boundary),
                   symmetry,
                   terms)
end
