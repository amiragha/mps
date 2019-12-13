function t1t2model(lx::Int,
                   t1::T, t2::T, mu::Float64;
                   boundary::Symbol=:OBC) where {T<:Number}

    hop1 = FermionQModelInteraction{1, 2, T}(
        t1,
        (1, 1),
        ((0,), (1,)), [FermionHop])
    hop2 = FermionQModelInteraction{1, 2, T}(
        t2,
        (1, 1),
        ((0,), (2,)), [FermionHop])

    chemical = FermionQModelInteraction{1, 1, Float64}(
        mu, (1,), (0,), [FermionNum])

    terms = FermionQModelInteraction[]
    t1 != 0 && push!(terms, hop1)
    t2 != 0 && push!(terms, hop2)

    UnitCellQModel(fermion,
                   QLattice(triangularunitcell, lx, boundary),
                   terms)
end

function triangularhopping(ls::Tuple{Int, Int},
                           t1::T,
                           t2::T,
                           t3::T,
                           mu::Float64;
                           boundary::Tuple{Symbol, Symbol}=(:PBC, :OBC)) where{T}

    ly, lx = ls

    fermion = Fermion{:f}()
    ft, f = fermionoperators(fermion)

    fhop = [(ft, f), (f, ft)]
    hop1 = FermionQModelInteraction{2, 2, T}(
        t1,
        (1, 1),
        ((0, 0), (0, 1)), fhop)

    hop2 = FermionQModelInteraction{2, 2, T}(
        t2,
        (1, 1),
        ((0, 0), (1, 0)), fhop)

    hop3 = FermionQModelInteraction{2, 2, T}(
        t3,
        (1, 1),
        ((0, 0), (1, -1)), fhop)

    terms = FermionQModelInteraction[]
    t1 != 0 && push!(terms, hop1)
    t2 != 0 && push!(terms, hop2)
    t3 != 0 && push!(terms, hop3)

    UnitCellQModel(fermion,
                   QLattice(triangularunitcell, (ly, lx), boundary),
                   terms)
end
