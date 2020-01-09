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
                   Z2Charge,
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
                   Z2Charge,
                   terms)
end

function honeycombhopping(ls::Tuple{Int, Int},
                          t1::T,
                          t2::T,
                          t3::T,
                          mu1::Float64=0.0,
                          mu2::Float64=0.0;
                          boundary::Tuple{Symbol, Symbol}=(:PBC, :OBC)) where{T}

    ly, lx = ls

    fermion = Fermion{:f}()
    ft, f = fermionoperators(fermion)

    fhop = [(ft, f), (f, ft)]
    hop1 = FermionQModelInteraction{2, 2, T}(
        t1,
        (1, 2),
        ((0, 0), (0, 0)), fhop)

    hop2 = FermionQModelInteraction{2, 2, T}(
        t2,
        (2, 1),
        ((0, 0), (1, 0)), fhop)

    hop3 = FermionQModelInteraction{2, 2, T}(
        t3,
        (2, 1),
        ((0, 0), (1, -1)), fhop)

    terms = FermionQModelInteraction[]
    t1 != 0 && push!(terms, hop1)
    t2 != 0 && push!(terms, hop2)
    t3 != 0 && push!(terms, hop3)

    UnitCellQModel(fermion,
                   QLattice(honeycombunitcell, (ly, lx), boundary),
                   Z2Charge,
                   terms)
end

function kagomestriphopping(lx::Int,
                            tl::T,
                            tc::T,
                            mu::Float64;
                            boundary::Symbol=:OBC) where{T}

    fermion = Fermion{:f}()
    ft, f = fermionoperators(fermion)

    fhop = [(ft, f), (f, ft)]
    hop1 = FermionQModelInteraction{1, 2, T}(
        tl,
        (1, 1),
        ((0,), (1,)), fhop)

    hop2 = FermionQModelInteraction{1, 2, T}(
        tl,
        (3, 3),
        ((0,), (1,)), fhop)

    hop3 = FermionQModelInteraction{1, 2, T}(
        tc,
        (1, 2),
        ((0,), (0,)), fhop)

    hop4 = FermionQModelInteraction{1, 2, T}(
        tc,
        (3, 2),
        ((0,), (0,)), fhop)

    hop5 = FermionQModelInteraction{1, 2, T}(
        tc,
        (1, 2),
        ((0,), (1,)), fhop)

    hop6 = FermionQModelInteraction{1, 2, T}(
        tc,
        (3, 2),
        ((0,), (1,)), fhop)


    terms = FermionQModelInteraction[]
    if tl != zero(T)
        push!(terms, hop1)
        push!(terms, hop2)
    end
    if tc != zero(T)
        push!(terms, hop3)
        push!(terms, hop4)
        push!(terms, hop5)
        push!(terms, hop6)
    end
    UnitCellQModel(fermion,
                   QLattice(kagomestripunitcell, lx, boundary),
                   Z2Charge,
                   terms)
end
