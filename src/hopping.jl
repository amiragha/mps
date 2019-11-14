function hopping_chain(lx         ::Int,
                       t          ::T=1.0,
                       mu         ::T=0.0;
                       boundary   ::Symbol=:OBC) where {T<:RLorCX}

    hopmatrix = diagm(0 => mu .* ones(T, lx)) +
        diagm(1 => -t .* ones(T, lx-1)) +
        diagm(-1 => conj(-t) .* ones(T, lx-1))

    if boundary == :OBC
        return hopmatrix
    elseif boundary == :PBC
        hopmatrix[1, lx] = -t
        hopmatrix[lx, 1] = conj(-t)
        return hopmatrix
    elseif boundary == :APBC
        hopmatrix[1, lx] = t
        hopmatrix[lx, 1] = conj(t)
        return hopmatrix
    else
        error("unrecognized boundary condition : ", boundary)
    end
    hopmatrix
end

function nnhoppingchain(lx::Int, t1::T, t2::T, mu::T=0.0;
                        boundary::Symbol=:OBC) where {T<:RLorCX}

    hopmatrix = diagm(0 => mu .* ones(T, lx)) +
        diagm(1 => -t1 .* ones(T, lx-1)) +
        diagm(-1 => conj(-t1) .* ones(T, lx-1)) +
        diagm(2 => -t2 .* ones(T, lx-2)) +
        diagm(-2 => conj(-t2) .* ones(T, lx-2))

    if boundary == :OBC
        return hopmatrix
    elseif boundary == :PBC
        hopmatrix[1, lx] = -t1
        hopmatrix[1, lx-1] = -t2
        hopmatrix[2, lx] = -t2
        hopmatrix[lx, 1] = conj(-t1)
        hopmatrix[lx-1, 1] = conj(-t2)
        hopmatrix[lx, 2] = conj(-t2)
        return hopmatrix

    elseif boundary == :APBC
        hopmatrix[1, lx] = t1
        hopmatrix[1, lx-1] = t2
        hopmatrix[2, lx] = t2
        hopmatrix[lx, 1] = conj(t1)
        hopmatrix[lx-1, 1] = conj(t2)
        hopmatrix[lx, 2] = conj(t2)
        return hopmatrix
    else
        error("unrecognized boundary condition : ", boundary)
    end
    hopmatrix
end

function generatebdg(model::UnitCellQFermionModel)
    typeof(model.qtype) == FermionType ||
        error("BdG Hamiltonian only for FermionType models!")
    model.lattice.bc == :OBC || error("Only OBC for now!")

    n_sites = prod(mode.lattice.sizes)
    H = zeros(n_sites, n_sites)
    for is in Iterators.product([1:l for l in model.lattice.sizes]...)
        for interaction in model.inters
            ns = interaction.ucidxs
            offs = interaction.offsets
            indexes = [sitelinearindex(model.lattice, ns[i], offs[i] .+ is)
                       for i in 1:support(interaction)]
            if all(indexes .> 0)
                if interaction.inters == FermionHop
                    i, j = indexes
                    H[i, j] += interaction.amp
                    H[j, i] += interaction.amp
                elseif interaction.inters == FermionNum
                    index = indexes[1]
                    H[index, index] += interaction.amp
                else
                    error("unsupported fermion interaction!")
                end
            end
        end
    end
    H
end
