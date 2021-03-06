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

function generatebdg(model::UnitCellQModel)
    typeof(model.qtype) <: Fermion ||
        error("BdG Hamiltonian only for FermionType models!")

    ft, f = fermionoperators(model.qtype)
    n_sites = prod(model.lattice.sizes) * model.lattice.unitc.n
    H = zeros(eltype(model.inters[1]), n_sites, n_sites)
    for is in Iterators.product([1:l for l in model.lattice.sizes]...)
        for interaction in model.inters
            if support(interaction) > 2
                println("dropping higher order term")
                continue
            end
            ns = interaction.ucidxs
            offs = interaction.offsets
            indexes = zeros(Int, support(interaction))
            crosses = []
            isinside = true
            for i in 1:support(interaction)
                index, crossings = sitelinearindex(model.lattice, ns[i], offs[i] .+ is)
                if isnothing(index)
                    isinside = false
                    break
                end
                indexes[i] = index
                push!(crosses, crossings)
            end
            # only two-fermionic terms
            if isinside && length(indexes) < 3
                sign = +1
                for n = 1:dimension(model)
                    if model.lattice.bcs[n] == :APBC
                        #println("$indexes, $crosses, $(interaction.amp)")
                        if isodd(crosses[2][n] - crosses[1][n])
                            sign *= -1
                        end
                    end
                end
                for term in interaction.terms
                    if term == (ft, f)
                        i, j = indexes
                        H[i, j] += sign * interaction.amp
                    elseif term == (f, ft)
                        j, i  = indexes
                        H[i, j] += sign * conj(interaction.amp)
                    else
                        error("term not supported yet!")
                    end
                end
            end
        end
    end
    H
end
