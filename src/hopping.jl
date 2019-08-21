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
        diagm(2 => -t2 .* ones(T, lx-1)) +
        diagm(-2 => conj(-t2) .* ones(T, lx-1))

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
