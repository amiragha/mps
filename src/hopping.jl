function hopping_chain(lx         ::Int64,
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
end
