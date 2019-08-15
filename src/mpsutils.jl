"""
    svdtrunc(A [;maxdim, tol])

Performs the singular value decomposition and truncated the singular
values accoding to the tolerance `tol` or `maxdim` specified,
whichever is lower.

Return value is the tuple of truncated U, S, Vt. Here S is a diagonal
matrix instead of a vector.

"""
function svdtrunc(A::Matrix{T}; maxdim::Int=200, tol::Float64=1.e-8) where{T<:Number}
    fact = svd(A, full=false)
    n = min(maxdim, sum(fact.S .> fact.S[1]*tol))

    fact.U[:, 1:n], Diagonal(fact.S[1:n]), fact.Vt[1:n, :]
end

"""
    entorpy(spectrum [; alpha])

calculate the entropy of a vector of numbers `spectrum`. The numbers
are assumed to be probabilities which means they are positive. The
numbers will be normalized so that they add up to 1. If `alpha=1` (the
default) it calculates the usual Shannon (Von-Neumann) entropy and if
`alpha > 1` calculates the Renyi entorpy.

If the input is a vector of vectors, then the entropy of each
individual vector is calculated.

"""
function entropy(spectrum::Vector{T};
                 alpha::Int64=1) where {T<:Number}

    s = spectrum ./ sum(spectrum)
    if alpha == 1
        return - sum(s .* log2.(s))
    else
        return log2(sum(s.^alpha))/(1-alpha)
    end
end

function entropy(spectrums::Vector{Vector{T}};
                 alpha::Int64=1) where {T<:Number}
    result = Vector{T}(undef, length(spectrums))
    for i in eachindex(spectrums)
        result[n] =  entropy(spectrums[i], alpha=alpha)
    end
    return result
end
