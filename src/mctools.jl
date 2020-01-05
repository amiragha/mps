"""
    random_ones(n, m)

returns a Vector{Int} of size `n` with `m` randomly chosen 1s and
the rest 0s.  TODO: explain how it works?

"""
function random_ones(n::Int, m::Int)
    @assert n > 0 && m <= n && m > 0

    # TODO: implement the m > n/2 case with random for zeros!
    output = zeros(Int, n)

    for pivot = (n-m+1):n
        pick = rand(1:n-m+1)
        if output[pick] == 1
            # pick is already chosen, choose the pivot.
            output[pivot] = 1
        else
            # choose pick.
            output[pick] = 1
        end
    end
    return output
end
