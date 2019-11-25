"""

The structure containing the set of neaerest-neighbor gates when
applied on the initial configuration `initconf` generates the
approximate MPS representation for a gaussian fermionic model. `lx` is
the size of the system. The left acting site of each gate is given in
the `positions` vector and each gate can be constructed using the
`θs` or possibly `ϕs`.

"""
struct FishmanGateSet
    lx :: Int64
    initconf :: Vector{Int64}
    positions :: Vector{Int64}
    θs :: Vector{Float64}

    function FishmanGateSet(lx::Int64,
                            initconf::Vector{Int64},
                            positions::Vector{Int64},
                            θs::Vector{Float64})
        @assert length(initconf) == lx
        @assert length(positions) == length(θs)

        new(lx, initconf, positions, θs)
    end
end

"""
    generate_fishmangates(corr_matrix, threshold, verbose)

For a given two-point correlation matrix `corr_matrix` called lambda
perform the approach explained in arxiv:1504.07701 and generate the
local (nearest-neighbor) gates.

"""

function generate_fishmangates(corr_matrix::Matrix{T},
                               threshold::Float64=1.e-8;
                               alternate::Bool=false) where {T<:RLorCX}
    T <: Complex && error("correlation matrix is complex!")
    lx = size(corr_matrix)[1]

    #initconf = Vector{Int64}(lx)
    initconf = zeros(Int64, lx)
    positions = Int64[]
    θs = Float64[]

    expected_next_evalue = 1
    for site=1:lx
        block_end = site
        evalue = corr_matrix[site, site]
        corr_block = [evalue]
        vals = corr_block

        if alternate
            delta = abs(expeted_next_evalue - evalue)
        else
            delta = min(abs(evalue), abs(1-evalue))
        end

        while (delta > threshold && block_end < lx)
            block_end += 1

            corr_block = corr_matrix[site:block_end, site:block_end]
            vals = eigvals(corr_block)

            if alternate
                delta = (expected_next_evalue == 0) ? abs.(minimum(vals)) : abs.(1-maximum(vals))
            else
                delta = min(abs.(minimum(vals)), abs.(1-maximum(vals)))
            end
        end

        if block_end > site
            block_size = block_end - site + 1
            block_range = site:block_end

            if alternate
                v_index = (expected_next_evalue == 0) ? indmin(abs.(vals)) : argmin(abs.(1-vals))
            else
                v_index = argmax(abs.(vals .- 0.5))
            end
            evalue = vals[v_index]
            v = eigvecs(corr_block)[:,v_index]

            # Find the set of (block_size-1) unitary gate that
            # diagonalize the current correlation block
            ugate_block = Matrix{T}(I, block_size, block_size)

            for p = block_size-1:-1:1 # p here stands for "pivot"

                θ = atan(v[p + 1]/ v[p])
                ugate = Matrix{Float64}(1.0I, block_size, block_size)
                ugate[p:p+1, p:p+1] =
                    [
                        cos(θ) -sin(θ);
                        sin(θ)  cos(θ)
                    ]
                push!(positions, site+p-1)
                push!(θs, θ)
                ## QQQ? acting on the vectors from left! Explain, why
                ## is this a better choice?
                v = transpose(transpose(v) * ugate)

                ugate_block = ugate_block * ugate
            end

            ugate_extended = Matrix{T}(I, lx, lx)
            ugate_extended[block_range, block_range] = ugate_block
            corr_matrix= Symmetric(ugate_extended' * corr_matrix * ugate_extended)
            # display(corr_matrix)
        end
        initconf[site] = round(evalue)
        expected_next_evalue = (round(evalue) == 1) ? 0 : 1
    end
    # display(corr_matrix)
    FishmanGateSet(lx, initconf, positions, θs)
end
