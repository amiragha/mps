"""
    fishman2mps(fsmset, maxdim)

Construct a matrix product state with a maximum bond dimension
`max_dim` from a set of Fishman gates `fsmset`.

Since the unitary gates orthogonalize the unitary matrix, if they are
applied in revese order to the configuration state (which is the state
in the occupation basis) they will generate the state in the original
spatial basis. So, a classical MPS from the configurations is made and
then the gates are applied in reverse order.

"""
function fishman2mps(fsmset::FishmanGateSet,
                     maxdim::Int64;
                     symmetry::Symbol=:NONE)

    if symmetry == :NONE
        mps = MatrixProductState{ComplexF64}(fsmset.lx, 2, fsmset.initconf)
    elseif symmetry == :Z2
        mps = MatrixProductState{ComplexF64}(fsmset.lx, 2, fsmset.initconf)
    end

    twobodyugate = Matrix{ComplexF64}(I, 4, 4)

    ## NOTE: the gates are applies in reverse order
    for n=length(fsmset.positions):-1:1
        site = fsmset.positions[n]
        θ = fsmset.θs[n]
        twobodyugate[2:3, 2:3] =
            [
                cos(θ) sin(θ);
                -sin(θ) cos(θ)
            ]

        ugatetensor = reshape(twobodyugate, 2, 2, 2, 2)

        # produce the counterclockwise convention starting from bottom left
        operator = permutedims(ugatetensor, [1,2,4,3])

        ## TODO: choose a better order of site or site+1 and push_to!
        move_center!(mps, site)
        apply_2siteoperator!(mps, site, operator, maxdim=maxdim, pushto=:R)
    end
    return mps
end
