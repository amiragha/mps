"""
    fishman2mps(fsmset, maxdim)

Construct a matrix product state with a maximum bond dimension
`max_dim` from a set of Fishman gates `fsmset`.

Since the unitary gates orthogonalize the correlation matrix, if they
are applied in revese order to the configuration state (which is the
state in the occupation basis) they will generate the state in the
original spatial basis. So, a classical MPS from the configurations is
made and then the gates are applied in reverse order.

"""

##TODO: think about how to incorporate the complex gates (that may
##come when the hoppings are complex numbers) into this type nicely!
# function fishman2mps(fsmset::FishmanComplexGateSet,
#                      maxdim::Int64;
#                      symmetry::Symbol=:NONE)

#     ##TODO: why are these defined complex here?!
#     if symmetry == :NONE
#         mps = MatrixProductState{Float64}(fsmset.lx, 2, fsmset.initconf)
#         _applyfishmangates!(mps, fsmset, maxdim)

#     elseif symmetry == :U1
#         mps = SymMatrixProductState{Float64}(fsmset.lx, 2, fsmset.initconf)
#         _applyfishmangates!(mps, fsmset, maxdim)
#     else
#         error( "Symmetry $symmetry is not defined!")
#     end
#     mps
# end

function fishman2mps(fsmset::FishmanGateSet,
                     maxdim::Int64;
                     symmetry::Symbol=:NONE)

    ##TODO: why are these defined complex here?!
    if symmetry == :NONE
        mps = MatrixProductState{Float64}(fsmset.lx, 2, fsmset.initconf)
        _applyfishmangates!(mps, fsmset, maxdim)

    elseif symmetry == :U1
        mps = SymMatrixProductState{Float64}(fsmset.lx, 2, fsmset.initconf)
        _applyfishmangates!(mps, fsmset, maxdim)
    else
        error( "Symmetry $symmetry is not defined!")
    end
    mps
end

function _applyfishmangates!(mps::MatrixProductState{T},
                             fsmset::FishmanGateSet,
                             maxdim::Int) where{T<:RLorCX}
    twobodyugate = Matrix{T}(I, 4, 4)
    ## NOTE: the gates are applies in reverse order
    for n=length(fsmset.positions):-1:1
        site = fsmset.positions[n]
        θ = fsmset.θs[n]
        twobodyugate[2:3, 2:3] = rotationmat(θ)

        ugatetensor = reshape(twobodyugate, 2, 2, 2, 2)

        # produce the counterclockwise convention starting from bottom left
        ##NOTE: this is needed just because I followed the convention,
        ##probably best to avoid id, but need to remain consistant
        operator = permutedims(ugatetensor, [1,2,4,3])

        ## TODO: choose a better order of site or site+1 and push_to!
        move_center!(mps, site)
        apply_2siteoperator!(mps, site, operator, maxdim=maxdim, pushto=:R, normalizeS=true)
    end
    nothing
end

function _applyfishmangates!(mps::SymMatrixProductState{Tv},
                             fsmset::FishmanGateSet,
                             maxdim::Int) where{Tv}

    # u is the two-body U gate matrix
    u = eye(Float64, 0, [0,1,2], [1,2,1])

    #   twobodyugate = Matrix{ComplexF64}(I, 4, 4)
    ## NOTE: the gates are applies in reverse order
    rlegs = (STLeg(+1, [0,1], [1,1]), STLeg(+1, [0,1], [1,1]))
    clegs = (STLeg(-1, [0,1], [1,1]), STLeg(-1, [0,1], [1,1]))

    for n=length(fsmset.positions):-1:1
        site = fsmset.positions[n]
        θ = fsmset.θs[n]
        change_nzblk!(u, (1,1), rotationmat(θ))

        uten = defuse_leg(defuse_leg(u, 2, clegs), 1, rlegs)

        ## TODO: choose a better order of site or site+1 and push_to!
        move_center!(mps, site)
        apply_2siteoperator!(mps, site, uten, maxdim=maxdim, pushto=:R, normalizeS=true)
    end
    nothing
end

function rotationmat(θ::Float64)
    R = [
        cos(θ) sin(θ);
        -sin(θ) cos(θ)
    ]
    R
end
