"""
    gmps2mps(gmps, maxdim)

Construct a matrix product state with a maximum bond dimension
`max_dim` from a set of GaussianMPS gates `gmps`.

Since the unitary gates orthogonalize the correlation matrix, if they
are applied in revese order to the configuration state (which is the
state in the occupation basis) they will generate the state in the
original spatial basis. So, a classical MPS from the configurations is
made and then the gates are applied in reverse order.

"""
function gmps2mps(gmps::GaussianMPS,
                  maxdim::Int64;
                  symmetry::Symbol=:NONE)

    ##TODO: why are these defined complex here?!
    if symmetry == :NONE
        mps = MPState{Float64}(gmps.lx, 2, gmps.initconf)
        _applygmpsgates!(mps, gmps, maxdim)

    elseif symmetry == :U1
        mps = U1MPState(Float64, gmps.lx, 2, gmps.initconf)
        _applygmpsgates!(mps, gmps, maxdim)
    else
        error( "Symmetry $symmetry is not defined!")
    end
    mps
end

function _applygmpsgates!(mps::MPState{T},
                          gmps::GaussianMPS,
                          maxdim::Int) where {T}
    twobodyugate = Matrix{T}(I, 4, 4)
    ## NOTE: the gates are applies in reverse order
    for n=length(gmps.xs):-1:1
        site = gmps.xs[n]
        θ = gmps.θs[n]
        twobodyugate[2:3, 2:3] = rotationmat(θ)

        ugatetensor = reshape(twobodyugate, 2, 2, 2, 2)

        # produce the counterclockwise convention starting from bottom left
        ##NOTE: this is needed just because I followed the convention,
        ##probably best to avoid id, but need to remain consistant
        operator = permutedims(ugatetensor, [1,2,4,3])

        ## TODO: choose a better order of site or site+1 and push_to!
        center_at!(mps, site)
        apply!(mps, operator, site, maxdim=maxdim, pushto=:R, svnormalize=true)
    end
    nothing
end

function _applygmpsgates!(mps::MPState{U1,T},
                          gmps::GaussianMPS,
                          maxdim::Int) where{T}

    Vd = U1Space(0=>1, 1=>1)
    # u is the two-body U gate matrix
    u = eye(Float64, fuse(false, Vd,Vd))

    #   twobodyugate = Matrix{ComplexF64}(I, 4, 4)
    ## NOTE: the gates are applies in reverse order

    for n=length(gmps.xs):-1:1
        site = gmps.xs[n]
        θ = gmps.θs[n]
        u[Sector{U1}(1,1)] =  rotationmat(θ)

        uten = splitleg(splitleg(u, 2, (dual(Vd), dual(Vd))), 1, (Vd,Vd))

        ## TODO: choose a better order of site or site+1 and push_to!
        center_at!(mps, site)
        apply!(mps, uten, site, maxdim=maxdim, pushto=:R, svnormalize=true)
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
