using Base.LinAlg.ARPACK

import Base.LinAlg: BlasInt

"""
    eigsfn(matvecA!, v0, sym; nev=6, ncv=max(20,2*nev+1), which=:LM, tol=0.0, maxiter=300, ritzvec=true)

This function uses the wrappers for the _aupd that is the
implementation of Implicitly restarted Arnoldi method (IRAM) and _eupd
for postprocessing. This is very similar functionality to usual
`LinAlg.eigs` except that instead of a matrix (linear map) an actual
function `matvecA!` is given.

Note that here the initial vector `v0` and bool 'sym' should always be
given to the function because parameters are chosen based on them.

Note that shift-and-invert method doesn't work here (because we don't
have the inverse or some useful factorization of the linear map), so
only `:LM, :LR, :Sr` are supported.

"""
function eigsfn(matvecA!::Function, v0::Vector{T}, sym::Bool;
                nev::Integer=6, ncv::Integer=max(20,2*nev+1), which::Symbol=:LM,
                tol=0.0, maxiter::Integer=300,
                ritzvec::Bool=true) where {T<:Union{Float64, Complex128}}

    n = length(v0)
    iscmplx = T <: Complex

    nevmax=sym ? n-1 : n-2
    if nevmax <= 0
        throw(ArgumentError("input matrix A is too small. Use eigfact instead."))
    end
    if nev > nevmax
        warn("Adjusting nev from $nev to $nevmax")
        nev = nevmax
    end
    if nev <= 0
        throw(ArgumentError("requested number of eigenvalues (nev) must be ≥ 1, got $nev"))
    end
    ncvmin = nev + (sym ? 1 : 2)
    if ncv < ncvmin
        warn("Adjusting ncv from $ncv to $ncvmin")
        ncv = ncvmin
    end
    ncv = BlasInt(min(ncv, n))

    # This has to be "I" because "G" for generalized is not supported
    # for this function (bc it requires factorization)
    bmat = "I"

    if (which != :LM && which != :LR && which != :SR)
        throw(ArgumentError("which must be :LM, :LR, :SR, got $(repr(which))"))
    end

    whichstr = "LM"
    if which == :LR
        whichstr = (!sym ? "LR" : "LA")
    end
    if which == :SR
        whichstr = (!sym ? "SR" : "SA")
    end

    matvecB = x -> x
    # regular mode (no shift-invert) or inverse
    mode       = 1
    solveSI = x -> x
    sigma = zero(T)

    # Compute the Ritz values and Ritz vectors
    (resid, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, TOL) =
        ARPACK.aupd_wrapper(T, matvecA!, matvecB, solveSI, n, sym, iscmplx, bmat, nev, ncv, whichstr, tol, maxiter, mode, v0)

    # Postprocessing to get eigenvalues and eigenvectors
    output = ARPACK.eupd_wrapper(T, n, sym, iscmplx, bmat, nev, whichstr, ritzvec, TOL,
                                 resid, ncv, v, ldv, sigma, iparam, ipntr, workd, workl, lworkl, rwork)

    nev = length(output[1])
    nconv = output[ritzvec ? 3 : 2]
    nev ≤ nconv || warn("not all wanted Ritz pairs converged. Requested: $nev, converged: $nconv")

    return output
end
