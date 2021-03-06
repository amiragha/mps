"""
    zipandgutzwiller(mps1, mps2, mode)

given two MPS, with the same number of sites and two physical
dimensions perform gutzwiller projection of the two by making the
matrices using iterative SVD. The final result is an MPS with 2
physical dimension where the final matrices are the tensor product of
the given ones according to the `mode` which can be :B14 for bosonic
or :F23 for fermionic.

"""

function zipandgutzwiller!(mps1::MPState{Y},
                           mps2::MPState{Y};
                           mode::Symbol=:B14,
                           maxdim::Int64=200) where {Y<:Array}
    if mode != :B14
        error("only Bosinic ↑↑ and ↓↓ is possible for no-symmetry MPSs!")
    end
    @assert mps1.d == mps2.d == 2
    lx = mps1.lx
    @assert mps2.lx == lx

    center_at!(mps1, 1)
    center_at!(mps2, 1)

    dims = ones(Int64, lx+1)
    matrices = Array{T, 3}[]

    # Gutzwiller projector for mode :B14
    gutzp = zeros(T, 2,2,2)
    gutzp[1,1,1] = 1
    gutzp[2,2,2] = 1

    E = ones(T, 1,1,1)
    for site=1:lx-1
        mat1 = mps1.matrices[site]
        mat2 = mps2.matrices[site]

        diml1, dimr1 = size(mat1, 1), size(mat1, 3)
        diml2, dimr2 = size(mat2, 1), size(mat2, 3)
        dimle = size(E, 1)

        @tensor C[le, dg, r2, r1] := ((E[le, l2, l1] * mat1[l1,d1,r1]) *
                                      gutzp[dg,d1,d2]) * mat2[l2,d2, r2]

        fact = svd(reshape(C, dimle*2, dimr2*dimr1), full=false)
        S, n, ratio = MPStateTools.truncate(fact.S, maxdim=maxdim)
        dims[site+1] = n
        U = fact.U[:, 1:n]
        push!(matrices, reshape(U, dimle, 2, n))
        E = reshape(Diagonal(S) * fact.Vt[1:n, :], n, dimr2, dimr1)
    end
    mat1 = mps1.matrices[lx]
    mat2 = mps2.matrices[lx]
    @tensor C[le, dg, r2, r1] := ((E[le, l2, l1] * mat1[l1,d1,r1]) *
                                  gutzp[dg,d1,d2]) * mat2[l2,d2,r2]
    push!(matrices, reshape(C, size(E, 1), 2, 1))

    return MPState{T}(lx, 2, dims, matrices, lx)
end

function zipandgutzwiller!(mps1::MPState{Y},
                           mps2::MPState{Y};
                           mode::Symbol=:F23,
                           maxdim::Int64=200) where {Y<:SymTensor}
    if mode==:B14
        return _zipandgutzwiller_B14!(mps1, mps2, maxdim=maxdim)
    elseif mode==:F23
        return _zipandgutzwiller_F23!(mps1, mps2, maxdim=maxdim)
    else
        error("mode not defined : ", mode)
    end
end

#bosonic version
function _zipandgutzwiller_B14!(mps1::MPState{Y},
                                mps2::MPState{Y};
                                maxdim::Int64=200) where {Y}
    #@assert mps1.d == mps2.d == 2
    lx = length(mps1)
    @assert length(mps2) == lx

    center_at!(mps1, 1)
    center_at!(mps2, 1)

    mps = MPState{Y}()

    ## NOTE: in order to make the gutzwiller projector respect the U1
    ## symmetry we need to do the follwoing. Assume that ↑↑
    ## corresponds to ↑ or 2 and ↓↓ corresponds or ↓ or 0.
    G = fill(one(T), 0,
             (STLeg(+1, [0,2], [1,1]),
              STLeg(-1, [0,1], [1,1]),
              STLeg(-1, [0,1], [1,1])))

    E = fill(one(T), 0, (STLeg(+1, [0], [1]),
                         STLeg(-1, [0], [1]),
                         STLeg(-1, [0], [1])))
    for l=1:lx-1
        A = mps1.matrices[l]
        B = mps2.matrices[l]

        ##NOTE:
        # contracting (E2, A1)
        # First tensor has EA = E1 E3 A2 A3
        # contracting (A2, G2)
        # Scond tensor has EAG = E1 E3 G1 A3 G3
        # contracting (E3, B1) and (G3, B2)
        # Final tensor has E1 G1 A3 B3
        C = contract(contract(contract(E, (1,-1, 2), A, (-1, 3, 4)),
                              (1, 2, -1, 4), G, (3, -1, 5)),
                     (1,-1, 2, 3,-2), B, (-1,-2, 4))


        u,s,v = svdtrunc(fuselegs(fuselegs(C, 1, 2), 2, 2), maxdim=maxdim)
        normalize!(s)
        push!(matrices,
              mapcharges(x->div(x,2),
                         unfuseleg(u, 1, (E.legs[1], G.legs[1]))))
        E = unfuseleg(s*v, 2, (A.legs[3], B.legs[3]))
    end
    A = mps1.As[lx]
    B = mps2.As[lx]
    C = contract(contract(contract(E, (1,-1, 2), A, (-1, 3, 4)),
                          (1, 2, -1, 4), G, (3, -1, 5)),
                 (1,-1, 2, 3,-2), B, (-1,-2, 4))

    push!(matrices, mapcharges(x->div(x,2), fuselegs(C, 3, 2)))

    return MPState{Y}(lx, 2, dims, matrices, lx)
end

function _zipandgutzwiller_F23!(mps1::MPState{Y},
                                mps2::MPState{Y};
                                truncation::Bool=true,
                                maxdim::Int64=200) where {Y}
    #@assert mps1.d == mps2.d == 2
    lx = length(mps1)
    @assert length(mps2) == lx

    center_at!(mps1, 1)
    center_at!(mps2, 1)

    mps = MPState{Y}()

    ## NOTE: in order to make the gutzwiller projector respect the U1
    ## symmetry we need to do the follwoing. Assume the first mps
    ## corresponds to ↑ or 1 and second mps to ↓ or -1.
    T = eltype(Y)
    S = vtype(Y)

    Vd = VectorSpace{S}(0=>1, 1=>1)
    G = fill(one(T), (dual(Vd), mapcharges(x->2*x-1, Vd), Vd))

    Vdummy = VectorSpace{S}(0=>1)
    E = fill(one(T), (Vdummy, dual(Vdummy), Vdummy))

    u,s,v = _svd_(eye(T, Vdummy))
    for l=1:lx-1
        A = mps1.As[l]
        B = dual(mps2.As[l], conjugate=false)

        fswap = fermionswapgate(A.space[2], B.space[1])
        ##NOTE:
        # contracting (E3, A1)
        # First tensor has EA = E1 E2 A2 A3
        # contracting with swap
        # second tensor has EAX = E1 XA2 E3 A3
        # contracting (A2, G2)
        # Third tensor has EAG = E1 G1 G3 E3 A3
        # contracting (E3, B1) and (G3, B2)
        # Final tensor has E1 G1 A3 B3
        C = contract(contract(contract(contract(E, (1, 2, -1), A, (-1, 3, 4)),
                                       (1, -1, -2, 4), fswap, (-1, 2, 3, -2)),
                              (1, -1, 4, 5), G, (-1, 2, 3)),
                     (1,2, -2, -1, 4), B, (-1,-2, 3))

        if truncation
            u,s,v = svdtrunc(SymMatrix(C, [1,2], [3,4]), maxdim=maxdim)
        else
            u,s,v = _svd_(SymMatrix(C, [1,2], [3,4]))
        end
        normalize!(s)

        fnl = x->div(x+l-1, 2)
        fnd = x->div(x+1, 2)
        fnr = x->div(x+l, 2)
        push!(mps,
              mapcharges((fnl,fnd,fnr),
                         splitleg(u, 1, (E.space[1], G.space[2]))))
        E = splitleg(s*v, 2, (B.space[3], A.space[3]))
    end
    A = mps1.As[lx]
    B =  dual(mps2.As[lx], conjugate=false)
    fswap = fermionswapgate(A.space[2], B.space[1])
    C = contract(contract(contract(contract(E, (1, 2, -1), A, (-1, 3, 4)),
                                   (1, -1, -2, 4), fswap, (-1, 2, 3, -2)),
                          (1, -1, 4, 5), G, (-1, 2, 3)),
                 (1,2, -2, -1, 4), B, (-1,-2, 3))

    if truncation
        u,s,v = svdtrunc(SymMatrix(C, [1,2], [3,4]), maxdim=maxdim)
    else
        u,s,v = _svd_(SymMatrix(C, [1,2], [3,4]))
    end
    normalize!(s)

    C = splitleg(u*s*v, 1, (E.space[1], G.space[2]))
    fnl = x->div(x+lx-1, 2)
    fnd = x->div(x+1, 2)
    fnr = x->div(x+lx, 2)
    push!(mps, mapcharges((fnl,fnd,fnr), C))

    # could the below alone be the issue?!
    # is this correct (need explanation and stuff)
    mps.center = lx
    mps
end

function gutzwillerexact(mps1::MPState{Y},
                         mps2::MPState{Y};
                         mode::Symbol=:F23) where {Y}
    _zipandgutzwiller_F23(mps1, mps2, truncation=false)
end

function _tensorproductzip!(mps1::MPState{Y},
                            mps2::MPState{Y};
                            truncation::Bool=true,
                            maxdim::Int64=200,
                            verbose::Bool=true) where {Y}

    #@assert mps1.d == mps2.d == 2
    lx = length(mps1)
    @assert length(mps2) == lx

    center_at!(mps1, 1)
    center_at!(mps2, 1)

    mps = MPState{Y}()

    ## NOTE: in order to make the gutzwiller projector respect the U1
    ## symmetry we need to do the follwoing. Assume the first mps
    ## corresponds to ↑ or 1 and second mps to ↓ or -1.
    T = eltype(Y)
    S = vtype(Y)

    Vdummy = VectorSpace{S}(0=>1)
    E = fill(one(T), (Vdummy, dual(Vdummy), Vdummy))

    u,s,v = _svd_(eye(T, Vdummy))
    for l=1:lx-1
        A = mps1.As[l]
        B = dual(mps2.As[l], conjugate=false)

        fswap = fermionswapgate(A.space[2], B.space[1])
        ##NOTE:
        # contracting (E3, A1)
        # First tensor has EA = E1 E2 A2 A3
        # contracting with swap
        # second tensor has EAX = E1 XA2 E2 A3
        # contracting (EX2, B1)
        # Final tensor has E1 XA2 B2 B3 A3
        C = contract(contract(contract(E, (1, 2, -1), A, (-1, 3, 4)),
                              (1, -1, -2, 4), fswap, (-1, 2, 3, -2)),
                     (1, 2, -1, 5), B, (-1, 3, 4))

        if truncation
            u,s,v = svdtrunc(SymMatrix(C, [1,2,3], [4,5]), maxdim=maxdim)
        else
            u,s,v = _svd_(SymMatrix(C, [1,2,3], [4,5]))
        end
        normalize!(s)

        push!(mps, fuselegs(
            splitleg(u, 1, (E.space[1], A.space[2], B.space[2])), 2, 2))

        E = splitleg(s*v, 2, (B.space[3], A.space[3]))
    end
    A = mps1.As[lx]
    B = dual(mps2.As[lx], conjugate=false)
    fswap = fermionswapgate(A.space[2], B.space[1])
    C = contract(contract(contract(E, (1, 2, -1), A, (-1, 3, 4)),
                          (1, -1, -2, 4), fswap, (-1, 2, 3, -2)),
                 (1, 2, -1, 5), B, (-1, 3, 4))

    push!(mps, fuselegs(fuselegs(C, 2, 2), 3, 2))
    mps.center = lx
    mps
end

function _applygutzwillerzip!(mps::MPState{Y}) where{Y}
    T = eltype(mps)
    G = Tensor(ones, zero(U1),
               (U1Space([-1=>1, 1=>1]),
                dual(U1Space([-1=>1, 0=>2, 1=>1]))
                ))

    mpsgutz = MPState{Y}()

    center_at!(mps, 1)
    C = contract(G, (2, -1), mps.As[1], (1, -1, 3))
    for l=1:lx-1
        u,s,v = svdtrunc(fuselegs(C, 1, 2))

        fnl = x->div(x+lx-1, 2)
        fnd = x->div(x+1, 2)
        fnr = x->div(x+lx, 2)
        push!(mpsgutz, mapcharges((fnl,fnd,fnr),
                                  splitleg(u, 1, (C.space[1], G.space[1]))))
        E = s*v
        C = contract(E, (1, -1),
                     contract(G, (2, -1), mps.As[l+1], (1, -1, 3)),
                     (-1,2,3))
    end
    push!(mpsgutz, C)

    mpsgutz
end

function _zipandgutzwiller_F23_analysis!(mps1::MPState{Y},
                                         mps2::MPState{Y};
                                         maxdim::Int64=200) where {Y}

    lx = length(mps1)
    @assert length(mps2) == lx

    center_at!(mps1, 1)
    center_at!(mps2, 1)

    mps = MPState{Y}()
    ## NOTE: in order to make the gutzwiller projector respect the U1
    ## symmetry we need to do the follwoing. Assume the first mps
    ## corresponds to ↑ or 1 and second mps to ↓ or -1.
    T = eltype(Y)
    S = vtype(Y)

    Vd = VectorSpace{S}(0=>1, 1=>1)
    G = fill(one(T), (dual(Vd), mapcharges(x->2*x-1, Vd), Vd))

    Vdummy = VectorSpace{S}(0=>1)

    ### Making the F tensors
    Vend1 = rightspace(mps1)
    Vend2 = rightspace(mps2)
    Fs = Vector{SymTensor{S,T,3}}(undef, lx)
    F = fill(one(T), (Vend1, dual(Vend2), fuse(true, dual(Vend1), Vend2)))
    for l=lx:-1:2
        A = mps1.As[l]
        B = dual(mps2.As[l], conjugate=false)
        fswap = fermionswapgate(A.space[2], B.space[1])
        C = contract(contract(contract(contract(A, (1, 2, -1), F, (-1, 3, 4)),
                                       (1, 2, -1, 5), B, (3, 4, -1)),
                              (1, -2, -1, 4, 5), fswap, (2, 3, -1, -2)),
                     (1, 2, -1, -2, 4), G, (-1, 3, -2))
        u,s,v = _svd_(SymMatrix(C, [1,2], [3,4]))
        normalize!(s)
        F = splitleg(u*s, 1, (C.space[1], C.space[2]))
        Fs[l] = F
    end

    error_ortho = Vector{Float64}()
    error_ets = Vector{Float64}()
    E = fill(one(T), (Vdummy, Vdummy, dual(Vdummy)))
    for l=1:lx-1
        A = mps1.As[l]
        B = dual(mps2.As[l], conjugate=false)

        fswap = fermionswapgate(A.space[2], B.space[1])
        ##NOTE:
        # contracting (E3, A1)
        # First tensor has EA = E1 E2 A2 A3
        # contracting with swap
        # second tensor has EAX = E1 XA2 E3 A3
        # contracting (A2, G2)
        # Third tensor has EAG = E1 G1 G3 E3 A3
        # contracting (E3, B1) and (G3, B2)
        # Final tensor has E1 G1 A3 B3
        C = contract(contract(contract(contract(E, (1, 2, -1), A, (-1, 3, 4)),
                                       (1, -1, -2, 4), fswap, (-1, 2, 3, -2)),
                              (1, -1, 4, 5), G, (-1, 2, 3)),
                     (1,2, -2, -1, 4), B, (-1,-2, 3))


        # actual zipgutz process
        u,s,v = svdtrunc(SymMatrix(C, [1,2], [3,4]), maxdim=maxdim)
        normalize!(s)
        println(norm(SymMatrix(C, [1,2], [3,4])))
        # The ortho center matrix
        CFmat = SymMatrix(contract(C, (1,2,-1,-2), Fs[l+1], (-2,-1, 3)), [1, 2], [3])
        println(norm(SymMatrix(Fs[l+1], [2,1], [3])))
        norm_CFmat = norm(CFmat)
        println(norm_CFmat)
        println()
        # othonormal truncation error
        _u,_s,_v = svdtrunc(CFmat, maxdim=maxdim)
        normalize!(_s)
        push!(error_ortho, dot(_u*_s*_v, CFmat) / norm_CFmat)

        # zipgutz truncation error
        CtrFmat = (u*s*v) * SymMatrix(Fs[l+1], [2,1], [3])
        norm_CtrFmat = norm(CtrFmat)
        push!(error_ets, dot(CtrFmat, CFmat) / (norm_CFmat*norm_CtrFmat))

        fnl = x->div(x+l-1, 2)
        fnd = x->div(x+1, 2)
        fnr = x->div(x+l, 2)
        push!(mps,
              mapcharges((fnl,fnd,fnr),
                         splitleg(u, 1, (E.space[1], G.space[2]))))
        E = splitleg(s*v, 2, (B.space[3], A.space[3]))
    end
    A = mps1.As[lx]
    B =  dual(mps2.As[lx], conjugate=false)
    fswap = fermionswapgate(A.space[2], B.space[1])
    C = contract(contract(contract(contract(E, (1, 2, -1), A, (-1, 3, 4)),
                                   (1, -1, -2, 4), fswap, (-1, 2, 3, -2)),
                          (1, -1, 4, 5), G, (-1, 2, 3)),
                 (1,2,-2,-1,4), B, (-1,-2, 3))

    u,s,v = svdtrunc(SymMatrix(C, [1,2], [3,4]), maxdim=maxdim)
    normalize!(s)

    push!(error_ortho, 1.0)
    push!(error_ets, 1.0)

    C = splitleg(u*s*v, 1, (E.space[1], G.space[2]))
    fnl = x->div(x+lx-1, 2)
    fnd = x->div(x+1, 2)
    fnr = x->div(x+lx, 2)
    push!(mps, mapcharges((fnl,fnd,fnr), C))

    mps.center = lx
    mps, error_ortho, error_ets
end
