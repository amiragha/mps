function generatempo(model::UnitCellQModel;
                     symmetry::Symbol=:AUTO,
                     verbose::Bool=false)
    if model.symmetry == Trivial
        return _generatempo_nosym(model, verbose=verbose)

    elseif model.symmetry == U1
        return _generatempo_sym(model, verbose=verbose)

    elseif eltype(model.inters) <: FermionQModelInteraction
        error("JW not yet implemented! Do it manually and input the spin model for the mpo!")
    else
        error("Dont' recognize $(eltype(model.inters)) for MPO generation!")
    end
end

function _generatempo_nosym(model::UnitCellQModel;
                            verbose::Bool=false)
    typeof(model.qtype) == SpinType ||
        error("Can generate MPO only for Bosonic models.")
    D = dimension(model)
    model.lattice.bcs[D] == :OBC ||
        @warn "MPO generated periodic in X direction! discouraged!"

    d = model.qtype.d
    n_sites = prod(model.lattice.sizes)
    T = eltype(model.inters[1])

    allterms = Vector{QInteraction}()
    for is in Iterators.product([1:l for l in model.lattice.sizes]...)
        for interaction in model.inters
            ns = interaction.ucidxs
            offs = interaction.offsets
            indexes = zeros(Int, support(interaction))
            isinside = true
            for i in 1:support(interaction)
                index, crossings = sitelinearindex(model.lattice, ns[i], offs[i] .+ is)
                if isnothing(index)
                    isinside = false
                    break
                end
                indexes[i] = index
            end
            if isinside
                perm = sortperm(indexes)
                verbose && println("Adding $(interaction.amp) between $(indexes[perm])")
                terms = [term[perm] for term in interaction.terms]
                push!(allterms,
                      QInteraction(interaction.amp, Tuple(indexes[perm]), terms))
            end
        end
    end
    sort!(allterms, by=x->x.sites[1])

    dims = ones(Int, n_sites+1)
    pointer=1
    tensors = Vector{Array{Float64, 4}}(undef, n_sites)
    ldim = 1
    n_pcols = 0
    pterms = Vector{QInteraction}()
    nextterms = Vector{QInteraction}()
    for n in 1:n_sites
        lterms = Vector{QInteraction}()
        sterms = Vector{QInteraction}()
        n_ncols = 0
        while pointer <= length(allterms) && allterms[pointer].sites[1] == n
            if support(allterms[pointer]) > 1
                push!(sterms, allterms[pointer])
                n_ncols += length(allterms[pointer].terms)
            else
                push!(lterms, allterms[pointer])
            end
            pointer += 1
        end

        n_pcols = 0
        for pterm in pterms
            if pterm.sites[1] > n || support(pterm) > 1
                n_pcols += length(pterm.terms)
            end
        end

        if n < n_sites
            rdim = 2 + n_pcols + n_ncols
        else
            n_ncols == 0 || error()
            n_pcols == 0 || error()
            rdim = 1
        end

        dims[n+1] = rdim
        # MPO matrix generation
        W = zeros(T, ldim, d, rdim, d)

        # better logic is needed for this!
        if n > 1
            W[1,:,1,:] = I(d)
        end
        if n < n_sites
            W[ldim, :, rdim, :] = I(d)
        end

        #W[ldim, :, 1, :] = sum(lterms).terms[1]
        row = 2
        col = 2
        for pterm in pterms
            if pterm.sites[1] == n
                if support(pterm) == 1
                    for i in eachindex(pterm.terms)
                        W[row, :, 1, :] = pterm.amp * pterm.terms[i][1]
                        row+=1
                    end
                else
                    for i in eachindex(pterm.terms)
                        W[row, :, col, :] = pterm.terms[i][1]
                        row+=1
                        col+=1
                    end
                    push!(nextterms, removehead(pterm))
                end
            else
                for i in eachindex(pterm.terms)
                    W[row, :, col, :] = I(d)
                    row+=1
                    col+=1
                end
                push!(nextterms, pterm)
            end
        end
        n_pcols = col - 2

        for lterm in lterms
            for i in eachindex(lterm.terms)
                W[ldim, :, 1, :] += lterm.amp * lterm.terms[i][1]
            end
        end

        for sterm in sterms
            for i in eachindex(sterm.terms)
                W[ldim, :, col, :] = sterm.terms[i][1]
                col += 1
            end
            push!(nextterms, removehead(sterm))
        end

        tensors[n] =  W
        pterms = nextterms
        nextterms = Vector{QInteraction}()
        ldim = rdim
    end
    MPOperator(n_sites, d, dims, tensors)
end

function _generatempo_sym(model::UnitCellQModel;
                          verbose::Bool=false)
    typeof(model.qtype) == SpinType ||
        error("Can generate MPO only for Bosonic models.")
    D = dimension(model)
    model.lattice.bcs[D] == :OBC ||
        @warn "MPO generated periodic in X direction! discouraged!"

    # typeof(model.inters[1]) <: SymQModelInteraction ||
    #     error("Only symmetric (U1) is allowed!")

    d = model.qtype.d
    d == 2 || error("only works for spin 1/2 for now!")
    n_sites = prod(model.lattice.sizes)
    T = eltype(model.inters[1])

    sz,sp,sm = spinoperators(1/2, symmetry=U1)
    dict = Dict{Symbol, SymMatrix{U1, T}}()
    dict[:Z] = sz
    dict[:P] = sp
    dict[:M] = sm
    dict[:I] = eye(T, U1Space(0=>1, 1=>1))

    allterms = Vector{QTerm{SymMatrix{U1,T}}}()
    for is in Iterators.product([1:l for l in model.lattice.sizes]...)
        for interaction in model.inters
            ns = interaction.ucidxs
            offs = interaction.offsets
            indexes = zeros(Int, support(interaction))
            isinside = true
            for i in 1:support(interaction)
                index, crossings = sitelinearindex(model.lattice, ns[i], offs[i] .+ is)
                if isnothing(index)
                    isinside=false
                    break
                end
                indexes[i] = index
            end
            if isinside
                perm = sortperm(indexes)
                for qaterm in interaction.terms
                    ops = [dict[qaterm.ops[i]] for i in eachindex(qaterm.ops)]
                    permutedops = ops[perm]
                    termops = (permutedops[1:end-1]..., qaterm.amp * interaction.amp * permutedops[end])
                    push!(allterms, QTerm(Tuple(indexes[perm]), termops))
                end
            end
        end
    end
    sort!(allterms, by=x->x.sites[1])

    pointer=1
    mpo = MPOperator{SymTensor{U1, T, 4}}()
    chrs = [zero(U1)]
    Vw = U1Space(0=>1)
    Vo = U1Space(0=>1, 1=>1)
    symI = eye(T, Vo)
    pchrs = Dict{U1Charge, Vector{QTerm}}()
    for n in 1:n_sites
        # find the new charge for each previous-term that still
        # continues to the next site
        nextchrs = U1[zero(U1)]
        for (pchr, pterms) in pchrs
            for pterm in pterms
                if pterm.sites[1] > n
                    push!(nextchrs, pchr)
                elseif support(pterm) > 1
                    push!(nextchrs, pchr + pterm.ops[1].charge)
                end
            end
        end

        # making local-terms and starting-terms and for starting terms
        # that continue find the charge the go to
        lterms = Vector{QTerm}()
        sterms = Vector{QTerm}()
        while pointer <= length(allterms) && allterms[pointer].sites[1] == n
            if support(allterms[pointer]) > 1
                push!(sterms, allterms[pointer])
                push!(nextchrs, allterms[pointer].ops[1].charge)
            else
                push!(lterms, allterms[pointer])
            end
            pointer += 1
        end

        sort!(nextchrs)
        chrs = U1[]
        chrdims = Int[]
        push!(chrs, nextchrs[1])
        push!(chrdims, 1)
        for c in nextchrs[2:end]
            if c == chrs[end]
                chrdims[end] += 1
            else
                push!(chrs, c)
                push!(chrdims, 1)
            end
        end

        nextterms = Dict(chr => Vector{QTerm}() for chr in chrs)
        # make the new leg
        if n < n_sites
            idx = searchsortedfirst(chrs, U1(0))
            chrdims[idx] += 1
            Vr = dual(U1Space([chrs[i]=>chrdims[i] for i in eachindex(chrs)]))
        else
            chrs == [U1(0)] || error("$chrs")
            chrdims == [1] || error("$chrdims")
            Vr = dual(U1Space(zero(U1)=>1))
        end

        W = fill(zero(T), U1(0), (Vw, Vo, Vr, dual(Vo)))

        l00, r00 = size(W.blocks[Sector{U1}(0,0,0,0)])[[1,3]]
        # better logic is needed for this!
        if n > 1
            W.blocks[Sector{U1}(0,0,0,0)][1,1,1,1] = one(T)
            W.blocks[Sector{U1}(0,1,0,1)][1,1,1,1] = one(T)
        end
        if n < n_sites
            W.blocks[Sector{U1}(0,0,0,0)][l00,1,r00,1] = one(T)
            W.blocks[Sector{U1}(0,1,0,1)][l00,1,r00,1] = one(T)
        end

        for lterm in lterms
            op = lterm.ops[1]
            op.charge == 0 || error()
            W.blocks[Sector{U1}(0,0,0,0)][l00, 1, 1, 1] += op[Sector{U1}(0,0)]
            W.blocks[Sector{U1}(0,1,0,1)][l00, 1, 1, 1] += op[Sector{U1}(1,1)]
        end

        rows = Dict(c=>1 for c in charges(Vw))
        rows[0] = 2
        cols = Dict(U1(c)=>1 for c in chrs)
        cols[0] = 2

        for (pchr, pterms) in pchrs
            for pterm in pterms
                if pterm.sites[1] == n
                    op = pterm.ops[1]
                    if support(pterm) == 1
                        # pchr goes to zero
                        row ,col = rows[pchr], cols[U1(0)]
                        for (s,blk) in op.blocks
                            c1, c2 = s.charges
                            sect = Sector{U1}(pchr, c1, 0, c2)
                            size(blk) == (1, 1) || error()
                            W.blocks[sect][row, 1, 1, 1] = blk[1,1]
                        end
                        rows[pchr] += 1
                    else
                        #pchr goes to pchr+op.charge
                        row ,col = rows[pchr], cols[pchr+op.charge]
                        for (s,blk) in op.blocks
                            c1, c2 = s.charges
                            sect = Sector{U1}(pchr, c1, pchr+op.charge, c2)
                            size(blk) == (1, 1) || error()
                            W.blocks[sect][row, 1, col, 1] = blk[1,1]
                        end
                        rows[pchr] += 1
                        cols[pchr+op.charge] += 1
                        push!(nextterms[pchr+op.charge], removehead(pterm))
                    end
                else
                    # pchr goes to pchr
                    row, col = rows[pchr], cols[pchr]
                    for (s,blk) in symI.blocks
                        c1, c2 = s.charges
                        sect = Sector{U1}(pchr, c1, pchr, c2)
                        size(blk) == (1, 1) || error()
                        W.blocks[sect][row, 1, col, 1] = blk[1,1]
                    end
                    rows[pchr] += 1
                    cols[pchr] += 1
                    push!(nextterms[pchr], pterm)
                end
            end
        end

        for sterm in sterms
            op = sterm.ops[1]
            row = rows[U1(0)]
            col = cols[op.charge]
            for (s,blk) in op.blocks
                c1, c2 = s.charges
                sect = Sector{U1}(0, c1, op.charge, c2)
                size(blk) == (1, 1) || error()
                W.blocks[sect][dim(Vw, U1(0)), 1, col, 1] = blk[1,1]
            end
            cols[op.charge] += 1
            push!(nextterms[op.charge], removehead(sterm))
        end

push!(mpo, W)
pchrs = nextterms
Vw = dual(Vr)
end
mpo
end

function _generatempo_infinite(model::UnitCellQModel;
                               verbose::Bool=false)
    typeof(model.qtype) == SpinType ||
        error("Can generate MPO only for Bosonic models!")
    D = dimension(model)
    lx = model.lattice.sizes[D]
    T = eltype(model.inters[1])
    model.lattice.bcs[D] in [:INF] ||
        error("boundary conitions inf not :INF!")

    typeof(model.inters[1]) <: SymQModelInteraction ||
        error("Only symmetric (U1) is allowed!")

    n_sites = prod(model.lattice.sizes[1:D-1])
    xrange = largestxrange(model)
    println(xrange)
    mpo = generatesymmpo(changeboundary(changesize(model, D, 2*(xrange+1)),
                                        D, :OBC))

    b = n_sites * xrange
    for i = 1:n_sites
        mpo.tensors[b+i] == mpo.tensors[b+i+n_sites] ||
            error("mpo not converged!")
    end

    tensors = Vector{SymTensor{T, 4}}()

    dims = Int[]
    for l = 1:lx
        for i=1:n_sites
            push!(dims, size(mpo.tensors[b+i], 1))
            push!(tensors, mpo.tensors[b+i])
        end
    end
    push!(dims, size(tensors[end], 3))
    return MPOperator(n_sites * lx, mpo.d, dims, tensors)
end
# function generateinfinitempo(model::UnitCellQModel)
#     typeof(model.qtype) == SpinType ||
#         error("Can generate MPO only for SpinType models.")
#     model.lattice.bc in [:PBCYAPBCX, :PBCX, :PBCYPBCX] ||
#         error("unrecognized boundary condition $boundary!")

#     typeof(model.inters[1]) <: SymQModelInteraction ||
#         error("Only symmetric (U1) is allowed!")

#
#     d = model.qtype.d
#     d == 2 || error("only works for spin 1/2 for now!")
#
#     T = eltype(model.inters[1])

#     ## notes: each term defined in the unit cell should be ordered by
#     ## indexes on the -inf to +inf index system. Then is should be
#     ## translated along the x-axis (last axis) so that the current
#     ## x-positions can be at any position from the beginnning to the
#     ## end, and all of those terms should be included in the infinte MPO
#     lattice = changeboundary(model.lattice, D, :INF)
#     allterms = Vector{QTerm{SymMatrix{T}}}()
#     for is in Iterators.product([1:l for l in model.lattice.sizes[1:D-1]]...)
#         for interaction in model.inters
#             ns = interaction.ucidxs
#             offs = interaction.offsets
#             indexes = zeros(Int, support(interaction))
#             isinside = true
#             for i in 1:support(interaction)
#                 index, crossings = sitelinearindex(lattice, ns[i], offs[i] .+ is)
#                 if isnothing(index)
#                     isinside = false
#                     break
#                 end
#                 indexes[i] = index
#             end
#             if isinside
#             perm = sortperm(indexes)
#             permindexes = indexes[perm]
#             if support(interaction) == 1
#                 for ops in interaction.terms
#                     n = mod1(indexes[1], n_sites)
#                     push!(lterms[n], interaction.amp * ops[1])
#                 end
#                 break
#             end
#             for ops in interaction.terms
#                 permutedops = ops[perm]
#                 termops = (permutedops[1:end-1]..., interaction.amp * permutedops[end])
#                 pointer = 1
#                 n = mod1(permindexes[pointer], n_sites)
#                 push!(sterms[n], termops[pointer])
#                 index = permindexes[pointer] + 1
#                 pointer += 1
#                 while pointer < support(interation)
#                     if index == permindexes[pointer]
#                         pointer += 1
#                     else
#                         push!(cterms[n], I)
#                     end
#                     index += 1
#                 end
#                 n = mod1(permindexes[pointer], n_sites)
#                 push!(eterms[n], termops[pointer])

#                         if n == permindexes[pointer]
#                             push!(cterms[n], termops[pointer])
#                             pointer +=1
#                         else
#                             # Add I
#                         end
#                     end
#                 end
#             end
#         end
#         end
#     end
#     sort!(allterms, by=x->x.sites[1])

#     dims = ones(Int, n_sites+1)
#     pointer=1
#     tensors = Vector{SymTensor{Float64, 4}}(undef, n_sites)
#     chrs = [0]
#     lleg = STLeg(+1, chrs, [1])
#     o1leg = STLeg(+1, [0,1], [1,1])
#     o2leg = STLeg(-1, [0,1], [1,1])
#     symI = eye(T, collect(0:d-1), ones(Int, d))
#     pchrs = Dict{Int, Vector{QTerm}}()
#     for n in 1:n_sites
#         # find the new charge for each previous-term that still
#         # continues to the next site
#         nextchrs = Int[0]
#         for (pchr, pterms) in pchrs
#             for pterm in pterms
#                 if pterm.sites[1] > n
#                     push!(nextchrs, pchr)
#                 elseif support(pterm) > 1
#                     push!(nextchrs, pchr + pterm.ops[1].charge)
#                 end
#             end
#         end

#         # making local-terms and starting-terms and for starting terms
#         # that continue find the charge the go to
#         lterms = Vector{QTerm}()
#         sterms = Vector{QTerm}()
#         while pointer <= length(allterms) && allterms[pointer].sites[1] == n
#             if support(allterms[pointer]) > 1
#                 push!(sterms, allterms[pointer])
#                 push!(nextchrs, allterms[pointer].ops[1].charge)
#             else
#                 push!(lterms, allterms[pointer])
#             end
#             pointer += 1
#         end

#         sort!(nextchrs)
#         chrs = Int[]
#         chrdims = Int[]
#         push!(chrs, nextchrs[1])
#         push!(chrdims, 1)
#         for c in nextchrs[2:end]
#             if c == chrs[end]
#                 chrdims[end] += 1
#             else
#                 push!(chrs, c)
#                 push!(chrdims, 1)
#             end
#         end

#         nextterms = Dict(chr => Vector{QTerm}() for chr in chrs)
#         # make the new leg
#         if n < n_sites
#             idx = searchsortedfirst(chrs, 0)
#             chrdims[idx] += 1
#             rleg = STLeg(-1, chrs, chrdims)
#         else
#             chrs == [0] || error("$chrs")
#             chrdims == [1] || error("$chrdims")
#             rleg = STLeg(-1, [0], [1])
#         end

#         dims[n+1] = fulldims(rleg)


#         W = fill(zero(T), 0, (lleg, o1leg, rleg, o2leg))

#         index0000 = index_sector(W, (0,0,0,0))
#         index0101 = index_sector(W, (0,1,0,1))
#         l00, r00 = size(W.blocks[index0000])[[1,3]]
#         # better logic is needed for this!
#         if n > 1
#             W.blocks[index0000][1,1,1,1] = one(T)
#             W.blocks[index0101][1,1,1,1] = one(T)
#         end
#         if n < n_sites
#             W.blocks[index0000][l00,1,r00,1] = one(T)
#             W.blocks[index0101][l00,1,r00,1] = one(T)
#         end

#         for lterm in lterms
#             op = lterm.ops[1]
#             op.charge == 0 || error()
#             W.blocks[index0000][l00, 1, 1, 1] += get_sector(op, (0,0))
#             W.blocks[index0101][l00, 1, 1, 1] += get_sector(op, (1,1))
#         end

#         rows = Dict(c=>1 for c in lleg.chrs)
#         rows[0] = 2
#         cols = Dict(c=>1 for c in chrs)
#         cols[0] = 2

#         for (pchr, pterms) in pchrs
#             for pterm in pterms
#                 if pterm.sites[1] == n
#                     op = pterm.ops[1]
#                     if support(pterm) == 1
#                         # pchr goes to zero
#                         row ,col = rows[pchr], cols[0]
#                         for i in 1:length(op.sects)
#                             c1, c2 = op.sects[i]
#                             idx = index_sector(W, (pchr, c1, 0, c2))
#                             size(op.blocks[i]) == (1, 1) || error()
#                             W.blocks[idx][row, 1, 1, 1] = op.blocks[i][1,1]
#                         end
#                         rows[pchr] += 1
#                     else
#                         #pchr goes to pchr+op.charge
#                         row ,col = rows[pchr], cols[pchr+op.charge]
#                         for i in 1:length(op.sects)
#                             c1, c2 = op.sects[i]
#                             idx = index_sector(W, (pchr, c1, pchr+op.charge, c2))
#                             size(op.blocks[i]) == (1, 1) || error()
#                             W.blocks[idx][row, 1, col, 1] = op.blocks[i][1,1]
#                         end
#                         rows[pchr] += 1
#                         cols[pchr+op.charge] += 1
#                         push!(nextterms[pchr+op.charge], removehead(pterm))
#                     end
#                 else
#                     # pchr goes to pchr
#                     row, col = rows[pchr], cols[pchr]
#                     for i in 1:length(symI.sects)
#                         c1, c2 = symI.sects[i]
#                         idx = index_sector(W, (pchr, c1, pchr, c2))
#                         size(symI.blocks[i]) == (1, 1) || error()
#                         W.blocks[idx][row, 1, col, 1] = symI.blocks[i][1,1]
#                     end
#                     rows[pchr] += 1
#                     cols[pchr] += 1
#                     push!(nextterms[pchr], pterm)
#                 end
#             end
#         end

#         for sterm in sterms
#             op = sterm.ops[1]
#             row, col = rows[0], cols[op.charge]
#             for i in 1:length(op.sects)
#                 c1, c2 = op.sects[i]
#                 idx = index_sector(W, (0, c1, op.charge, c2))
#                 size(op.blocks[i]) == (1, 1) || error()
#                 W.blocks[idx][getdim(lleg, 0), 1, col, 1] = op.blocks[i][1,1]
#             end
#             cols[op.charge] += 1
#             push!(nextterms[op.charge], removehead(sterm))
#         end

# tensors[n] =  W
# pchrs = nextterms
# lleg = STLeg(+1, rleg.chrs, rleg.dims)
# end
# MPOperator(n_sites, d, dims, tensors)
# end
