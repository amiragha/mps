function generatempo(model::UnitCellQModel)
    typeof(model.qtype) == SpinType || error("Can generate MPO only for SpinType models.")
    model.lattice.bc == :OBC || error("Can only generate MPO for :OBC boundary condition!")
    D = dimension(model)
    d = model.qtype.d
    n_sites = prod(model.lattice.sizes)
    T = eltype(model.inters[1])

    allterms = Vector{QInteraction}()
    for is in Iterators.product([1:l for l in model.lattice.sizes]...)
        for interaction in model.inters
            ns = interaction.ucidxs
            offs = interaction.offsets
            indexes = [sitelinearindex(model.lattice, ns[i], offs[i] .+ is)
                       for i in 1:support(interaction)]
            if all(indexes .> 0)
                perm = sortperm(indexes)
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
    MatrixProductOperator(n_sites, d, dims, tensors)
end

function generatesymmpo(model::UnitCellQModel)
    typeof(model.qtype) == SpinType || error("Can generate MPO only for SpinType models.")
    model.lattice.bc == :OBC || error("Can only generate MPO for :OBC boundary condition!")
    typeof(model.inters[1]) <: SymQModelInteraction || error("Only symmetric (U1) is allowed!")

    D = dimension(model)
    d = model.qtype.d
    d == 2 || error("only works for spin 1/2 for now!")
    n_sites = prod(model.lattice.sizes)
    T = eltype(model.inters[1])

    allterms = Vector{QTerm{SymMatrix{T}}}()
    for is in Iterators.product([1:l for l in model.lattice.sizes]...)
        for interaction in model.inters
            ns = interaction.ucidxs
            offs = interaction.offsets
            indexes = [sitelinearindex(model.lattice, ns[i], offs[i] .+ is)
                       for i in 1:support(interaction)]
            if all(indexes .> 0)
                perm = sortperm(indexes)
                for ops in interaction.terms
                    permutedops = ops[perm]
                    termops = (permutedops[1:end-1]..., interaction.amp * permutedops[end])
                    push!(allterms, QTerm(Tuple(indexes[perm]), termops))
                end
            end
        end
    end
    sort!(allterms, by=x->x.sites[1])

    dims = ones(Int, n_sites+1)
    pointer=1
    tensors = Vector{SymTensor{Float64, 4}}(undef, n_sites)
    chrs = [0]
    lleg = STLeg(+1, chrs, [1])
    o1leg = STLeg(+1, [0,1], [1,1])
    o2leg = STLeg(-1, [0,1], [1,1])
    symI = eye(T, collect(0:d-1), ones(Int, d))
    pchrs = Dict{Int, Vector{QTerm}}()
    for n in 1:n_sites
        # find the new charge for each previous-term that still
        # continues to the next site
        nextchrs = Int[0]
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
        chrs = Int[]
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
            idx = searchsortedfirst(chrs, 0)
            chrdims[idx] += 1
            rleg = STLeg(-1, chrs, chrdims)
        else
            chrs == [0] || error("$chrs")
            chrdims == [1] || error("$chrdims")
            rleg = STLeg(-1, [0], [1])
        end

        dims[n+1] = fulldims(rleg)


        W = fill(zero(T), 0, (lleg, o1leg, rleg, o2leg))

        index0000 = index_sector(W, (0,0,0,0))
        index0101 = index_sector(W, (0,1,0,1))
        l00, r00 = size(W.nzblks[index0000])[[1,3]]
        # better logic is needed for this!
        if n > 1
            W.nzblks[index0000][1,1,1,1] = one(T)
            W.nzblks[index0101][1,1,1,1] = one(T)
        end
        if n < n_sites
            W.nzblks[index0000][l00,1,r00,1] = one(T)
            W.nzblks[index0101][l00,1,r00,1] = one(T)
        end

        for lterm in lterms
            op = lterm.ops[1]
            op.charge == 0 || error()
            W.nzblks[index0000][l00, 1, 1, 1] += get_sector(op, (0,0))
            W.nzblks[index0101][l00, 1, 1, 1] += get_sector(op, (1,1))
        end

        rows = Dict(c=>1 for c in lleg.chrs)
        rows[0] = 2
        cols = Dict(c=>1 for c in chrs)
        cols[0] = 2

        for (pchr, pterms) in pchrs
            for pterm in pterms
                if pterm.sites[1] == n
                    op = pterm.ops[1]
                    if support(pterm) == 1
                        # pchr goes to zero
                        row ,col = rows[pchr], cols[0]
                        for i in 1:length(op.sects)
                            c1, c2 = op.sects[i]
                            idx = index_sector(W, (pchr, c1, 0, c2))
                            size(op.nzblks[i]) == (1, 1) || error()
                            W.nzblks[idx][row, 1, 1, 1] = op.nzblks[i][1,1]
                        end
                        rows[pchr] += 1
                    else
                        #pchr goes to pchr+op.charge
                        row ,col = rows[pchr], cols[pchr+op.charge]
                        for i in 1:length(op.sects)
                            c1, c2 = op.sects[i]
                            idx = index_sector(W, (pchr, c1, pchr+op.charge, c2))
                            size(op.nzblks[i]) == (1, 1) || error()
                            W.nzblks[idx][row, 1, col, 1] = op.nzblks[i][1,1]
                        end
                        rows[pchr] += 1
                        cols[pchr+op.charge] += 1
                        push!(nextterms[pchr+op.charge], removehead(pterm))
                    end
                else
                    # pchr goes to pchr
                    row, col = rows[pchr], cols[pchr]
                    for i in 1:length(symI.sects)
                        c1, c2 = symI.sects[i]
                        idx = index_sector(W, (pchr, c1, pchr, c2))
                        size(symI.nzblks[i]) == (1, 1) || error()
                        W.nzblks[idx][row, 1, col, 1] = symI.nzblks[i][1,1]
                    end
                    rows[pchr] += 1
                    cols[pchr] += 1
                    push!(nextterms[pchr], pterm)
                end
            end
        end

        for sterm in sterms
            op = sterm.ops[1]
            row, col = rows[0], cols[op.charge]
            for i in 1:length(op.sects)
                c1, c2 = op.sects[i]
                idx = index_sector(W, (0, c1, op.charge, c2))
                size(op.nzblks[i]) == (1, 1) || error()
                W.nzblks[idx][getdim(lleg, 0), 1, col, 1] = op.nzblks[i][1,1]
            end
            cols[op.charge] += 1
            push!(nextterms[op.charge], removehead(sterm))
        end

tensors[n] =  W
pchrs = nextterms
lleg = STLeg(+1, rleg.chrs, rleg.dims)
end
SymMatrixProductOperator(n_sites, d, dims, tensors)
end
