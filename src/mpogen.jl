function generatempo(model::UnitCellQModel)
    typeof(model.qtype) == SpinType || error("Can generate MPO only for SpinType models.")
    model.lattice.bc == :OBC || error("Can only generate MPO for :OBC boundary condition!")
    D = dimension(model)
    d = model.qtype.d
    n_sites = prod(model.lattice.sizes)
    T = eltype(model.inters[1])

    allterms = Vector{QInteraction}()
    for is in Iterators.product([1:l for l in model.lattice.sizes]...)
        #println(is)
        for interaction in model.inters
            ns = interaction.ucidxs
            offs = interaction.offsets
            indexes = [sitelinearindex(model.lattice, ns[i], offs[i] .+ is)
                       for i in 1:support(interaction)]
            #println(indexes)
            if all(indexes .> 0)
                perm = sortperm(indexes)
                terms = [term[perm] for term in interaction.terms]
                push!(allterms,
                      QInteraction(interaction.amp, Tuple(indexes[perm]), terms))
            end
        end
    end
    sort!(allterms, by=x->x.sites[1])

    #println(allterms)
    # for term in allterms
    #     println(term)
    # end

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

        dims[n] = rdim
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
                        W[row, :, col, :] = pterm.amp * pterm.terms[i][1]
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
    eltype(model.inters) == SymQModelInteraction || error("Only symmetric (U1) is allowed!")

    D = dimension(model)
    d = model.qtype.d
    n_sites = prod(model.lattice.sizes)
    T = eltype(model.inters[1])

    allterms = Vector{SymQInteraction}()
    for is in Iterators.product([1:l for l in model.lattice.sizes]...)
        #println(is)
        for interaction in model.inters
            ns = interaction.ucidxs
            offs = interaction.offsets
            indexes = [sitelinearindex(model.lattice, ns[i], offs[i] .+ is)
                       for i in 1:support(interaction)]
            #println(indexes)
            if all(indexes .> 0)
                perm = sortperm(indexes)
                terms = [term[perm] for term in interaction.terms]
                push!(allterms,
                      SymQInteraction(interaction.amp, Tuple(indexes[perm]), terms))
            end
        end
    end
    sort!(allterms, by=x->x.sites[1])

    #println(allterms)
    # for term in allterms
    #     println(term)
    # end

    # to continue from here!
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

        dims[n] = rdim
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
                        W[row, :, col, :] = pterm.amp * pterm.terms[i][1]
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
