function generatempo(model::UnitCellQModel)
    typeof(model.qtype) == SpinType || error("Can generate MPO only for SpinType models.")
    model.lattice.bc == :OBC || error("Can only generate MPO for :OBC boundary condition!")
    D = dimension(model)
    d = model.qtype.d

    # We need a map between sites and their linear ordering
    #cartesian coordinates of (uc.n, ly, lx)

    # traverse through the lattice generate all terms and sort them
    allterms = Vector{QInteraction}()
    for is in Iterators.product([1:l for l in model.lattice.sizes]...)
            for interaction in model.inters
                ns = interaction.sites
                offs = interaction.offsets
                indexes = [sitelinearindex(model.lattice, ns,offs .+ is)
                           for i in 1:N]
                perm = sortperm(indexes)
                terms = [permute(term, perm) for term in interaction.terms]
                allterms.push!(QInteraction(interaction.amp, indexes[perm], terms))
            end
    end
    sort!(allterm, by=x->x.sites[1])

    pointer=1
    tensors = Vector{Array{Float64, 4}}(undef, n_sites)
    ldim = 1
    fterms = Vector{QInteraction}()
    cterms = Vector{QInteraction}()
    pterms = Vector{QInteraction}()
    for n in 1:n_sites
        # traves allterms to see how many sterms we have
        lterms = Vector{QInteraction}()
        sterms = Vector{QInteraction}()
        while pointer <= length(allterms) && allterms[pointer].sites[1] == n
            if support(allterms[i]) > 1
                push!(sterms, allterms[pointer])
            else
                push!(lterms, allterms[pointer])
            end
            pointer += 1
        end
        rdim = 2 + n_sterms + n_cterms

        # MPO matrix generation
        W = zeros(T, ldim, d, rdim, d)
        W[1,:,1,:] = I(d)
        W[ldim, :, rdim, :] = I(d)

        #W[ldim, :, 1, :] = sum(lterms).terms[1]
        for t in 1:n_pterms
            if pterms[t].sites[1] == n
                if support(pterms[t], 1)
                    W[t+1, :, 1, :] = pterm[t].amp * pterms[t].terms
                else
                    W[t+1, :, c+1, :] = pterms[t].terms[i][1]
                    push!(nexterms, removehead(pterms[t]))
                    c += 1
                end
            else
                W[t+1, :, c+1, :] = I(d)
                push!(nexterms, pterms[t])
                c += 1
            end
        end

        for t in 1:n_lterms
            W[ldim, :, 1, :] += lterms[t].terms[i][1]
        end

        for t in 1:n_sterms
            W[ldim, :, c+1, :] = lterms[t].terms[i][1]
            push!(nexterms, removehead(pterms[t]))
            c += 1
        end

        tensors.push!(W)
    end
    ldim = rdim
end
# return MPO
