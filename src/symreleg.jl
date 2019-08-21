function uniquesorted(vec::Vector{T}) where{T}
    unqvec = T[]
    refs = Int[]
    push!(unqvec, vec[1])
    idx = 1
    push!(refs, idx)
    for i=2:length(vec)
        if vec[i] != unqvec[end]
            push!(unqvec, vec[i])
            idx += 1
            push!(refs, idx)
        else
            push!(refs, idx)
        end
    end
    unqvec, refs
end

function inv_perm(perm::Vector{Int})
    inv = zeros(Int, length(perm))
    for i in eachindex(perm)
        inv[perm[i]] = i
    end
    inv
end

unzip(a) = map(x->getfield.(a,x), fieldnames(eltype(a)))

# fuse n consequative legs `l` and `l+n-1` to a new leg with sign `sign`
# function fuse_conseqlegs(sten::SymTensor{T, N},
#                          sign::Int, l::Int, n::Int=2; debug=false) where{T<:Number, N}
#     @assert abs(sign) == 1
#     #@assert 0 < l < N+n-2
#     @assert 0 < l < N+2-n

#     fusedsectors = NTuple{N-n+1, Int}[]
#     patranges = UnitRange{Int}[]

#     csects = sten.sects
#     fchrdict = Dict{Int, FusedCharge{n}}()
#     signs = Tuple(sten.legs[c].sign for c in l:l+n-1)
#     for i in eachindex(csects)
#         sect = [csects[i]...]

#         pat = Tuple(sect[l:l+n-1])
#         spat = sign .* signs .* pat
#         fchr = sign * sum(signs .* [sect[c] for c in l:l+n-1])
#         patrange = UnitRange{Int}(1,1)
#         if haskey(fchrdict, fchr)
#             if haskey(fchrdict[fchr].pats, spat)
#                 patrange = fchrdict[fchr].pats[spat]
#             else
#                 fdim = prod([getdim(sten.legs[c+l-1], pat[c]) for c in 1:n])
#                 dim = fchrdict[fchr].dim
#                 fchrdict[fchr].dim += fdim
#                 patrange = dim+1:dim+fdim
#                 fchrdict[fchr].pats[spat] = patrange
#             end
#         else
#             fdim = prod([getdim(sten.legs[c+l-1], pat[c]) for c in 1:n])
#             patrange = 1:fdim
#             fchrdict[fchr] = FusedCharge(fchr, fdim, Dict(spat => patrange))
#         end

#         fsect = vcat( sect[1:l-1], fchr, sect[l+n:end])
#         push!(fusedsectors, Tuple(fsect))
#         push!(patranges, patrange)
#     end

#     fusedsectperm = _sectors_sortperm(fusedsectors)
#     new_sectors, refs = uniquesorted(fusedsectors[fusedsectperm])
#     refs = refs[inv_perm(fusedsectperm)]

#     # the infos are: new_sectors, refs, patranges

#     # make the new leg
#     fleg = STLeg(sign, unzip([(k, fchrdict[k].dim) for k in sort(collect(keys(fchrdict)))])...)

#     new_legs = Tuple(vcat([sten.legs[1:l-1]...], fleg, [sten.legs[l+n:end]...]))

#     new_nzblks = Array{T, N-n+1}[]
#     for sect in new_sectors
#         push!(new_nzblks,Array{T, N-n+1}(undef, [getdim(new_legs[c], sect[c]) for c in 1:N-n+1]...))
#     end
#     #indexing = [1:end for i=1:N-n+1]...
#     for i in eachindex(sten.nzblks)
#         s = size(sten.nzblks[i])
#         new_nzblks[refs[i]][[1:s[i] for i=1:l-1]..., patranges[i],[1:s[i] for i=l+n:N]...] =
#             reshape(sten.nzblks[i], s[1:l-1]..., prod(s[l:l+n-1]), s[l+n:end]...)
#     end

#     SymTensor(sten.charge, new_legs, new_sectors, new_nzblks)
# end

# fuse n consequative legs `l` and `l+n-1` to a new leg with sign `sign`
function fuselegs(sten::SymTensor{Tv, N},
                  sign::Int,
                  l::Int,
                  n::Int=2;
                  debug=false) where{Tv<:Number, N}
    @assert abs(sign) == 1
    @assert 0 < l < N+2-n

    fsects = NTuple{N-n+1, Int}[]
    patranges = UnitRange{Int}[]
    signs = Tuple(sten.legs[i].sign for i in l:l+n-1)
    fleg, fcdict = fuse(sign, sten.legs[l:l+n-1])

    csects = sten.sects
    for i in eachindex(csects)
        csect = csects[i]
        pat = Tuple(csect[l:l+n-1])
        spat = sign .* signs .* pat
        fc = sum(spat)
        patrange = fcdict[fc].pats[spat]
        fsect = (csect[1:l-1]..., fc, csect[l+n:end]...)
        push!(fsects, fsect)
        push!(patranges, patrange)
    end

    fsectperm = _sectors_sortperm(fsects)
    newsects, refs = uniquesorted(fsects[fsectperm])
    refs = refs[inv_perm(fsectperm)]

    newlegs = Tuple(vcat([sten.legs[1:l-1]...], fleg, [sten.legs[l+n:end]...]))

    newnzblks = Array{Tv, N-n+1}[]
    for sect in newsects
        push!(newnzblks,zeros(Tv, [getdim(newlegs[c], sect[c])
                                  for c in 1:N-n+1]...))
    end

    for i in eachindex(sten.nzblks)
        s = size(sten.nzblks[i])
        newnzblks[refs[i]][[1:s[i] for i=1:l-1]..., patranges[i],[1:s[i] for i=l+n:N]...] =
            reshape(sten.nzblks[i], s[1:l-1]..., prod(s[l:l+n-1]), s[l+n:end]...)
    end

    SymTensor(sten.charge, newlegs, newsects, newnzblks)
end

function delinsert(tuple::NTuple{N, T}, items::NTuple{M, T}, index::Int) where{T, N, M}
    (tuple[1:index]..., items..., tuple[index+1:N])
end

# defuses a leg into some given legs
function defuse_leg(sten::SymTensor{T, N},
                    l::Int,
                    legs::NTuple{M, STLeg}) where {T<:Number, N, M}
    @assert 0 < l <= N
    #println(sten.legs)
    #println(legs)
    if length(legs) == 1
        if legs[1].sign == sten.legs[l].sign
            sten.legs[l] == legs[1] || error("Not the same legs", sten.legs[l], legs[1])
            return sten
        else
            stenleg = change_sign(sten.legs[l])
            stenleg == legs[1] || (error("Not the same legs", stenleg, legs[1]))
            return change_legsign(sten, l)
        end
    end

    sects = NTuple{N+M-1, Int}[]
    nzblks = Array{T, N+M-1}[]

    sign = sten.legs[l].sign
    csects = sten.sects
    fchrdict = Dict{Int, Tuple{Vector{NTuple{M,Int}}, Vector{NTuple{M,Int}}}}()
    for i in eachindex(csects)
        charge = sign * csects[i][l]
        if haskey(fchrdict, charge)
            pats, patdims = fchrdict[charge]
        else
            pats, patdims = _possible_fuse_patterns(charge, legs)
            #println(charge, pats, patdims)
            ##TODO: is sorting required or can be avoided with a better design?!
            patperms = _sectors_sortperm(pats)
            pats, patdims = pats[patperms], patdims[patperms]
            #println(pats, patdims)
            fchrdict[charge] = (pats, patdims)
        end
        old_nzblock =  sten.nzblks[i]
        s = size(old_nzblock)
        pivot = 0
        for patidx in eachindex(pats)
            sl = prod(patdims[patidx])
            push!(sects, (csects[i][1:l-1]..., pats[patidx]..., csects[i][l+1:N]...))
            push!(nzblks, reshape(old_nzblock[[1:s[n] for n=1:l-1]...,pivot+1:pivot+sl,[1:s[n] for n=l+1:N]...],
                                  s[1:l-1]...,patdims[patidx]...,s[l+1:N]...))
            pivot += sl
        end
    end


    # now check to see if the new legs can fuse into the original leg
    ###TODO: there should be an options to remove this step for
    ###efficientcy when it is not needed!
    for charge in keys(fchrdict)
        cidx = findall(sten.legs[l].chrs .== sign * charge)[1]
        sten.legs[l].dims[cidx] == sum([prod(patdims) for patdims in fchrdict[charge][2]]) ||
            error("For defuse, dimensions don't match for charge " , sign*charge)
        #, " ", sten.legs[l].dims[cidx], " == ", fchrdict[charge]

    end

    #println(fchrdict)
    new_legs = (sten.legs[1:l-1]...,legs...,sten.legs[l+1:end]...)
    #println(sten.sects)
    #println(sects)

    ##NOTE: at the end of defuse the sectors may no longer be sorted!
    sectperm = _sectors_sortperm(sects)
    SymTensor(sten.charge, new_legs, sects[sectperm], nzblks[sectperm])
end
