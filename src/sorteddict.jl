struct SortedDict{K,V} <: AbstractDict{K,V}
    keys :: Vector{K}
    values :: Vector{V}
end

SortedDict{K,V}() where{K,V} = SortedDict(Vector{K}(), Vector{V}())
function SortedDict{K,V}(pairs::Vector{Pair{K,V}}) where {K,V}
    if !issorted(ps, by=first)
        pairs = sort(ps, by=first)
    end
    SortedDict{K,V}(map(first, ps), map(last, ps))
end

SortedDict(ps::Pair...)= SortedDict(ps)
SortedDict{K,V}(ps::Pair{K,V}...) where {K,V} = SortedDict{K,V}(ps)
function SortedDict{K,V}(kv) where {K,V}
    d = SortedDict{K,V}()
    sizehint!(d, length(kv))
    for (k,v) in kv
        push!(d, k=>v)
    end
    return d
end

@inline Base.keytype(d::SortedDict{K,V}) where{K,V} = K
#@inline Base.valuetype(d::SortedDict{K,V}) where{K,V} = V

Base.length(d::SortedDict) = length(d.keys)
function Base.sizehint!(d::SortedDict, sz)
    sizehint!(d.keys, sz)
    sizehint!(d.values, sz)
    d
end

@inline function findkey(d::SortedDict, k)
    key = convert(keytype(d), k)
    if !isequal(key, k)
        throw(ArgumentError("$(limitrepr(k)) is not a valid key for type $K"))
    end
    i = searchsortedfirst(d.keys, key)
    @inbounds i<=length(d.keys) && d.keys[i] == k && return i, true
    i, false
end

@inline function haskey(d::SortedDict, k)
    i, exactfound = findkey(d, k)
    exactfound
end

@inline function get(d::SortedDict, k, _default)
    i, exactfound = findkey(d, k)
    if !exactfound
        _default
    end
    d.values[i]
end

@inline function Base.getindex(d::SortedDict, k)
    i, exactfound = findkey(d, k)
    if !exactfound
        throw(KeyError(k))
    end
    d.values[i]
end


function Base.setindex!(d::SortedDict, v, k)
    i, exactfound = findkey(d, k)
    if exactfound
        d.values[i] = v
    else
        insert!(d.keys, i, k)
        insert!(d.values, i, v)
    end
    d
end

function Base.iterate(d::SortedDict, i = 1)
    @inbounds if i > length(d)
        return nothing
    else
        return (d.keys[i] => d.values[i]), i+1
    end
end
