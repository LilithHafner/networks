struct OrderedCountVector{K <: Any}
    keycounts :: Vector{Tuple{K,Int}}
    key2ind :: Dict{K,Int}
end
OrderedCountVector(K) = OrderedCountVector(Tuple{K,Int}[], Dict{K,Int}())

import Base.push!
function push!(v::OrderedCountVector{K}, k::K) where K
    if k in keys(v.key2ind)
        i = v.key2ind[k]
        count = v.keycounts[i][1] + 1
        while i > 1 && counts[i-1] < count
            v.key2ind[i] = v.key2ind[keys[i-1]]
            v.keycounts[i] = v.keycounts[i-1]
            i -= 1
        end
        v.key2ind[k] = i
        v.keycounts[i] = (k, count)
    else
        push!(v.keycounts, (k,1))
        push!(v.key2ind, (k => length(v.keycounts)))
    end
end

import Base.keys
keys(v::OrderedCountVector{K}) where K = getindex.(v.keycounts, 1)
