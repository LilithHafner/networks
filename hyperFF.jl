# The hyperFF paper:
# Evolution of Real-world Hypergraphs: Patterns and Models without Oracles
# by Yunbum Kook, Jihoon Ko, Kijung Shin
# https://arxiv.org/pdf/2008.12729.pdf
#
# How I found hyperFF https://arxiv.org/abs/2104.14332

using DataStructures: Queue, enqueue!, dequeue!
include("ordered_count_vector.jl")

function hyperFF(n::Int, p::Real, q::Real)
    hyperFF(n, (v) -> length(v) == 0 ? true : rand() < (length(v) == 1 ? p : q))
end
function hyperFF(n::Int, Θ::Function)
    edges = []
    adjacency = [OrderedCountVector(Int) for _ ∈ 1:n]
    for i ∈ 2:n#i is the new node
        #Could set P(Θ([_])) ≠ 1 for a not fully connected graph
        !Θ([i]) && continue

        # Could select initial ambassador differently.
        for ambassador ∈ burn((n, _) -> Θ([i, n]), adjacency, rand(1:i-1))#Select nodes to create edges to
            push!(adjacency[i], ambassador)#Only include narrow edge
            push!(adjacency[ambassador], i)#Could include expanded as well
            push!(edges, collect(burn((n, b) -> Θ(vcat(n, collect(b))), adjacency, ambassador, i)))
        end
    end
    edges
end

function burn(Θ, adjacency, ambassador)
    q = Queue{typeof(ambassador)}()
    enqueue!(q, ambassador)
    burn!(Set(ambassador), Θ, adjacency, q)
end
function burn(Θ, adjacency, ambassador, origin)
    q = Queue{typeof(ambassador)}()
    enqueue!(q, ambassador)
    burn!(Set([origin, ambassador]), Θ, adjacency, q)
end
function burn!(burned, Θ, adjacency, q)
    while !isempty(q)
        ambassador = dequeue!(q)
        for neighbor ∈ keys(adjacency[ambassador])
            !Θ(neighbor, burned) && break
            if neighbor ∉ burned
                push!(burned, neighbor)
                enqueue!(q, neighbor)
            end
        end
    end
    burned
end

display(hyperFF(10, .5, .5))
