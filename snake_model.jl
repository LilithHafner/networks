function snake(collection::AbstractVector{T}, size::Integer,
               width::Integer, speed::Real) where T <: Any
    #O(size) constant factor = 12.5 ns for
    #    collection size 1–typemax(Int), size 100–10^9, width 1–typemax(Int), speed -Inf–Inf

    out = Vector{T}(undef, size)
    k0,k1 = length(collection)/speed,width/(1-speed)#*

    for i in 1:min(width, size)
        out[i] = rand(collection)
    end
    for i in width+1:size
        r = rand()#*
        out[i] = r < speed ? collection[T(floor(r*k0)+1)] : out[i-1-T(floor((r-speed)*k1))]#*
    #    out[i] = rand() < speed ? rand(collection) : out[rand(i-width:i-1)] #This simplification costs a 13–20% increase in runtime
    end
    out
end

function snake_model(hyperdegree::Integer, nodes::Integer, edges::Integer,
                     width::Integer, speed::Real, duplicates::Bool)
    source = snake(1:nodes, hyperdegree*edges, width, speed)
    if duplicates
        [source[hyperdegree*(e-1)+1:hyperdegree*e] for e in 1:edges]
    else
        collect.(collect(Set(Set(source[hyperdegree*(e-1)+1:hyperdegree*e])) for e in 1:edges))
    end
end


#=
#Basic inspection
display(transpose(reshape(snake(10:99, 40, 5, .5), 5, 8)))
display(snake_model(5, 99, 8, 5, .5))
display(snake_model(3, 3, 60, 5, .5))

#Benchmarking
using BenchmarkTools
mean((@belapsed snake(1:Int(10)^18, $(Int(floor(10^i))), 200, .5) seconds=1)/Int(floor(10^i))*1e9 for i in 1:5)

#Displaying graphs
using LightGraphs
using GraphPlot
function graph_plot(node_count, edges)
    graph = SimpleGraph(node_count)
    for (a,b) in edges
        add_edge!(graph, a,b)
    end
    gplot(graph)
end
graph_plot(60, snake_model(2, 60, 60, 5, .5))
=#
