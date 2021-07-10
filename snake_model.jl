function snake(collection::AbstractVector{T}, size::Integer,
               width::Integer, speed::Real) where T <: Any
    out = Vector{T}(undef, size)
    for i in 1:width
        out[i] = rand(collection)
    end
    for i in width+1:size
        out[i] = rand(rand() < speed ? collection : out[i-width:i-1])
    end
    out
end

function snake_model(hyperdegree::Integer, nodes::Integer, edges::Integer,
                     width::Integer, speed::Real)
    source = snake(1:nodes, hyperdegree*edges, width, speed)
    [source[hyperdegree*(e-1)+1:hyperdegree*e] for e in 1:edges]
end

#=
display(transpose(reshape(snake(10:99, 40, 5, .5), 5, 8)))
display(snake_model(5, 99, 8, 5, .5))


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
