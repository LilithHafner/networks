##Helpers
max_node(edge_list) = maximum(maximum.(edge_list))
function adjacency_matrix(edge_list)::BitMatrix
    nodes = max_node(edge_list)
    #Convert hyper-edges to cliques
    out = falses(nodes, nodes)
    for hyper_edge in edge_list
        for a in hyper_egde
            for b in hyper_egde
                out[a,b] = true
            end
        end
    end
    out
end
function adjacency_table(edge_list)::Vector{Set{Int}}
    #Convert hyper-edges to cliques
    out = [Set() for node in 1:max_node(edge_list)]
    for hyper_edge in edge_list
        for a in hyper_edge
            for b in hyper_edge
                push!(out[a], b)
            end
        end
    end
    out
end

##Statistics
clique_count(max_k, graph::Vector{Set{T}}) where T <: Integer =
    clique_count(max_k, graph, Set(1:length(graph)))
function clique_count(max_k, graph::Vector{Set{T}}, mask) where T <: Integer
    #WARNING Exponential runtime in k. Perhaps θ(n^k*k^2)?
    #[(number of k-cliques in graph ∩ mask) for k in 1:max_k]
    if max_k == 0 return [] end
    if max_k == 1 return [length(mask)] end
    if isempty(mask) return fill(0,max_k) end
    pushfirst!(sum(clique_count(max_k-1, graph, graph[node] ∩ mask)
        for node in mask).÷max_k, length(mask))
end

##Analysis
#TODO Automate discrimination searching analysis
#TODO interpret that in an application minded way
#     i.e. literally using a more sociology & GWSS definition of discrimination

using Plots
unzip(array) = [getindex.(array, i) for i in 1:length(array[1])]
function plot_clique_counts(models::Vector{T}, n; kw...) where T <: Function
    graphss = [[model() for i in 1:n] for model in models]
    graphss, plot_clique_counts(graphss; kw...)
end
function plot_clique_counts(graphss; max_k=4, labels=:auto, kw...)
    if labels == :auto
        labels = 1:length(graphss)
    end
    out = scatter()
    for (graphs, label) in zip(graphss, labels)
        k,x = unzip(reduce(vcat, collect.(
            enumerate.(clique_count.(max_k, adjacency_table.(graphs))))))
        scatter!(k,x; label=label, kw...)
    end
    out
end

##Usage
include("grasshop_populate.jl")
include("coinflip_populate.jl")
include("hyper_stochastic_block.jl")
function compare_clique_count(populate!s, n, args...; kw...)
    models = [(() -> hyper_stochastic_block_edge_count(populate!, args...))
              for populate! in populate!s]
    plot_clique_counts(models, n; kw...)[2]
end

function f(gen_k=3, group_sizes=[20,10,5], edges=100, show_k=4)
    compare_clique_count(
        (grasshop_populate!, coinflip_populate!), 10,
        gen_k, group_sizes, edges;
        max_k = show_k, labels = ("grasshop", "coinflip"),
        markeralpha=.1)
end

f()

#TODO test all statistics functions (esp. helpers)
#TODO use a better output format
