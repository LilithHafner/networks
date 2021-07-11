#TODO include anything worth keeping from statistics.jl

function wedge_distribution(edges::AbstractVector{E}) where {N <: Any, E <: AbstractVector{N}}
    #This could but does not use concepts of relationship to self.
    #This is very slow O(k^2*Σdegree^2)
    edges_by_node = Dict{N,Set{E}}()
    nodes_by_node = Dict{N,Set{N}}()
    for edge ∈ edges
        for node ∈ edge
            if haskey(edges_by_node, node)
                push!(edges_by_node[node], edge)
                union!(nodes_by_node[node], edge)
            else
                push!(edges_by_node, node=>Set([edge]))
                push!(nodes_by_node, node=>Set(edge))
            end
        end
    end

    #out = fill(0, length(edges)-1) This was wrong and would be too slow.
    open,closed = 0,0
    for edge1 ∈ edges
        for apex ∈ edge1
            for edge2 ∈ edges_by_node[apex]
                if edge1 ≠ edge2
                    #closures = 0
                    for node1 ∈ edge1
                        if node1 ≠ apex
                            for node2 ∈ edge2 ∩ nodes_by_node[node1]
                                if node2 ≠ apex
                                    closed += 1
                                    @goto escape_label
                                end
                            end
                        end
                    end
                    open += 1
                    @label escape_label
                    #out[closures+1] += 1
                end
            end
        end
    end
    #resize!(out,findlast(0 .≠ out))
    open+closed, closed/(open+closed)#, nodes_by_node, edges_by_node
end

##Example usage (see also snake_model_plots.jl)
#=
include("snake_model.jl")
g = snake_model(3, 30, 60, 5, .5)
wd = wedge_distribution(g)
lcc_nonzero = (sum(wd)-wd[1])/sum(wd)
lcc_mean = sum((c-1)*wd[c] for c ∈ eachindex(wd))/sum(wd)
=#
