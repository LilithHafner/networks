#TODO include anything worth keeping from statistics.jl

##Fast O(size)
function clustering_coefficients(edge_list::Vector{Vector{T}}) where T <: Any
    #Θ(size)
    #Computes a nearly uniform and not independent sample.
    #Non-uniformity is at most 1 node.
    #Approximates the average local and the global clustering coefficients

    size = 0
    n2e = Dict{T,Vector{Set{T}}}() #Θ(size)
    for edge ∈ edge_list
        size += length(edge)
        for n ∈ edge
            if n ∈ keys(n2e)
                push!(n2e[n], Set(edge))
            else
                push!(n2e, n=>[Set(edge)])
            end
        end
    end


    #Each triangle contains 3 closed wedges
    fuel = size + 100
    wedges, closed_wedges = 0,0#For global clustering coefficient
    lcc_nan = []#For local clustering coefficient
    nodes = collect(keys(n2e))
    for n ∈ shuffle(nodes)# Θ(fuel) = θ(size)
        fuel <= 0 && break
        local_closed_wedges, local_wedges = 0,0
        edges = shuffle(collect(setdiff(e, n) for e in n2e[n]))
        for (i,e1) ∈ enumerate(edges)
            for e2 ∈ edges[i+1:length(edges)]
                local_wedges += 1
                fuel -= length(e2)
                for n1 ∈ e1
                    for e11 ∈ n2e[n1]
                        if !isdisjoint(e11, e2)#TODO change documentation for isdisjoint in Julia.
                            local_closed_wedges += 1
                            @goto break_label
                        end
                    end
                end
                @label break_label
            end
        end
        wedges += local_wedges
        closed_wedges += local_closed_wedges
        push!(lcc_nan, local_closed_wedges/local_wedges)
    end
    lcc = filter!((!) ∘ isnan, lcc_nan)
    Dict{Symbol,Any}(
        :gcc => closed_wedges/wedges,
        :lcc => lcc,
        :lcc_nan => lcc_nan,
        :acc => sum(lcc)/length(lcc),
        :wedges => wedges*length(nodes)/length(lcc_nan),
        :triangles => closed_wedges*length(nodes)/length(lcc_nan)/3,
        :sample_nodes => length(lcc_nan),
        :sample_proportion => length(lcc_nan)/length(nodes),
    )
end

##Slow
function wedge_distribution(edges::AbstractVector{E}) where {N <: Any, E <: AbstractVector{N}}
    #O(k^2*Σdegree^2) ~ O(k^2dm)
    #This could but does not use concepts of relationship to self.
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
