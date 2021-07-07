include("grasshop_populate.jl")
include("hyper_stochastic_block.jl")
grasshop(args...) = hyper_stochastic_block(grasshop_populate!, args...)
grasshop_edge_count(args...) =
    hyper_stochastic_block_edge_count(grasshop_populate!, args...)

using BenchmarkTools
using Plots
using Random

function meter_stick(k, edges)
    e = []
    for i in 1:edges
        push!(e,rand(1:100,k))
    end
    e
end
function plot_meter_stick()
    plots = []
    ks = shuffle(1:2:50)
    @time kt = [minimum([@elapsed meter_stick(k, 10^4)
            for _ in 1:10]) for k in ks] ./ 10^4
    push!(plots, scatter(ks,kt.*10^9,
        ylabel="ns per edge", xlabel="dimension", title="Meter Stick",
        legend=false, ylims=(0,Inf), xlims=(0,Inf), ywiden=true, xwiden=true))

    es = shuffle(800:800:30000)
    @time et = [minimum([@elapsed meter_stick(5, e)
            for _ in 1:10]) for e in es]
    push!(plots, scatter(es,et*10^3,
        ylabel="total ms", xlabel="edges",
        legend=false, ylims=(0,Inf), xlims=(0,Inf), ywiden=true, xwiden=true))

    plot(plots...,layout=(length(plots),1))
end

time_in_meters(args...; kwargs...) = time_in_meters(grasshop, args...; kwargs...)
time_in_meters_edge_count(args...; kwargs...) = time_in_meters(grasshop_edge_count, args...; kwargs...)
function time_in_meters(f::Function, args...; n=:auto, runtime=.1)
    k = args[1]
    if n == :auto
        time_guess = @elapsed f(args...)
        n = Int(ceil(runtime/time_guess))
    end

    edge_counts = []

    times = [@elapsed push!(edge_counts, length(f(args...)))
        for i in 1:n]
    _, time, edge_count = minimum(zip(times./edge_counts, times, edge_counts))
    meter_stick_time = minimum(@elapsed meter_stick(k, edge_count)
        for i in 1:n)
    (meters=time/meter_stick_time, seconds=time, edges=edge_count, samples=n)
end

function myscatter(xs, fs...; labels=nothing, title="title", kwargs...)
    plot = scatter()
    if :legend ∉ keys(kwargs)
        kwargs = merge(kwargs, Dict(:legend => labels != nothing))
    end
    labels = labels != nothing ? labels : fill("auto", length(fs))
    time = @elapsed begin
        xs = shuffle(xs)
        for (f,l) in zip(fs,labels)
            ys = [f(xi) for xi in xs]
            scatter!(xs,ys; label=l,
                ylims=(0,Inf), xlims=(0,Inf), ywiden=true, xwiden=true,
                kwargs...)
        end
    end
    title!(title*string(Int(floor(time*1000)))*"ms")
    plot
end

function benchmark(quick=false)
    samples = quick ? 1 : 3
    disperse = quick ? 4 : 1
    plots = [
        myscatter(1000:disperse*10^5:10^6, e->1000*
                minimum(@elapsed grasshop_edge_count(3,fill(1000,10), e)
                    for i in 1:samples), e->1000*
                minimum(@elapsed meter_stick(3, e)
                    for i in 1:samples),
            xlabel="Edges", ylabel="ms",
            title="Hyperdegree 3\n10 groups of 1000 nodes\n10^5 edges\n",
            labels=["Actual", "Theoretical Minimum"],
            legend=(.12,.85)),
        myscatter(1000:disperse*5*10^3:10^5, n->1000*
                minimum(@elapsed grasshop_edge_count(3,fill(n÷10,10), 10^5)
                    for i in 1:samples), n->1000*
                minimum(@elapsed meter_stick(3, 10^5)
                    for i in 1:samples),
            xlabel="Total nodes", ylabel="ms",
            title="Hyperdegree 3\n10 groups\n10^5 edges\n"),
        myscatter(1:disperse*3:50, g->1000*
                minimum(@elapsed grasshop_edge_count(3,fill(10^4÷g,g), 3*10^5)
                    for i in 1:samples), g->1000*
                minimum(@elapsed meter_stick(3, 3*10^5)
                    for i in 1:samples),
            xlabel="# of groups", ylabel="ms",
            title="Hyperdegree 3\n10^4 total nodes\n3*10^5 edges\n"),
        myscatter([2:disperse:7...], k->1000*
                minimum(@elapsed grasshop_edge_count(k,fill(100,10), 10^5)
                    for i in 1:samples), k->1000*
                minimum(@elapsed meter_stick(k, 10^5)
                    for i in 1:samples),
            xlabel="Hyperdegree", ylabel="ms",
            title="10 groups of 100 nodes\n10^5 edges\n"),
    ]
    plot(plots..., size=(700,700))
    #plots
end
##Old
function populate_edges_to_probability(group_sizes, m, edges)
    size = prod(group_sizes[m])/
           prod(factorial(count(x->x==i, m)) for i in 1:length(group_sizes))
    probability = edges/size
    group_sizes, m, probability
end
function populate_time_per_edge_probability(args... ; n=100)
    lengths = []
    times = [@elapsed push!(lengths,length(populate(args...)))
            for i in 1:n]
    t,l = minimum(zip(times./lengths,lengths))
    l, t*1e9
end
function populate_time_per_edge(args... ; n=100)
    populate_time_per_edge_probability(populate_edges_to_probability(args...)..., n=n)
end
