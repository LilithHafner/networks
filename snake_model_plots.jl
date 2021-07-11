include("analysis.jl")
include("snake_model.jl")
using BenchmarkTools
using Plots
using Random: shuffle

cache = Dict()
function sm(args...; trial=1)
    key = (args, trial)
    if haskey(cache, key)
        return cache[key]
    end
    result_tuple = @timed snake_model(args...)
    result = Dict(zip(keys(result_tuple), result_tuple))
    result[:graph] = result[:value]
    @time result[:wedges],result[:lcc] = wedge_distribution(result[:graph])
    #result[:lcc_nonzero] = (sum(wd)-wd[1])/sum(wd)
    #result[:lcc_mean] = sum((c-1)*wd[c] for c ∈ eachindex(wd))/sum(wd)
    cache[key] = result
end

x,y = [],[]
for trial in 1:1
    for speed in shuffle(0:.025:1)
                  # k, node, edge, width, sp, dupli, trial
        result = sm(4, 1000, 1000, 20, speed, false, trial=trial)
        push!(x, speed)
        push!(y, result)
    end
end
s1 = scatter()
scatter!(s1, x, [i[:lcc] for i in y], label="lcc")
#scatter!(s1, x, [i[:lcc_mean] for i in y], label="llc mean")
s2 = scatter(x, [i[:time] for i in y]*1e3, label="time (ms)")
#g = grasshop_populate([100], [1,1,1,1], 10^3/binomial(100+4-1,4))
