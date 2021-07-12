include("analysis.jl")
include("snake_model.jl")
using BenchmarkTools
using Plots
using Random: shuffle

cache = Dict()
function generate(f, args...; trial=1)
    key = (f, args, trial)
    if haskey(cache, key)
        return cache[key]
    end
    result_tuple = @timed f(args...)
    result = Dict(zip(keys(result_tuple), result_tuple))
    result[:graph] = result[:value]
    merge!(result, clustering_coefficients(result[:graph]))
    cache[key] = result
end


    x,y = [],[]
    for trial in 1:2
        for width in shuffle(1:128)
            speed = .2
            sz = 40
            result = generate(snakemodel2,
                sz*[9,16,50,120,20],
                sz*[5,10,30,30,30,20,20,10,4,1],
                width, speed, trial=trial)
            push!(x, width)
            push!(y, result)
        end
    end
    xy = zip(x,y)
    xs = sort(collect(Set(x)))

#scatter(xlim=(0,1), ylim=(0,1))
#scatter(ylim=(0,100))
scatter()

scatter!(x, [y[:acc] for y in y], label="acc", ma=.1)
scatter!(x, [y[:gcc] for y in y], label="gcc", ma=.1)
plot!(xs, [mean(y[:acc] for (x,y) in xy if x == xi) for xi in xs], label="acc mean")
plot!(xs, [mean(y[:gcc] for (x,y) in xy if x == xi) for xi in xs], label="gcc mean")

#end

#scatter!(s1, x, [i[:lcc_mean] for i in y], label="llc mean")
#s2 = scatter(x, [i[:time] for i in y]*1e3, label="time (ms)")
#g = grasshop_populate([100], [1,1,1,1], 10^3/binomial(100+4-1,4))

## Time
#=

#for width in [1,2,4,8,16,32,64,128]
#println(width)
width = 20
speed = .2
x,y = [],[]
#runs = reduce(vcat, [[(trial,sz) for trial in 1:500÷sz] for sz in 1:10:500])
#for (trial, sz) in shuffle(runs)
for trial in 1:5
    for sz in shuffle(1:4:50)#0:.025:1)
                                     # k, node, edge, width, sp, dupli, trial
        #result = generate(snakemodel1, 4, 1000, 1000, 20, speed, false, trial=trial)
        result = generate(snakemodel2,
            sz*[9,16,50,120,20],
            sz*[5,10,30,30,30,20,20,10,4,1],
            width, speed, trial=trial)
        push!(x, sum(length.(result[:graph])))
        push!(y, result[:time]*1000)
    end
end
xy = zip(x,y)
sx = sort(collect(Set(x)))

scatter(xlim=(0,1), ylim=(0,1))
scatter(ylim=(0,100))
scatter()
plot!(sx,
    [mean(y for (x,y) in xy if x == xi) for xi in sx],
    label="mean")
my = [minimum(y for (x,y) in xy if x == xi) for xi in sx]
plot!(sx, my, label="minimum", legend=false)
plot!([0,maximum(sx)], [0, maximum(sx)*sum(my.*sx)/sum(sx.*sx)], label="minimum", legend=false)

plot!(sx,
    [median(y for (x,y) in xy if x == xi) for xi in sx],
    label="median")
scatter!(x, y, label="time (ms)", smooth=true)

#end

s1
#scatter!(s1, x, [i[:lcc_mean] for i in y], label="llc mean")
#s2 = scatter(x, [i[:time] for i in y]*1e3, label="time (ms)")
#g = grasshop_populate([100], [1,1,1,1], 10^3/binomial(100+4-1,4))
=#
