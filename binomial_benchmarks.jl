using Plots
using BenchmarkTools
using Random

function measure2(n,k,samples=30, evals=3)
    g() = binomial(n,k)
    measure(g,samples,evals)
end
function measure(f, samples, evals)
    minimum([@elapsed (for i in 1:evals f() end)
        for j in 1:samples])/evals
end
bexp(x) = BigInt(floor(exp(log(10)*BigFloat(x))))

function to_color(x)
    x = [isnan(i) ? 1 : i for i in x]
    x = max.(min.(x,1e4),1e-4)
    x ./= exp(mean(log.(x)))
    x = 1 ./ (1 .+ 1 ./ x)
    [RGB{Float32}(i, .2, 1-i) for i in x]
end
function my_plot(n, k, t, f, g, title, xlabel, alpha)
    scatter(x, t, markercolor=to_color(g.(n,k)),
        markersize=2, markeralpha=alpha, markerstrokewidth = 0,
        smooth=true, linecolor=:black, linestyle=:dash,
        legend=false, yformatter=:scientific,
        title=title, ylabel="runtime (s)", xlabel=xlabel)
end

i1 = 0:10^-2.5:1
n1 = shuffle([bexp(x*100) for x in i1])
k1 = shuffle([BigInt(floor(x*1000)) for x in i1])
t1 = [measure2(i,j) for (i,j) in zip(n1,k1)]
large = my_plot(n1, k1, t1,
    (n,k) -> log(max(binomial(n,k),1))/log(2)*k, (n,k)-> log(n)/k,
    "For larger cases\n n<1e100, k<1e3\nO(output_size*k)",
    "output_size*k", .4)

i2 = 0:10^-3:1
n2 = shuffle([Int(floor(exp(x*3*log(10)))) for x in i2])
k2 = shuffle([Int(floor(x*6+1)) for x in i2])
filtered = [(n,k) for (n,k) in zip(n2,k2) if n>2*k]
n2,k2 = [getindex.(filtered,i) for i in 1:2]
t2 = [measure2(i,j,300,30) for (i,j) in zip(n2,k2)]
small = my_plot(n2, k2, t2,
    (n,k) -> k, (n,k)-> log(n)/k,
    "For small cases\nbinomial(n,k) < 2^64\nO(k)",
    "k", .1)

plot(small, large, layout=(1,2))

#x1 = log.(max.(binomial.(n1,k1),1))./log(2).*k1
#print(sum(t1.*x1)/sum(x.*x))
#print(sum(t2.*k2)/sum(k2.*k2))
