using Plots

function unrank2d(n)
    # Roughly O(1) for n :: Int64.
    # Constant: 2ns = sqrt runtime +5% -0%
    # Max input: 2^(NN-2)-1 for IntNN
    a = Int(round(sqrt(2n)))
    b = n-(a*(a-1))>>1
    #Depending on benchmark,
    # >>1 makes unrank2d 1.1-5x faster than ÷2.
    a,b
end

### unrank2d Tests
#=
println(collect(unrank2d.(1:10)))
@assert collect(unrank2d.(1:10)) == [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2),
                                     (3, 3), (4, 1), (4, 2), (4, 3), (4, 4)]
@assert unrank2d(87109723841) == (417396, 222131)

#@benchmark unrank2d(9720290439024)

function test(n)
    try
        @assert typeof(n) <: Integer
        @assert n >= 1
        a,b = unrank2d(n)
        @assert a >= 1
        @assert b >= 1

        @assert b <= a
        @assert a*(a-1)÷2+b == n
    catch
        println(typeof(n))
        println(n)
    end
end

EVALS = 1000
SAMPLES = 30
function plt(x)
    x = shuffle(collect(x))
    t = [minimum([@elapsed (for i in 1:EVALS unrank2d(xi) end)
        for j in 1:SAMPLES])/EVALS for xi in x]
    for xi in x test(xi) end
    scatter(x,t)
end

plot(plt(1:10^2),
     plt(2^62-1:-(2^62-1)÷10^2:1),
    layout=(2,1))
=#
