using StaticArrays
using BenchmarkTools

unrank1d(n) = n # 0.05ns ≈ .05ns * d^2
function unrank2d(n)
    # O(1) for n :: Int64.
    # Constant: 1.8ns = sqrt runtime +5% -0% (VERY FAST) ≈ .5ns * d^2
    # Max input: 2^(NN-2)-1 for IntNN
    b = Int(round(sqrt(2n)))
    a = n - (b*(b-1))>>1
    #Depending on benchmark,
    # >>1 makes unrank2d 1.1-5x faster than ÷2.
    b,a
end
function unrank3d(n)
    # O(1) for n :: Int64.
    # Constant: 80ns ≈ 50x unrank2d ≈ 5x memory read or write ≈ 9ns * d^2
    # Max input: ??? (>= 10^8)
    c = Int(floor((6n+(6n)^(1/3))^(1/3)))
    n -= (c-1)*c*(c+1)÷6
    b = Int(round(sqrt(2n)))
    a = n - (b*(b-1))>>1
    c,b,a
end
function unrank4d(n)
    # O(1) for n :: Int64.
    # Constant: 100ns ≈ 1.3x unrank3d ≈ 6x memory read or write ≈ 7ns * d^2
    # Max input: ??? (>= 10^8)
    d = Int(floor(sqrt(sqrt(24n+2.5sqrt(24n)))-.5))
    n -= (d-1)*d*(d+1)*(d+2)÷24
    c = Int(floor((6n+(6n)^(1/3))^(1/3)))
    n -= (c-1)*c*(c+1)÷6
    b = Int(round(sqrt(2n)))
    a = n - (b*(b-1))>>1
    d,c,b,a
end
function unrank5d(n)
    # O(1) for n :: Int64.
    # Constant: 200ns ≈ 2x unrank4d ≈ 11x memory read or write ≈ 8ns * d^2
    # Max input: ??? (>= 10^8)
    e1 = (120n)^(1/5)
    e = Int(floor((120n+5*e1^3+11*e1)^(1/5)-1))
    n -= (e-1)*e*(e+1)*(e+2)*(e+3)÷120
    d = Int(floor(sqrt(sqrt(24n+2.5sqrt(24n)))-.5))
    n -= (d-1)*d*(d+1)*(d+2)÷24
    c = Int(floor((6n+(6n)^(1/3))^(1/3)))
    n -= (c-1)*c*(c+1)÷6
    b = Int(round(sqrt(2n)))
    a = n - (b*(b-1))>>1
    e,d,c,b,a
end

function swap(array,i,j)
    temp = array[i]
    array[i] = array[j]
    array[j] = temp
end
function unrank_permutation!(array,d,n)
    #Algorithm source:
    #Wendy Myrvold & Frank Ruskey
    #Ranking and unranking permutations in linear time (2001)
    #http://webhome.cs.uvic.ca/~ruskey/Publications/RankPerm/MyrvoldRuskey.pdf
    if d > 0
        swap(array, d, n%d+1)
        unrank_permutation!(array, d-1, n÷d)
    else
        array
    end
end
function unrank_permutation(d,n)
    unrank_permutation!(MVector{d,typeof(n)}(1:d),d,n)
end

function rootfactorial(d)
    out = 1.0#factorial(min(n,20))^(1/n)
    for i in 2:d
        out *= i^(1/d)
    end
    out
end
function inverse_binomial(d, n)
    #Returns the lowest a such that binomial(a,d) ≥ n
    a = max(d,Int(ceil(rootfactorial(d)*n^(1/d))))
    n_a = binomial(a,d)
    itterations = 0
    while n_a < n
        itterations += 1
        a += 1
        n_a *= a
        n_a ÷= a-d
    end
    @assert itterations < d
    a
end
function unrank_a(d, n)
    inverse_binomial(d,n)+1-d
end
function unrank!(array, d, n)
    array[d] = unrank_a(d,n)
    if d == 1
        return array
    end
    unrank!(array, d-1, n-binomial(array[d]+d-2,d))
end
function unrank(d, n)
    #O(d^2)
    #runtime ≈ 3,000ns + 15ns*d^2
    #1500x unrank2d, 10x-40x unrank3d-unrank5d
    unrank!(MVector{d}(ones(typeof(d),d)), d, n)
end

#A simple iterative approach for testing
start(d) = MVector{d}(ones(Int,d))
next!(vector) = next!(vector, length(vector))
function next!(vector, i)
    if i > 1 && vector[i] == vector[i-1]
        vector[i] = 1
        return next!(vector, i-1)
    else
        vector[i] += 1
        return vector
    end
end

#Testing
function test!(f,d,n=10^5)
    j = start(d)
    for i in 1:n
        try
            @assert f(i) == j
        catch exeption
            if i > 1
                println("f(",i-1,") = ",f(i-1)," = ",f(i-1))
            end
            println("f(",i,") = ",f(i)," ≠ ",j)
            throw(exeption)
        end
        next!(j)
    end
end

test!(vcat ∘ unrank1d,1)
test!(collect ∘ unrank2d,2)
test!(collect ∘ unrank3d,3)
test!(collect ∘ unrank4d,4)
test!(collect ∘ unrank5d,5)
test!(i -> reverse(unrank(3,i)),3,3*10^3)
test!(i -> reverse(unrank(17,i)),17,10^3)
println("This typically doens't print in atom ","Pass")

#= unrank benchmarks
using Random
using Plots
x = shuffle([1:9...,10:3:29...,30:20:90...,200])
t = [@belapsed unrank($(xi),4870923245) seconds=.1
    for xi in x]
scatter(x,t)
for (x,t) in zip(x,t)
    println(x,",",t*1e9)
end
=#

### unrank2d benchmarks, empirical runtime analysis, and old tests
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
