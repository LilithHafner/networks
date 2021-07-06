# TODO
# alpha renaming to conform to x = binomial(n,k): d -> k, n -> x, a -> n.
# create cohesive interface (perhaps export only a varient of unrank!)

#Small special case unranking
unrank0d(n) = nothing
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

#General unranking
function rootfactorial(d)
    out = 1.0 #factorial(min(n,20))^(1/n) #This seems to make unrank slower
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
function unrank!(array, d, n)
    if d == 0
        return array
    end
    array[d] = inverse_binomial(d,n)+1-d
    unrank!(array, d-1, n-binomial(array[d]+d-2,d))
end
function unrankr!(array, d, n, limit, i, v)
    if d == 0
        return array
    end
    if d == 1#This special case gives ~2x speedup
        array[i] = v + n - 1
        return array
    end#TODO extend this?
    total = binomial(limit+d-v,d)
    ib = inverse_binomial(d,total-n+1)
    array[i] = limit-ib+d
    unrankr!(array, d-1, n-total+binomial(ib,d), limit, i+1, array[i])
end
unrankr!(array, d, n, limit, i) = unrankr!(array, d, n, limit, i, array[i])
unrankr!(array, d, n, limit)    = unrankr!(array, d, n, limit, 1, 1)

function unrank(d, n)
    @assert d >= 0
    @assert n > 0
    #O(d^2)
    #for d < 10 runtime ≈ 65ns + 6ns * d + 28ns*d^2
    #for large d runtime ≈ 15ns * d^2
    #{2000, 120, 4, 5, 4}x the runtime of unrank{1:5}d
    unrank!(ones(typeof(n),d), d, n)
    #unrank!(MVector{d}(ones(typeof(n),d)), d, n)
    # (This would require d to be const and using StaticArrays)
end

#A simple iterative approach for testing
start(d) = ones(Int,d)
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
function quicktest!()
    test!(vcat ∘ unrank1d,1)
    test!(collect ∘ unrank2d,2)
    test!(collect ∘ unrank3d,3)
    test!(collect ∘ unrank4d,4)
    test!(collect ∘ unrank5d,5)
    test!(i -> reverse(unrank(1,i)),1,3*10^2)
    test!(i -> reverse(unrank(3,i)),3,3*10^3)
    test!(i -> reverse(unrank(17,i)),17,10^3)

    println("unrank passed quicktest (unrankr! not tested)")
end

#= general unrank benchmarks
using BenchmarkTools
using Random
using Plots

x = shuffle([1:9...,10:3:29...,30:20:90...,200])
t1 = [@belapsed unrank($(xi), 4870923245) seconds=.1
    for xi in x]
t2 = [@belapsed unrank!($(ones(Int,xi)), $(xi), 4870923245) seconds=.1
    for xi in x]
scatter(x,t1)
scatter!(x,t2)
for (x,t1,t2) in zip(x,t1,t2)
    println(x,",",t1*1e9,",",t2*1e9)
end
=#
