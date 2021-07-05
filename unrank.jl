using Plots
using StaticArrays
using BenchmarkTools

unrank1d(n) = n
function unrank2d(n)
    # O(1) for n :: Int64.
    # Constant: 1.6ns = sqrt runtime +5% -0% (VERY FAST)
    # Max input: 2^(NN-2)-1 for IntNN
    b = Int(round(sqrt(2n)))
    a = n-(b*(b-1))>>1
    #Depending on benchmark,
    # >>1 makes unrank2d 1.1-5x faster than ÷2.
    b,a
end
function unrank3d(n)
    # O(1) for n :: Int64.
    # Constant: 80ns ≈ 50x unrank2d ≈ 5x memory read or write
    # Max input: ??? (>= 10^8)
    c = Int(floor((6n+(6n)^(1/3))^(1/3)))
    n -= (c-1)*c*(c+1)÷6
    b = Int(round(sqrt(2n)))
    a = n- (b*(b-1))>>1
    c,b,a
end
function unrank4d(n)
    # O(1) for n :: Int64.
    # Constant: 105ns ≈ 1.3x unrank3d ≈ 6x memory read or write
    # Max input: ??? (>= 10^8)
    d = Int(floor(sqrt(sqrt(24n+2.5sqrt(24n)))-.5))
    n -= (d-1)*d*(d+1)*(d+2)÷24
    c = Int(floor((6n+(6n)^(1/3))^(1/3)))
    n -= (c-1)*c*(c+1)÷6
    b = Int(round(sqrt(2n)))
    a = n- (b*(b-1))>>1
    d,c,b,a
end

function unrank(d,n)
    #=if d == 1
        return [n]
    elseif d == 2
        return unrank2d(n)
    elseif d == 3
        return unrank3d(n)
    end=#
    unrank!(zeros(MVector{d,typeof(n)}), d,n)
    #unrank!(Vector{typeof(n)}(undef,d), d,n)
end
function unrank!(array, d,n)
    if d == 1
        array[d] = n
        return array
    end
    if d == 2
        a = Int(round(sqrt(2n)))
        array[d] = a
        return unrank!(array, 1, n -(a*(a-1))>>1)
    end
    if d == 3
        a = Int(floor((6n)^(1/3)))
        array[d] = a
        return unrank!(array, 2, n - (a-1)*a*(a+1)÷6)
    end
    if d == 4
        a = Int(floor((24n+12)^(1/4)-.5))
        array[d] = a
        return unrank!(array, 3, n - (a-1)*a*(a+1)*(a+2)÷24)
    end
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
