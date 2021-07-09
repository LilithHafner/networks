# TODO
# create cohesive interface (perhaps export only a varient of unrank!)

#Small special case unranking
unrank0k(x) = nothing
unrank1k(x) = x # 0.05ns ≈ .05ns * k^2
function unrank2k(x)
    # O(1) for x :: Int64.
    # Constant: 1.8ns = sqrt runtime +5% -0% (VERY FAST) ≈ .5ns * k^2
    # Max input: 2^(NN-2)-1 for IntNN
    b = Int(round(sqrt(2x)))
    a = x - (b*(b-1))>>1
    #Depending on benchmark,
    # >>1 makes unrank2k 1.1-5x faster than ÷2.
    b,a
end
function unrank3k(x)
    # O(1) for x :: Int64.
    # Constant: 80ns ≈ 50x unrank2k ≈ 5x memory read or write ≈ 9ns * k^2
    # Max input: ??? (>= 10^8)
    c = Int(floor((6x+(6x)^(1/3))^(1/3)))
    x -= (c-1)*c*(c+1)÷6
    b = Int(round(sqrt(2x)))
    a = x - (b*(b-1))>>1
    c,b,a
end
function unrank4k(x)
    # O(1) for x :: Int64.
    # Constant: 100ns ≈ 1.3x unrank3k ≈ 6x memory read or write ≈ 7ns * k^2
    # Max input: ??? (>= 10^8)
    d = Int(floor(sqrt(sqrt(24x+2.5sqrt(24x)))-.5))
    x -= (d-1)*d*(d+1)*(d+2)÷24
    c = Int(floor((6x+(6x)^(1/3))^(1/3)))
    x -= (c-1)*c*(c+1)÷6
    b = Int(round(sqrt(2x)))
    a = x - (b*(b-1))>>1
    d,c,b,a
end
function unrank5k(x)
    # O(1) for x :: Int64.
    # Constant: 200ns ≈ 2x unrank4k ≈ 11x memory read or write ≈ 8ns * k^2
    # Max input: ??? (>= 10^8)
    e1 = (120x)^(1/5)
    e = Int(floor((120x+5*e1^3+11*e1)^(1/5)-1))
    x -= (e-1)*e*(e+1)*(e+2)*(e+3)÷120
    d = Int(floor(sqrt(sqrt(24x+2.5sqrt(24x)))-.5))
    x -= (d-1)*d*(d+1)*(d+2)÷24
    c = Int(floor((6x+(6x)^(1/3))^(1/3)))
    x -= (c-1)*c*(c+1)÷6
    b = Int(round(sqrt(2x)))
    a = x - (b*(b-1))>>1
    e,d,c,b,a
end

#General unranking
function rootfactorial(k)
    out = 1.0 #factorial(min(x,20))^(1/x) #This seems to make unrank slower
    for i in 2:k
        out *= i^(1/k)
    end
    out
end
function inverse_binomial(k, x)
    #Returns the lowest n such that binomial(n,k) ≥ x
    n = max(k,Integer(ceil(rootfactorial(k)*x^(1/k))))
    n_a = binomial(n,k)
    itterations = 0
    while n_a < x && itterations < k
        itterations += 1
        n += 1
        n_a *= n
        n_a ÷= n-k
        if itterations >= k
            throw(OverflowError)#TODO confirm that this is the only plausible cause
        end
    end
    n
end

function unrank!(array, k, x)
    if k == 0
        return array
    end
    array[k] = inverse_binomial(k,x)+1-k
    unrank!(array, k-1, x-binomial(array[k]+k-2,k))
end
function unrank(k, x)
    @assert k >= 0
    @assert x > 0
    #O(k^2)
    #for k < 10 runtime ≈ 65ns + 6ns * k + 28ns*k^2
    #for large k runtime ≈ 15ns * k^2
    #{2000, 120, 4, 5, 4}x the runtime of unrank{1:5}k
    unrank!(ones(typeof(x),k), k, x)
    #unrank!(MVector{k}(ones(typeof(x),k)), k, x)
    # (This would require k to be const and using StaticArrays)
end

function unrank_lex!(array, k, x, limit, i, v)
    #sets the entries of array[i:i+k] to
    #  the lexicographically xth element of [v:limit]^k
    if k == 0
        return array
    end
    if k == 1#This special case gives ~2x speedup
        array[i] = v + x - 1
        return array
    end#TODO extend this?
    total = binomial(limit+k-v,k)
    ib = inverse_binomial(k,total-x+1)
    array[i] = limit-ib+k
    unrank_lex!(array, k-1, x-total+binomial(ib,k), limit, i+1, array[i])
end
unrank_lex!(array, k, x, limit, i) = unrank_lex!(array, k, x, limit, i, array[i])
unrank_lex!(array, k, x, limit)    = unrank_lex!(array, k, x, limit, 1, 1)
unrank_lex(k, args...) = unrank_lex!(ones(Int, k), k, args...)

#A simple iterative approach for testing
start(k) = ones(Int,k)
function start0(k)
    out = start(k)
    out[k] -= 1
    out
end
next!(vector) = next!(vector, length(vector))
function next!(vector, k)
    if k > 1 && vector[k] == vector[k-1]
        vector[k] = 1
        return next!(vector, k-1)
    else
        vector[k] += 1
        return vector
    end
end
#TODO consistant style w.r.t. order of main method and aliases
next_lex!(vector, limit) = next_lex!(vector, limit, length(vector))
function next_lex!(vector, limit, k)
    if vector[k] == limit
        if k == 1
            throw(OverflowError("next_lex!("*
                string(vector)*", "*string(limit)*") overflows"))
        end
        next_lex!(vector, limit, k-1)
        vector[k] = vector[k-1]
    else
        vector[k] += 1
    end
    vector
end

#Testing
function test!(next!,unrank,k,max_x=10^5)
    i = start0(k)
    for x in 1:max_x
        next!(i)
        try
            @assert unrank(x) == i
        catch
            if x > 1
                println("unrank(",x-1,") = ",unrank(x-1)," = ",unrank(x-1))
            end
            println("unrank(",x,") = ",unrank(x)," ≠ ",i)
            rethrow()
        end
    end
end
function quicktest!()
    test!(next!, vcat ∘ unrank1k,1,10^4)
    test!(next!, collect ∘ unrank2k,2,10^4)
    test!(next!, collect ∘ unrank3k,3,10^4)
    test!(next!, collect ∘ unrank4k,4,10^4)
    test!(next!, collect ∘ unrank5k,5,10^4)

    test!(next!, x -> reverse(unrank(1,x)),1,3*10^2)
    test!(next!, x -> reverse(unrank(3,x)),3,3*10^3)
    test!(next!, x -> reverse(unrank(17,x)),17,10^3)

    test!(i -> next_lex!(i,300), x -> unrank_lex(1,x,300),1,300)
    test!(i -> next_lex!(i,30), x -> unrank_lex(3,x,30),3,3*10^3)
    test!(i -> next_lex!(i,5), x -> unrank_lex(17,x,5),17,10^3)

    #Confirm that unrank_lex and unrank_lex! are consistant
    for size in [3,10,100,1000]
        for _ in 1:100
            start = rand(-size:size,size)
            i = rand(1:size)
            v = rand(-size:size)
            range = rand(1:size)
            limit = v + range - 1
            max_k = Int(floor(log(typemax(Int))/(log(range)+1)))#Keep limit^k < typemax(Int)
            k = rand(1:min(max_k, size-i+1))
            max_x = binomial(range+k-1, k)
            for x in [1, max_x, rand(1:max_x,5)...]
                copy1 = unrank_lex!(Vector(start), k, x, limit, i, v)
                copy2 = Vector(start)
                copy2[i:i+k-1] = unrank_lex(k, x, range) .+ (v-1)
                try
                    @assert copy1 == copy2
                catch
                    println("unrank_lex!'s copy = ",copy1," ≠ ")
                    println("unrank_lex's copy  = ",copy2)
                    println("for parameter values:")
                    println("start = ",start)
                    println("    i = ",i)
                    println("    k = ",k)
                    println("    v = ",v)
                    println("limit = ",limit," (range = ",range,")")
                    println("    x = ",x)
                    rethrow()
                end
            end
        end
    end

    #println("unrank passed quicktest")
end
quicktest!()#This adds ~50ms to startup times

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
