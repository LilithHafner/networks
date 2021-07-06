include("unrank.jl")

#unrank!(vector, index, dimensions, rank)
#increments the entries of vector starting at index and extending for dimensions
#entries to have a new local rank equal to rank relative to
#all values equal to vector[index] = rank of 1.

#TODO support total node count up to at least 2^30, ideally 2^60
#TODO benchmark and speed up. Right now ~14ns/byte. Possible 10x~100x speedup
#TODO support seamless transition to bigint for larger
function populate!(edges, group_sizes, m, probability)
    # WARNING entire region must fit in an Int
    @assert prod(BigInt.(group_sizes[m])) < typemax(Int)
    #O(output size)
    # ~ 500ns per edge
    # ~ 35ns per edge*dimension
    # ~ .6ns per bit @ 2^60
    # ~ 5ns per bit @ 100
    d = length(m)
    k = 1/log(1-probability)
    Σ_group_sizes = cumsum(group_sizes)
    upper_limit = Σ_group_sizes[m]
    prepend!(Σ_group_sizes, 0)
    lower_limit = Σ_group_sizes[m].+1
    current_edge = Vector(lower_limit)
    current_edge[length(current_edge)] -= 1 #Start before 1
    #=function limit(i)
        if i == d || m[i] != m[i+1]
            return limit[i] #TODO benchmark factor out?
        else
            return current_edge[i+1]
        end
    end=#

    while true
        increment = Int(floor(k*log(1-rand())))+1 #~20.ns 17%

        #current_edge[length(current_edge)] += 1 #1ns

        #= this special case can have a final runtime impact of
        #  ~ 15% faster — 5% slower
        if last(current_edge)+increment <= last(upper_limit)
            current_edge[length(current_edge)] += increment
            push!(edges, Vector(current_edge))
            continue
        end=#
        #Should we special case 2D as well? (Guess: not here, but in unrankr!)

        td = 1 #trailing dimensions
        for i in length(m):-1:1#We start from the end for lexegraphic order
            room = binomial(upper_limit[i]-current_edge[i]+td, td)
            #if we can finish progressing to the next edge placement, finish.
            if room > increment
                unrankr!(current_edge, td, increment+1, upper_limit[i], i)
                push!(edges, Vector(current_edge))
                break
            end

            #if we got to the very first digit and still can't finish, this edge
            #is out of bounds, and we return all accumulated edges
            if i == 1
                return edges
            end

            #if we're at a break between groups, unrank.
            if m[i-1] != m[i]
                #implicitly pull all subsequent nodes in the group to minimum
                #defined by group and compensate by increasing increment.
                all_room = binomial(upper_limit[i]-lower_limit[i]+td,td)
                increment += all_room-room
                #allocate the increment between this group and previous
                unrankr!(current_edge, td, increment%all_room+1, upper_limit[i], i, lower_limit[i])
                increment ÷= all_room
                td = 1
            else
                #implicitly pull all subsequent nodes in the group to minimum
                #defined by previouss value and compensate by increasing increment.
                increment += binomial(upper_limit[i]-current_edge[i-1]+td,td)-room
                #quietly continue
                td += 1
            end
            i -= 1
        end
    end
end
populate(args...) = populate!([], args...)

edges = []
group_sizes = [2,4,3]
m = [1,2,2,2]
probability = .999
display(populate!(edges, group_sizes, m, probability))

using BenchmarkTools
function time_per_edge_probability(args... ; n=100)
    lengths = []
    times = [@elapsed push!(lengths,length(populate(args...)))
            for i in 1:n]
    t,l = minimum(zip(times./lengths,lengths))
    #l = mean(lengths)
    l, t*1e9
end
function time_per_edge(group_sizes, m, edges ; n=10)
    size = prod(group_sizes[m])/prod(factorial(count(x->x==i, m)) for i in 1:length(group_sizes))
    probability = edges/size
    lengths = []
    times = [@elapsed push!(lengths,length(populate(group_sizes, m, probability)))
            for i in 1:n]
    t,l = minimum(zip(times./lengths,lengths))
    #l = mean(lengths)
    l, t*1e9
end
