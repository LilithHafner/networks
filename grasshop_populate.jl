include("unrank.jl")


#TODO support total node count up to at least 2^30, ideally 2^60
# hyper_stochastic_block_edge_count(3,fill(10000,10), 10^5); works
# hyper_stochastic_block_edge_count(3,fill(100000,10), 10^5); hangs
# Theory: prod(BigInt.(group_sizes[m])) > typemax(Int)
#   is too generous an upper limit
# Ideally support seamless transition to bigint for larger
#TODO identify bottleneck for high dimension (e.g. >5)
#TODO benchmark and speed up constant factors
#   Possible 2x~5x speedup for low dimensions
#   Unknkown bounds for higher dimensions.
#TODO test increasing dimension with constant probability
function grasshop_populate!(edges, group_sizes, m, probability)
    if probability == 0 return edges end
    # WARNING this is not quite strict enough:
    if prod(BigInt.(group_sizes[m])) > typemax(Int)
        throw(OverflowError("Entire region must fit in an Int"))
    end
    if probability < 0 || probability > 1
        throw(ArgumentError("edge probability "*string(probability)*" ∉ [0,1]"))
    end
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
        #Should we special case 2D as well? (Guess: not here, but TODO in unrank_lex!)

        td = 1 #trailing dimensions
        for i in length(m):-1:1#We start from the end for lexegraphic order
            room = binomial(upper_limit[i]-current_edge[i]+td, td)
            #if we can finish progressing to the next edge placement, finish.
            if room > increment
                unrank_lex!(current_edge, td, increment+1, upper_limit[i], i)
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
                unrank_lex!(current_edge, td, increment%all_room+1, upper_limit[i], i, lower_limit[i])
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
grasshop_populate(args...) = grasshop_populate!([], args...)

@assert grasshop_populate([2,3,7], [1,2,2], 1) ==
   [[1, 3, 3], [1, 3, 4], [1, 3, 5], [1, 4, 4], [1, 4, 5], [1, 5, 5],
    [2, 3, 3], [2, 3, 4], [2, 3, 5], [2, 4, 4], [2, 4, 5], [2, 5, 5]]
