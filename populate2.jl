include("unrank.jl")

#unrank!(vector, index, dimensions, rank)
#increments the entries of vector starting at index and extending for dimensions
#entries to have a new local rank equal to rank relative to
#all values equal to vector[index] = rank of 1.


function populate!(edges, group_sizes, m, probability)
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
        increment = Int(floor(k*log(1-rand()))+1)

        #=TODO benchmark is this special case helpful?
        if last(current_edge)+increment <= upper_limit[i]
            current_edge[length(current_edge)] += increment
            push!(edges, vector(current_edge))
            continue
        end #TODO benchmark should we special case 2D as well?=#
        td = 1 #trailing dimensions
        i = 1
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

edges = []
group_sizes = [2,4,3]
m = [1,2,2,2]
probability = .999
display(populate!(e, group_sizes, m, probability))
