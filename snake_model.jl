function snake(collection::AbstractVector{T}, size::Integer,
               width::Integer, speed::Real) where T <: Any
    #O(size) constant factor = 12.5 ns for
    #    collection size 1–typemax(Int), size 100–10^9, width 1–typemax(Int), speed -Inf–Inf

    out = Vector{T}(undef, size)
    k0,k1 = length(collection)/speed,width/(1-speed)#*

    for i in 1:min(width, size)
        out[i] = rand(collection)
    end
    for i in width+1:size
        r = rand()#*
        out[i] = r < speed ? collection[T(floor(r*k0)+1)] : out[i-1-T(floor((r-speed)*k1))]#*
    #    out[i] = rand() < speed ? rand(collection) : out[rand(i-width:i-1)] #This simplification costs a 13–20% increase in runtime
    end
    out
end

function snakemodel1(hyperdegree::Integer, nodes::Integer, edges::Integer,
                     width::Integer, speed::Real, duplicates::Bool)
    source = snake(1:nodes, hyperdegree*edges, width, speed)
    if duplicates
        [source[hyperdegree*(e-1)+1:hyperdegree*e] for e in 1:edges]
    else
        collect.(collect(Set(Set(source[hyperdegree*(e-1)+1:hyperdegree*e])) for e in 1:edges))
    end
end

#=  Meta TODO create a testng infrastructure to automatically run tests?
        This would be more work to maintain
        keep a library of functional examples availible
        result in a higher codebase cleanliness?
        reduce the intruduction of new bugs while modifying
            This is the primary reason to have testing.
            I also don't think it is fully applicible here.
        Temporary answer: no.
=#

#=
# Mildly out of date:
#Basic inspection
display(transpose(reshape(snake(10:99, 40, 5, .5), 5, 8)))
display(snakemodel1(5, 99, 8, 5, .5))
display(snakemodel1(3, 3, 60, 5, .5))

#Benchmarking
using BenchmarkTools
mean((@belapsed snake(1:Int(10)^18, $(Int(floor(10^i))), 200, .5) seconds=1)/Int(floor(10^i))*1e9 for i in 1:5)

#Displaying graphs
using LightGraphs
using GraphPlot
function graph_plot(node_count, edges)
    graph = SimpleGraph(node_count)
    for (a,b) in edges
        add_edge!(graph, a,b)
    end
    gplot(graph)
end
graph_plot(60, snakemodel1(2, 60, 60, 5, .5))
=#

##
function distribute(distribution::Vector{I}) where I <: Integer
    i = zero(I)
    out = Vector{I}(undef, sum(distribution))
    for (j,count) ∈ enumerate(distribution)
        out[i+1 : (i+=count)] .= j
    end
    out
end
using Random: shuffle!
function snakemodel2(hyperdegree_distribution::Vector{I},
                          degree_distribution::Vector{I},
                                        width::Integer,
                                        speed::Real;
                                      maxsize::Any = 10^9,
         silently_correct_degree_distribution::Bool = false,
                                     verywide::Bool = false,
                             suppresswarnings::Bool = false,
                                             ) where I <: Integer
    #O(size) constant factor = 150~600 ns
    #Allows duplicates

    #Validate input
    size, degree_size = [sum(deg*count for (deg,count) in enumerate(dist))
            for dist in [hyperdegree_distribution, degree_distribution]]
    if !suppresswarnings
        if max(size, degree_size) > maxsize
            #Atom often omits the first argument to println. This is a hack for that.
            #       \/
            println("","WARNING: size = 10^"*string(round(10*log10(max(size, degree_size)))/10)*" > max_size = 10^"
                *string(round(10*log10(maxsize)/10))*"; this may take a while")
        end
        if !silently_correct_degree_distribution && size ≠ degree_size
            println("","WARNING: hyperdegree_distribution indicates "*string(size)
                *" node-edge interfaces while degree_distribution indicates "*string(degree_size)
                *"; actual degree distribution will be adjusted into alignment.")
        end
        if !verywide && width > 100*max(size,10)
            println("","WARNING: width = "*string(width)
                *" >> # of node-edge interfaces = "*string(size)
                *"; runtime = O(width).")
        end
        if speed < 0 || speed > 1
            println("", "WARNING: speed = "*string(speed)*" ∉ [0,1]")
        end
    end

    #Distribute distributions
    hyperdegrees = distribute(hyperdegree_distribution)
    degrees = distribute(degree_distribution)
    node_bucket = distribute(degrees)

    if size ≠ degree_size
        while length(node_bucket) < size
            node = rand(1:length(degrees))
            push!(node_bucket, node)
            degrees[node] += 1
        end
    end

    shuffle!(hyperdegrees)#Style: where to put these shuffles
    shuffle!(node_bucket)

    exclusions = fill(zero(I), sum(degree_distribution))#Tweak point
    past_i = zero(typeof(width))
    past = rand(node_bucket, width)#Tweak point
    k1, k2 = width/(1-speed), 1 - speed*width/(1-speed)
    function get_node()
        r = rand()
        while true
            if r >= speed
                node = past[I(floor(r*k1+k2))]
                if degrees[node] > 0
                    exclusions[node] += 1
                    degrees[node] -= 1
                    past_i = past_i%width+1
                    return past[past_i] = node
                end
            end
            node = pop!(node_bucket)
            if exclusions[node] > 0
                exclusions[node] -= 1
                r = r < speed ? r/speed : rand()
            else
                degrees[node] -= 1
                past_i = past_i%width+1
                return past[past_i] = node
            end
        end
    end

    [[get_node() for _ in 1:k] for k in hyperdegrees]
end

#Basic inspection
display(snakemodel2([1,3,2], [3,3,0,1], 2, .5))

#Benchmarking
using BenchmarkTools
function tm2(hdd, args...)
    t = @belapsed snakemodel2($hdd, $args...) seconds=1
    t/sum(distribute(hdd))*1e9, sum(distribute(hdd))
end

println("","(constant factor (ns), size)")
for i in 1:5
    println("",tm2([7,40,20+3*10^i,2*10^i,4*10^(i-1)],[0,20,19+3*10^i,10+2*10^i,2+4*10^(i-1)],20000,.1))
end

#Displaying graphs
using LightGraphs
using GraphPlot
function graph_plot(node_count, edges)
    graph = SimpleGraph(node_count)
    for (a,b) in edges
        add_edge!(graph, a,b)
    end
    gplot(graph)
end
graph_plot(60, snakemodel2([0,60], [20,15,10,10], 60, .5))
