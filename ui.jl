include("coinflip_populate.jl")
include("grasshop_populate.jl")
include("hyper_stochastic_block.jl")


#This is the classic way to opperate the model:
hyperdegree = 3#Nodes per edge
group_sizes = [3,2]#The sum of groups sizes is the total number of nodes
#This is a function that takes in a block index (which is a nondecreasing
#sequence of hyperdegree integers between 1 and the number of groups) and returns the
#probability with which each edge in that block will be generated
probability = (m) -> all(m[1] .== m) ? .9 : .1
graph1 = hyper_stochastic_block(grasshop_populate!, hyperdegree, group_sizes, probability)
display(graph1)


hyperdegree = 6
group_sizes = [3,2]
probability = .02# <- Notice that `Real`s are allowed as well as `Function`s
graph2 = hyper_stochastic_block(grasshop_populate!, hyperdegree, group_sizes, probability)
display(graph2)

#      Select population algorithm with this argument
#                                     \/
graph3 = hyper_stochastic_block(coinflip_populate!, hyperdegree, group_sizes, probability)
display(graph3)

# Want a certian number of edges but don't want to
# calculate the requisite probability? Here's how:
hyperdegree = 20
group_sizes = [4,4,4,4]
edge_count = 20# <- Specify edges as an Integer
graph3 = hyper_stochastic_block(grasshop_populate!, hyperdegree, group_sizes, edge_count)
display(graph3)

# Want lots of nodes with huge block sizes
# (i.e. on the order of or larger than typemax(Int))?
# Specify inputs as BitInts and BigFloats:
try
    local g4 = hyper_stochastic_block(grasshop_populate!, 6,Int[400,300000,400],6)
    display(g4)
catch OverflowError
    println("Got an overflow or inexact error? try providing input in a bigger type (e.g. BigInt & BigFloat)")
end

g4 = hyper_stochastic_block(grasshop_populate!, 6,BigInt[400,300000,400],6)
display(g4)

# Have more questions? Write to me at <hafnerda@grinnell.edu>!
