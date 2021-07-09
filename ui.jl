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
probability = .02# <- Notice that constants are allowed as well as functions
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
edge_count = 20# <- Specify edges
graph3 = hyper_stochastic_block_edge_count(grasshop_populate!, hyperdegree, group_sizes, edge_count)
display(graph3)#               ^^^^^^^^^^^
#                                  /\
# I was not well aquainted with the finer points of multiple dispatch when I
# wrote this, so I was hesitant to trust it. Hence the name change.

# Have more questions? Write to me at <hafnerda@grinnell.edu>!
