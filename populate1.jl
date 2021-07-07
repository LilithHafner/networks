include("unrank.jl")

###################################################################
#                                                                 #
#  The function generate_edges_in_cell generates a random         #
#    selection of valid edges in a specified region. It           #
#    supports rectangular regions of any rank.                    #
#                                                                 #
#                    Runtime = O(output size)                     #
#                                                                 #
#  R is the list of group sizes                                   #
#  m is the multi-index of the sub-tensor we are selectring from  #
#  p is the edge density                                          #
#                                                                 #
###                                                             ###
function generate_edges_in_cell(group_sizes, m, probability)
  nodes_from_each_group = [count(==(i), m) for i=1:length(group_sizes)]
  option_count = binomial.(group_sizes.+nodes_from_each_group.-1, nodes_from_each_group)
  Π_option_count = vcat([1],cumprod(option_count))
  Σ_group_sizes = vcat([0],cumsum(group_sizes))
  [vcat((unrank(nodes_from_each_group[i],
             (n-1)÷Π_option_count[i]%option_count[i]+1)
         .+ Σ_group_sizes[i] for i = 1:length(group_sizes))...)
   for n = grasshopper(last(Π_option_count),probability)]
end


###################################################################
#                                                                 #
#  grasshopper(n,p) is equivalent to filter(rand() < p, 1:n)      #
#     in output probability distribution, but runs faster.        #
#                                                                 #
#                       Runtime = O(output size)                  #
#               Average runtime = O(n*p)                          #
#                                                                 #
#  n is the size of the list we select from                       #
#  p is the expected fraction of elements to keep                 #
#                                                                 #
###                                                             ###
function grasshopper(n,probability)
  if probability*n == 0 return [] end
  @assert n > 0
  @assert probability > 0
  @assert probability <= 1
  out = Array{typeof(n)}(undef, min(n,ceil(Int,
    n*probability+(n*probability)^(3/4))))
  count = 0
  k = 1/log(1-probability)
  current = floor(k*log(1-rand()))+1
  while current <= n
    count += 1
    if count > length(out)
      resize!(out, min(n,length(out)*2))
    end
    out[count] = current
    current += floor(k*log(1-rand()))+1
  end
  resize!(out, count)
  out
end

group_sizes = [2,4,3]
m = [1,2,2,2]
probability = .2
display(generate_edges_in_cell(group_sizes, m, probability))
