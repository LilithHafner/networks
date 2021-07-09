#TODO? utilize args... IDK if this works as well with Atom's editor hints, tho.
function hyper_stochastic_block!(edges, populate!, k, group_sizes, probability::Function)
    groups = length(group_sizes)
    m = start0(k)
    for i in 1:binomial(length(group_sizes)+k-1,k)
        next_lex!(m, groups)
        populate!(edges, group_sizes, m, probability(m))
    end
    edges
end
function hyper_stochastic_block(populate!, k, group_sizes, probability::Function)
    hyper_stochastic_block!(typeof(group_sizes)[], populate!, k, group_sizes, probability)
end
function hyper_stochastic_block(populate!, k, group_sizes, probability::Real; kw...)
     hyper_stochastic_block(populate!, k, group_sizes, m->probability; kw...)
end
function hyper_stochastic_block(populate!, k, group_sizes, edge_count::Integer; kw...)
     hyper_stochastic_block(populate!, k, group_sizes,
        min(1,edge_count/binomial(sum(group_sizes)+k-1,k)); kw...)
end

#include("grasshop_populate.jl")
#include("coinflip_populate.jl")
#display(hyper_stochastic_block(grasshop_populate!, 3,[10,10,10],a -> a[1]== a[2]==a[3] ? .01 : .001))
#display(hyper_stochastic_block(coinflip_populate!, 3,[10,10,10],a -> a[1]== a[2]==a[3] ? .01 : .001))
