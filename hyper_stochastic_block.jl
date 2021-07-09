function hyper_stochastic_block(populate!, k, group_sizes, probability_function; edges=:auto)
    if edges == :auto
        edges = []
    end
    groups = length(group_sizes)
    m = start0(k)
    for i in 1:binomial(length(group_sizes)+k-1,k)
        next_lex!(m, groups)
        populate!(edges, group_sizes, m, probability_function(m))
    end
    return edges
end
function hyper_stochastic_block(populate!, k, group_sizes, probability::Real; kw...)
     hyper_stochastic_block(populate!, k, group_sizes, m->probability; kw...)
end
function hyper_stochastic_block_edge_count(populate!, k, group_sizes, edge_count::Integer; kw...)
     hyper_stochastic_block(populate!, k, group_sizes,
        min(1,edge_count/binomial(BigInt(sum(group_sizes)+k-1),k)); kw...)
end
#display(hyper_stochastic_block(3,[10,10,10],a -> a[1]== a[2]==a[3] ? .01 : .001))
