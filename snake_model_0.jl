function hyper_snake(hyperdegree, length, nodes, bp1, bpn)
    edges = Vector{Vector{typeof(nodes)}}(undef, length)
    bpo = bp1
    k = 1/log(1-bpn)
    for i in 1:length
        edges[i] = [begin
                r = rand();
                rand(if r < bp1-bpo#Forward or backward

                    edges[i-1-typeof(nodes)(floor(k*log(1-r/bp1)))]
                    #r = 0 → i-1-I(floor(k*log(1))) = i-1
                    #r = bp1-bpo → i-1-I(floor(k*log(1-(bp1-bpo)/bp1))) =
                    # = i-1-I(floor(log(bpo/bp1)/log(1-bpn))) =
                    # = i-1-I(floor(log((1-bpn)^(i-1))/log(1-bpn))) =
                    # = i-1-I(floor(i-1)) = 0
                    # R = [0,pb1-bp0) → [i-1, 0) → [i-1, 1]
                else
                    1:nodes#Insert degree correction here
                end)
            end for _ in 1:hyperdegree]
        bpo *= 1-bpn
    end
    edges
end
display(hyper_snake(4, 8, 100, .9, .9))
