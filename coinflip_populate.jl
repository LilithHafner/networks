nondecreasing(e) = all(map(<=, e[1:length(e)-1], e[2:length(e)]))

coinflip_populate!(edges, group_sizes, m, probability) =
    append!(edges, coinflip_populate(group_sizes, m, probability))
function coinflip_populate(group_sizes, m, probability)
    offset = (cumsum(group_sizes) .- group_sizes)[m]
    filter(nondecreasing, [offset.+Tuple(i) for i in
        findall(rand(group_sizes[m]...) .< probability)])
end

## Tests
@assert nondecreasing(1:10)         == true
@assert nondecreasing([1,3,2])      == false
@assert nondecreasing([])           == true
@assert nondecreasing((5,2))        == false
@assert nondecreasing((-2,-2,5,5,5))== true
@assert nondecreasing(10:-3:2 )     == false

@assert sort(coinflip_populate([2,3,7], [1,2,2], 1)) ==
   [[1, 3, 3], [1, 3, 4], [1, 3, 5], [1, 4, 4], [1, 4, 5], [1, 5, 5],
    [2, 3, 3], [2, 3, 4], [2, 3, 5], [2, 4, 4], [2, 4, 5], [2, 5, 5]]
