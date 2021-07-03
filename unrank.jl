function unrank2d(n)
    # Roughly O(1) for n :: Int64.
    # Constant: 2ns
    a = Int(round(sqrt(2n+1)-1))
    b = n-(a*(a+1))>>1
    a,b
end

#println(collect(unrank2d.(0:10)))
#@benchmark unrank2d(9720290439024)
