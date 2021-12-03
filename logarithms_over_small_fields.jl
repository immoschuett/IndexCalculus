#We here define functions that aim to find logarithms in cyclic groups of smaller orders.
#Solve g^x=a for x in the cyclic group G=<g>

using Hecke

function bsgs(G, g, a) #baby-step giant-step algorithm
    n = fmpz(length(G))
    m = ceil(root(n, 2))
    g_inv = inv(g)
    B = [(a*g_inv^fmpz(r), fmpz(r)) for r = 0:(Int64(m)-1)] #baby-steps
    for (i,j) in B
        if i == one(G)   
            @debug g^j == a ? true : @error("logarithm wrong")           # a*g^(-r) = 1 -> a = g^r
            return j
        end
    end
    d = g^m
    q = fmpz(1)
    while true
        @debug q <= n ? true : @error("no logarithm found") 
        for (k,l) in B
            if d^q == k                                              # giant-steps d^q
                @debug g^(q*m+l) == a ? true : @error("logarithm wrong(2)")
                return q*m + l
            end
        end
        q+=1
    end
end


Z2017 = ResidueRing(ZZ, 2017)
g = Z2017(5)
a = Z2017(3)
bsgs(Z2017, g, a)  #solution should be x = 22*45+40 = 1030

function test()
    for i = 1:1000
        randomnumber = Z2017(rand(2:2016))
        ans = bsgs(Z2017, g, randomnumber)
        if g^ans != randomnumber
            return false 
        end 
    end 
    return true 
end 

test()

ENV["JULIA_DEBUG"] = "" # enable debugging