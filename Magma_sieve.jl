
##
Fachpraktikum = UInt32(1) #v0x00.000.001
##

#WARNING all implementations are quick and naive (slow but mostly correct i hope)
# improvement will be done later
#TODO insert abstract type control in all functions
##
using Combinatorics,Hecke,Revise,Nemo
revise()
##

function erat(b::Int)                                                             #TODO output should be fmpz
    #basic Sieve of Erastosthenes
    L1 = [2] # only even prime
    L2 = [i for i in 3:2:b]  # only get uneven numbers  -> use O(n/2) memory
    while length(L2) != 0
        nextp = L2[1]
        L1 = vcat(L1,[nextp])   #add p according to Erastosthenes
        for i in L2
            if i%nextp == 0
                k = findfirst(isequal(i),L2)
                deleteat!(L2,k)
            end
        end
    end
    return L1
end

function issmoth_v(n::fmpz,B::Array)
    # for p in B do. n = n/p as long as possible.
    # remaining  n is 1 -> true
    # if remaining n != 1  -> false
    r = deepcopy(n)
    r != 0 || error("zero_error")
    for p in B
        while r%p == 0
            r = divexact(r,p)
        end
        r != 1 || return true
    end
    return false
end

function pm1(n,B)
    # pollard p-1 methode
    a = 2
    for p in my_primesupto(Int(B))
        prod = p
        while prod < B
            a = powermod(a,p,n)
            prod  = prod * p # such that p^i is b smooth
        end
        g = gcd(a-1,n)
        if 1 < gcd(a-1,n) < n
            return (g,true)
        end
    end
    return (0,false)
    #error("failed pollard p-1 ")
end

function prho(n,f)
    (x,z,d) = (1,1,1)
    while true
        x = evaluate(f,x)
        z = evaluate(f,evaluate(f,z))
        d = gcd(lift(z-x),n)
        d != n || return (0,false)
        abs(d) == 1 || return (d,true)
    end
end

function prhofactor(n)                                                            #TODO check this function with prho again
    # this should give all prime factors.           #TODO vielfachheiten der Primzahlen
    Zy,y = PolynomialRing(ResidueRing(ZZ,n),"y")
    f = y^2 + 1
    f1,bool = prho(n,f)    # if bool == false n should be prime
    bool || return [n]
    f2 = divexact(n,f1)
    return vcat(prhofactor(f1),prhofactor(f2))
end

function alldivisors(n)
    pdivisors = prhofactor(n)
    A = [1]
    for vector in powerset(pdivisors,1)# gives iterator and filter emptyset
        factor = prod(vector)
        (factor in A) || (A = vcat(A,factor))
    end
    return A
end

function primitive_elem(K,first)
    #returns a primitive element s.t. lift(a) is prime in Z
    #additional first==true return first primitive element  prime in Z (ordered by its lifting in Z)
    p = Int(length(K))
    Fact = alldivisors(p-1)[1:end-1] # factor #K-1 = p-1 without p-1
    if first
        while true
            for y in K
                A = [y^exp for exp in Fact]
                one(K) in A  || ((isprime(lift(y)) ? (return y) : nothing))
            end
        end
    end
    # random primitive elem
    while true
        y = rand(K)
        A = [y^exp for exp in Fact]
        one(K) in A  || ((isprime(lift(y)) ? (return y) : nothing))
    end
end

function verify_primitive_elem(elem,K)
    p = Int(length(K))
    for i = 1:p-2
        elem^i != 1 || return false
    end
    return true
end

##
# implementation of the Magma_Sieve
##

function Sieve(K, qlimit, climit, ratio)
    p = characteristic(K) #(p = Int(length(K)))
    H = floor(root(p,2))+1
    J = H^2 - p

    # factorbase up to qlimit
    fb_primes = erat(qlimit)
                                                                                    # TODO faster way here
    log2 = log(2.0);
    logqs = [log(q)/log2 for q in fb_primes] #real logarithms for sieve

    alpha = primitive_elem(GF(103),true) # first primitive elem, prime in Z i.e small
    FB = vcat([lift(alpha)],deleteat!(fb_primes,findfirst(isequal(lift(alpha)),fb_primes)))
    FBs = FactorBase(FB)
    #FB = Factorbase([fmpz(i) for i in FB)
    #fb_primes is factorbase with lift of alpha on first place
    #A = zeros(Int64,length(fb_primes),length(fb_primes))
    A = sparse_matrix(ZZ) #TODO

    for c1 = 1:climit
        nrows(A)/length(FB) < ratio || break
        sieve = zeros(climit)
        den = H + c1;                # denominator of relation
        num = -(J + c1*H)            # numerator
        for i=1:length(fb_primes)
            q = fb_primes[i]
            qpow = q
            nextqpow = qpow
            logq = logqs[i]
            while qpow < qlimit      # qlimit-smooth
                den % qpow != 0 || break
                c2 = num * invmod(den, fmpz(qpow))  % qpow
                (c2 != 0) || (c2 = qpow)
                nextqpow = qpow*q    #increase q_power
                while c2 < c1   #c2>=c1 to remove redundant relations of (c1,c2) tuples, just increase c2
                    c2+=qpow
                end
                while c2 < length(sieve)
                    sieve[Int(c2)] += logq
                    if nextqpow > qlimit
                        prod = (J + (c1 + c2)*H + c1*c2) % p
                        nextp = nextqpow
                        while rem(prod,nextp) == 0
                            sieve[c2] += logq
                            nextp = nextp*q
                        end
                    end
                    c2 += qpow
                end
                qpow = nextqpow
            end
        end
        #include relations / check sieve for full factorizations.
        rel = den * (H + 1)
        relinc = H + c1       # add to relation to get next relation
        idx = 0
        for c2 in 1:length(sieve)
            n = rel % p
            if abs(sieve[c2] - floor(log(n)/log2)) < 1
                # TODO insert default factorbase algorithm
                #FBs = FactorBase(FB) #generate Factorbas from updated FBs with new c_i´s
                if issmooth(FBs,fmpz(n))
                    dict_factors = Hecke.factor(FBs,fmpz(n))
                    #Include each H + c_i in extended factor basis.
                    (H + c1) in FB || (FB = vcat(FB,[H + c1]))
                    (H + c2) in FB || (FB = vcat(FB,[H + c2]))
                    #Include relation (H + c1)*(H + c2) = fact.
                    row = nrows(A) + 1 # insert new relation (new row)to sparse_matrix
                    v = values(dict_factors)
                    cnt = 0
                    J1 = Vector{Int}([])
                    V = []
                    for entry in dict_factors
                        cnt+=1
                        (entry[2] == 0) || ((J1,V)=(vcat(J1,[cnt]),vcat(V,entry[2])))
                    end
                    if c1 == c2
                        h = findfirst(isequal(H + c1),FB)
                        (J1,V)=(vcat(J1,[h]),vcat(V,fmpz(-2)))
                    else
                        h = findfirst(isequal(H + c1),FB)
                        (J1,V)=(vcat(J1,[h]),vcat(V,fmpz(-1)))
                        h = findfirst(isequal(H + c2),FB)
                        (J1,V)=(vcat(J1,[h]),vcat(V,fmpz(-1)))
                    end
                    #println(J1,V)
                    push!(A,sparse_row(ZZ, J1, V))
                    #setindex!(A, new_sparserow, row) #TODO doesnt work as wanted.
                end
            end
            rel += relinc
        end
    end
    return A,FBs,FB
end

K = Sieve(GF(103),35,27,1.1)
println(K[3])

# returns a 30x69 matrix
# in FB still are redundant primes... ?!
