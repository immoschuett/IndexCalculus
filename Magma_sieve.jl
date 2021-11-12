
##
Fachpraktikum = UInt32(1) #v0x00.000.002
##

#WARNING all implementations are quick and naive (slow but mostly correct i hope)
# improvement will be done later
#TODO insert abstract type control in all functions
##

using Hecke,Nemo,Revise
ENV["JULIA_DEBUG"] = "all" # enable debugging
revise()
##
function primitive_elem(K::Nemo.GaloisField,first::Bool) #TODO implement compatible for Abstract field or Nemo.GaloisfmpzField
    #returns a (the first) generator alpha of K° s.t. lift(alpha) is prime in ZZ
    p = length(K)
    Fact = divisors(fmpz(p-1))[1:end-1]
    while true # alpha exists
        for y in K
            if !first y = rand(K) end
            A = [y^exp for exp in Fact]
            one(K) in A  || ((isprime(lift(y)) ? (return y) : nothing))
        end
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
    fb_primes = [i for i in PrimesSet(1,qlimit)]
    log2 = log(2.0);
    logqs = [log(q)/log2 for q in fb_primes] #real logarithms for sieve

    alpha = primitive_elem(K,true) # first primitive elem, prime in Z i.e small
    FB = vcat([lift(alpha)],deleteat!(fb_primes,findfirst(isequal(lift(alpha)),fb_primes))) # tauschen a[1] = a[2] , a[2] = [1]
    FBs = deepcopy(FactorBase(FB))
    l = length(FB)
    Indx = Dict(zip(FB,[i for i=1:length(FB)])) #Index in a dictionary
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
                while c2 <= length(sieve)
                    sieve[Int(c2)] += logq
                    if nextqpow > qlimit
                        prod = (J + (c1 + c2)*H + c1*c2) % p
                        nextp = nextqpow
                        while rem(prod,nextp) == 0
                            sieve[Int(c2)] += logq
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
                    r = length(Indx)+1
                    if !((H + c1) in keys(Indx))
                        push!(FB,H + c1)
                        push!(Indx,(H+c1) => r)
                    end#(FB = vcat(FB,[H + c1])) #push!(FB,wert)
                    r = length(Indx)+1
                    if !((H + c2) in keys(Indx))
                        push!(FB,H + c2)
                        push!(Indx,(H+c2) => r)
                    end#(FB = vcat(FB,[H + c2]))
                    #Include relation (H + c1)*(H + c2) = fact.
                    #row = nrows(A) + 1 # insert new relation (new row)to sparse_matrix
                    J1 = Vector{Int}([])
                    V = Vector{fmpz}([])
                    for (prime,power) in dict_factors
                        if !(power == 0)
                            push!(J1,Indx[prime])
                            push!(V,power)
                        end
                    end
                    if c1 == c2
                         push!(J1,Indx[H+c1])
                         push!(V,fmpz(-2))
                    else
                         push!(J1,Indx[H+c1])
                         push!(J1,Indx[H+c2])
                         push!(V,fmpz(-1))
                         push!(V,fmpz(-1))
                    end
                    push!(A,sparse_row(ZZ, J1, V))
                end
            end
            rel += relinc
        end
    end
    @debug @assert check_relation_mat(K,A,FB) "Relation Matrix wrong"
    return A,FBs,FB
end

function check_relation_mat(K,A,FB)
    #getindex(A::SMat, i::Int, j::Int)
    p = Int(length(K))
    for k = 1:nrows(A)
        prod = 1
        for i = 1:length(FB)
            ex = getindex(A,k,i)
            if ex < 0
                prod *= invmod(FB[i],fmpz(p))^abs(ex)
            else
                prod  *= FB[i]^ex
            end
        end
        if !(prod % p  == 1)  return false end
    end
    return true
end



K,K2,K3 = Sieve(GF(107),35,28,1.1)

#Hecke.subset
#flog,clog,iroot

#badge smoothneth test -> glatttest (factorisieren)
#psi lower, psi upper in Hecke (Schranken), rho-funktion (dickmann_rho)

G = GF(1000000007)
A,B,C = Sieve(G, 250, 80, 1.1)
