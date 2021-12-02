
##
Fachpraktikum = UInt32(2) #v0x00.000.002
##

#WARNING all implementations are quick and naive (slow but mostly correct i hope)
# improvement will be done later
#TODO insert abstract type control in all functions
##

using Hecke,Nemo,Revise
ENV["JULIA_DEBUG"] = "all" # enable debugging
revise()
##

abstract type DLP end

mutable struct FField<:DLP
	K::FinField
	a::FinFieldElem #primitive element
end
mutable struct Sparam<:DLP
	qlimit::Int64
	climit::Int64
	ratio::Float64
	inc::Tuple{Int64, Int64} #for increasing of Sparam
end


##
function primitive_elem(K::FinField,first::Bool) #TODO implement compatible for Abstract field or Nemo.GaloisfmpzField
    #returns a (the first) generator alpha of K° s.t. lift(alpha) is prime in ZZ
    p = length(K)
    Fact = prime_divisors(fmpz(p-1))
    while true # alpha exists
        for y in K
            if !first y = rand(K) end
			if isprime(lift(y))
				if !(one(K) in [y^(div(fmpz(p-1),i)) for i in Fact])
					#WARNING note we still have too many tests here (to avoid inner loop)
					#@debug verify_primitive_elem(y, K) ? true : @error "primitive element wrong"
            		return y
				end
			end
        end
    end
end

function verify_primitive_elem(elem,K::Nemo.GaloisField)
	!iszero(elem) || return false
    p = Int(length(K))
    for i = 1:p-2
        elem^i != 1 || return false
    end
    return true
end

##
# implementation of the Magma_Sieve
##
function sieve_params(p,eps::Float64,ratio::Float64)
	#TODO more analysis and optimization of Sieve Params

	# assymtotic bounds by Coppersmith, Odlyzko, and Schroeppel L[p,1/2,1/2]# L[n,\alpha ,c]=e^{(c+o(1))(\ln n)^{\alpha }(\ln \ln n)^{1-\alpha }}}   for c=1
	qlimit = exp((0.5* sqrt(log(p)*log(log(p)))))
	qlimit *= log(qlimit) # since aproximately    #primes
	climit = exp((0.5+eps)*sqrt(log(p)*log(log(p))))

	qlimit = ceil(max(qlimit,30))
	climit = ceil(max(climit,35))
	steps = (ceil(0.15 * sqrt(p)/qlimit),ceil(0.15 * sqrt(p)/(log(climit)*climit)))
	return Sparam(qlimit,climit,ratio,steps)
end


function Sieve(F::FField,sieveparams::Sparam)

    p = characteristic(F.K) #(p = Int(length(A.K)))
    H = floor(root(p,2))+1
    J = H^2 - p

    # factorbase up to qlimit
    fb_primes = [i for i in PrimesSet(1,sieveparams.qlimit)]
    log2 = log(2.0);
    logqs = [log(q)/log2 for q in fb_primes] #real logarithms for sieve   #WARNING in source code this is more exact (LongFloat ?)


    #FB[findfirst(isequal(lift(alpha))] FB[1] = lift(alpha), FB[]
    FB = vcat([lift(F.a)],deleteat!(fb_primes,findfirst(isequal(lift(F.a)),fb_primes))) # tauschen a[1] = a[2] , a[2] = [1]
    FBs = deepcopy(FactorBase(FB))
    l = length(FB)
    Indx = Dict(zip(FB,[i for i=1:l])) #Index in a dictionary
    A = sparse_matrix(ZZ)

    for c1 = 1:sieveparams.climit
        nrows(A)/length(FB) < sieveparams.ratio || break
        sieve = zeros(sieveparams.climit)
        den = H + c1;                # denominator of relation
        num = -(J + c1*H)            # numerator
        for i=1:length(fb_primes)
            q = fb_primes[i]
            qpow = q
            nextqpow = qpow   #WARNING inserted, because of some error with nextqpow
            logq = logqs[i]
            while qpow < sieveparams.qlimit      # qlimit-smooth
                den % qpow != 0 || break
                c2 = num * invmod(den, fmpz(qpow))  % qpow
                (c2 != 0) || (c2 = qpow)
                nextqpow = qpow*q    #increase q_power
                while c2 < c1   #c2>=c1 to remove redundant relations of (c1,c2) tuples, just increase c2
                    c2+=qpow
                end
                while c2 <= length(sieve)
                    sieve[Int(c2)] += logq
                    if nextqpow > sieveparams.qlimit
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
    @debug !check_relation_mat(F.K,A,FB) ? (@error "Relation Matrix wrong") : nothing
	if nrows(A)/length(FB) < sieveparams.ratio
		#TODO global counter here for optimization of sieveparams
		sieveparams.qlimit += sieveparams.inc[1]
		sieveparams.climit += sieveparams.inc[2]
		return Sieve(F,sieveparams)
	end
    return A,FBs,FB
end

function check_relation_mat(K,A,FB)
    p = Int(length(K))
    for k = 1:nrows(A)
        prod = 1
        for i = 1:length(FB)
            ex = getindex(A,k,i)#getindex(A::SMat, i::Int, j::Int)
            if ex < 0
                prod *= invmod(FB[i],fmpz(p))^abs(ex) % p
            else
                prod  *= FB[i]^ex % p
            end
        end
        if !(prod % p  == 1)  return false end
    end
    return true
end


#Hecke.subset
#flog,clog,iroot

#badge smoothneth test -> glatttest (factorisieren)
##Examples:

B = FField(GF(1000000007),primitive_elem(GF(1000000007),true))
A,Q,C = Sieve(B, sieve_params(1000000007,0.02,1.1))

B = FField(GF(103),primitive_elem(GF(103),true))
A,Q,C = Sieve(B, sieve_params(103,0.02,1.1))

primitive_elem(GF(103),true)
