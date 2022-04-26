##########################################################################################################################################
revise()
#Types
@doc Markdown.doc"""
    abstract type DLP
special type for `Discrete Logarithm Problem`
"""
abstract type DLP end
@doc Markdown.doc"""
    FField<:DLP
A FFIELD $F$ with Fields
```jldoctest
F.K::Nemo.GaloisField `finite Field`
F.a::gfp_elem `primitive element` 
```
"""
mutable struct FField<:DLP
	K::Nemo.GaloisField
	a::gfp_elem #primitive element
end
@doc Markdown.doc"""
    BigFField<:DLP
A FFIELD $F$ with Fields
```jldoctest
F.K::Nemo.GaloisFmpzField `finite Field`
F.a::gfp_fmpz_elem `primitive element` 
```
"""
mutable struct BigFField<:DLP
	K::Nemo.GaloisFmpzField
	a::gfp_fmpz_elem #primitive element
end
@doc Markdown.doc"""
    Sparam<:DLP
Sieveparameters $S$ with Fields
```jldoctest
S.qlimit::Int6
S.climit::Int64
S.ratio::Float64
S.inc::Tuple{Int64, Int64} `step increasing parameters` 
```
generated with `sieve_params(p,eps::Float64,ratio::Float64)`
"""
mutable struct Sparam<:DLP
	qlimit::Int64
	climit::Int64
	ratio::Float64
	inc::Tuple{Int64, Int64} #for increasing of Sparam
end
##########################################################################################################################################
#Types 2
function sieve_params(p,eps::Float64,ratio::Float64)
	#TODO more analysis and optimization of Sieve Params

	# assymtotic bounds by Coppersmith, Odlyzko, and Schroeppel L[p,1/2,1/2]# L[n,\alpha ,c]=e^{(c+o(1))(\ln n)^{\alpha }(\ln \ln n)^{1-\alpha }}}   for c=1
	qlimit = exp((0.5* sqrt(log(p)*log(log(p)))))
	qlimit *= log(qlimit) # since aproximately    #primes
	climit = exp((0.5+eps)*sqrt(log(p)*log(log(p))))

	qlimit = ceil(0.5*max(qlimit,30))
	climit = ceil(max(climit,35))
	steps = (0,100)
	return Sparam(qlimit,climit,ratio,steps)
end
function sieve_params2(q_ratio,p,eps::Float64,ratio::Float64)
	#TODO more analysis and optimization of Sieve Params

	# assymtotic bounds by Coppersmith, Odlyzko, and Schroeppel L[p,1/2,1/2]# L[n,\alpha ,c]=e^{(c+o(1))(\ln n)^{\alpha }(\ln \ln n)^{1-\alpha }}}   for c=1
	qlimit = exp((0.5* sqrt(log(p)*log(log(p)))))
	qlimit *= log(qlimit) # since aproximately    #primes
	climit = exp((0.5+eps)*sqrt(log(p)*log(log(p))))

	qlimit = ceil(q_ratio*max(qlimit,30))
	climit = ceil(max(climit,35))
	steps = (0,300)
	return Sparam(qlimit,climit,ratio,steps)
end

function set_sieve_params(q::Int64,c::Int64,r::Float64,i::Tuple{Int64, Int64})
    #later some stuff here
    return Sparam(q,c,r,i)
end 
##########################################################################################################################################
#Sieve
@doc Markdown.doc"""
    Sieve(F::FField,sieveparams::Sparam) -> A,FB? TODO work in progress
Sieve a relation matrix $A$ s.t `ker(A)` includes all `discrete logarithms` of $FB$-elements"""
function Sieve(F::BigFField,sieveparams::Sparam)
    @debug @info "SIEVE: starting Sieve"
    p = characteristic(F.K) #(p = Int(length(A.K)))
    H = floor(root(p,2))+1
    J = H^2 - p

    # factorbase up to qlimit
    fb_primes = [i for i in PrimesSet(1,sieveparams.qlimit)]


    #FB[findfirst(isequal(lift(alpha))] FB[1] = lift(alpha), FB[]
    indx = findfirst(isequal(lift(F.a)),fb_primes)
    FB = vcat([lift(F.a)],deleteat!(fb_primes,indx)) # tauschen a[1] = a[2] , a[2] = [1] 
    # use shift! / unshift! here...
    log2 = log(2.0);
    logqs = [log(Int(q))/log2 for q in FB] #real logarithms for sieve 
    FBs = deepcopy(FactorBase(FB))
    l2 = length(FB)
    l = deepcopy(l2)
    Indx = Dict(zip(FB,[i for i=1:l])) #Index in a dictionary
    A = sparse_matrix(ZZ)
    #IDEA: dont add logs. add INT counter, then add cnt * log in the end. ??
    ##########################################################################################################################################
    # Sieve for ci
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
        ##########################################################################################################################################
        #include relations / check sieve for full factorizations.
        rel = den * (H + 1)
        relinc = H + c1       # add to relation to get next relation
        idx = 0
        for c2 in 1:length(sieve)
            n = rel % p
            if abs(sieve[c2] - (nbits(n)-1)) < 1 #TODO nbits(n)+1
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
                    J1 = Vector{Int64}([])
                    V = Vector{Int64}([])
                    for (prime,power) in dict_factors
                        if !(power == 0)
                            push!(J1,Indx[prime])
                            push!(V,Int(power))
                        end
                    end
                    if c1 == c2
                         push!(J1,Indx[H+c1])
                         push!(V,-2)
                    else
                         push!(J1,Indx[H+c1])
                         push!(J1,Indx[H+c2])
                         push!(V,-1)
                         push!(V,-1)
                    end
                    push!(A,sparse_row(ZZ, J1, V))
                end
            end
            rel += relinc
        end
    end
    ##########################################################################################################################################
    #increase Sieve 
    # IDEA: reuse already given Sieve vec. and only sieve for new c_i (otherwise exp waste of time)
    # still TODO
    if nrows(A)/length(FB) < sieveparams.ratio
        @debug @info "SIEVE: increase sieveparams"
		#TODO global counter here for optimization of sieveparams
		sieveparams.qlimit += sieveparams.inc[1]
		sieveparams.climit += sieveparams.inc[2]
		return Sieve(F,sieveparams)
	end
    @debug check_relation_mat(F.K,A,FB) ? (@info "SIEVE: all Relations in Matrix correct") : (@error "SIEVE: Relation Matrix wrong")
    @info "SIEVE_params:", sieveparams
    return A,FBs,FB,l
end
##########################################################################################################################################
## Aux functions
function check_relation_mat(K,A,FB)
    p = BigInt(length(K))
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
function primitive_elem(K::FinField,first::Bool) #TODO implement compatible for Abstract field or Nemo.GaloisfmpzField
    #returns a (the first) generator alpha of K° s.t. lift(alpha) is prime in ZZ
    p = length(K)
    Fact = prime_divisors(fmpz(p-1))
    while true # alpha exists
        for y in K
			bool = true
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
function verify_primitive_elem(elem,K)
	!iszero(elem) || return false
    p = Int(length(K))
    for i = 1:p-2
        elem^i != 1 || return false
    end
    return true
end
##########################################################################################################################################
# Comments
#TODO instead invcrease steps. 
#use @goto increase c limit 
# dont sieve again all this stuff.
#dynamic sieving. 
#write expanding iterator 