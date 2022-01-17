#Pohlig Hellman Algorithm
include("wiedemann.jl"),include("Magma_sieve.jl")

function pohlig_hellman(g,h,F)
    #assume n = (p-1)/2 * 2   only factorization  
    #Input. A cyclic group G of order n with generator g, an element h G, and a prime factorization n=\prod _{i=1}^{r}p_{i}^{e_{i}}}
    #assume n = (p-1)/2 * 2   only factorization   
    #TODO later: n of unknown factorization
    #Output. The unique x s.t.  g^x=h   i.e x = disclog_g(h)
    p = length(F.K)
    L = [2,div(p-1,2)]
    Sol = [] 
    @label retour

    for i in L 
        g_i,h_i = g^div(p,i),h^div(p,i) 
        @debug isprime(lift(g_i)) ? nothing : (@warn "g_i not prime in Z") # NOTE we wont have this wl we will have g_i = g^2 or other powers, since g prime already
        F_i = FField(GF(i),g_i)
         #compute x_i such that g_i ^ x_i = h_i
        if i != 2 
            SP = sieve_params(p,0.02,1.1)
            RELMat,FB,FB_x = Sieve(F,SP)
            v = wiedemann(RELMat,i)
            lambda = invmod(v[1],i-1)
            v = lambda.*v
            FB_x2,new_logs = new_base_with_logs(F_i,FB_x,v,length(FB))
            #TODO maybe insert some H+ci to compensate wrong primes in FB
            randomexp = fmpz(rand(1:i-1))
            while !issmooth(FB_x2,a^randomexp*h_i)    # TODO NF-sieve or Q-sieve to find l
                randomexp = fmpz(rand(1:i-1))
            end  
            dict_factors2 = Hecke.factor(FBs,h_i*a^randomexp)
            log_hi = -randomexp + sum([exp*new_logs(prime) for (prime,exp) in dict_factors2])
            @debug (g_i^log_hi == h_i) ? nothing : (@error "subcalculation of log_hi failed")
        else
            #bruteforce
            #TODO Babystep/giantstep
        end 
        x_i = 0 #TODO
        if false #catched divisor of p-1 TODO 
            #update L 
            @goto retour 
        else 
        end
        push!(Sol,x_i)
    end 
    #use the CRT
end 

function new_base_with_logs(F::FField,FB_x,v,l)
    #input Field and FB + wiedemann kernel vector
    #return new FB and Dict with q_i and logq_i 
    #mult out a^log_a(v[i]) == v ?
    mask = [F.a for i = 1:length(FB_x)].^v .== FB_x
    @debug isszero(sum(mask)) ? (@error "all FB_x logs wrong") : nothing
    newFB_x = deepcopy(FB_x[mask]) #TODO inplace operation
    v_new = deepcopy(v[mask])
    logdict = Dict(zip(newFB_x[1:(l-sum(mask[1:l]))],v_new[l-sum(mask[1:l])]))
    return FactorBase(newFB_x[1:(l-sum(mask[1:l]))]),logdict
end


function crs_fmpz(N::Array{fmpz},P::Array{fmpz},positive=true)
    product = fmpz(1)
    for prime in P
        product = product*prime
    end
    R = fmpz(0)
    for i in 1:length(P)
        npi = divexact(product,P[i])  #wich is an Int because p | primeproduct
        mpi = invmod(npi,P[i])
        R += (N[i]*mpi*npi)
        R = rem(R,product)
    end
    positive ? R = R%product : return R%product
    while R<0 R+=product end
    return R
end


B = FField(GF(200087),primitive_elem(GF(200087),true))
A,Q,C = Sieve(B, sieve_params(200087,0.02,1.1))
h = B.a^(123)
t,z,u = pohlig_hellman(B.a,h,B)


(200087-1)/2

B = FField(GF(100043),GF(100043)(25))
Sieve(B, sieve_params(100043,0.02,1.1))