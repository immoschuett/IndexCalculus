#Pohlig Hellman Algorithm
include("wiedemann.jl"),include("Magma_sieve.jl")

function pohlig_hellman()
    #assume n = (p-1)/2 * 2   only factorization  
    #Input. A cyclic group G of order n with generator g, an element h G, and a prime factorization n=\prod _{i=1}^{r}p_{i}^{e_{i}}}
    #assume n = (p-1)/2 * 2   only factorization   
    #TODO later: n of unknown factorization
    #Output. The unique x s.t.  g^x=h   i.e x = disclog_g(h)
    L = [2,div(p-1,2)]
    Sol = [] 
    @label retour

    for i in L 
        [g_i,h_i] .= [g,h].^div(n,i)
        @debug isprime(lift(g_i)) ? nothing : (@error "g_i not prime in Z")
        F_i = FField(GF(i),g_i)
         #compute x_i such that g_i ^ x_i = h_i
        if i != 2 
            SP = sieve_params(i,0.02,1.1)
            RELMat,FB,FB_x = Sieve(F_i,SP)
            v = wiedemann(RELMat,i-1)
            lambda = invmod(v[1],i-1)
            v = lambda.*v
            FB_x2 = new_base_logs(F_i,FB_x,v,length(FB))
            #TODO maybe insert some H+ci to compensate wrong primes in FB
            while true 
                randomexp = rand(1:i-1)
                # search for l s.t a^l h_i is  Q_new smoth 
                # TODO NF-sieve or Q-sieve to find l
                if !issmooth(FBs,fmpz(h_i*a^randomexp)) break end 
                dict_factors2 = Hecke.factor(FBs,h_i*a^randomexp)
                # some sol log = sum of log * exponents push!(Sol,sum())

            end 
            #

        else 
            #bruteforce
            #TODO Babystep/giantstep
        end 

        x_i = 0 #TODO
        if false #catched divisor of p-1 TODO 
            #update L 
            @goto retour 
        end 
        push!(Sol,x_i)
    end 
    #use the CRT
end 

function new_base_logs(F::FField,FB_x,v,l)
    #mult out a^log_a(v[i]) == v ?
    mask = [F.a for i = 1:length(FB_x)].^v .== FB_x
    @debug isszero(sum(mask)) ? (@error "all FB_x logs wrong") : nothing
    newFB_x = FB_x[mask]
    return FactorBase(newFB_x[1:(l-sum(mask[1:l]))])
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


function subdisclog()

end
