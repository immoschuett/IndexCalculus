
#Input. A cyclic group G of order n with generator g, an element h G, and a prime factorization n=\prod _{i=1}^{r}p_{i}^{e_{i}}}
#assume n = (p-1)/2 * 2   only factorization   
#TODO later: n of unknown factorization with pohlig hellmann
#Output. The unique x s.t.  g^x=h   i.e x = disclog_g(h)

function disc_logs(F,h)
    p = length(F.K)
    @debug issprime((p-1)/2) ? nothing : (@warn "(p-1)/2 not prime")

    #precalculate disc log relations
    SP = sieve_params(p,0.02,1.1)
    RELMat,FB,FB_x = Sieve(F,SP)

    #find q_logs mod p-1/2
    modul = fmpz((p-1)/2)
    g_temp1 = F.a^2
    v = wiedemann(RELMat,modul)
    v = (invmod(lift(v[1]),modul)).*v #norm
    return FB
    FB_2,FB_x2,new_logs = new_base_with_logs(FField(GF(modul),g_temp1),FB_x,v,length(FB))
    
    #find log_q mod 2 
    g_temp2 = F.a^((p-1)/2)
    Solmod2 = g_temp == FB_2
    
    #reconstruct logs mod p-1 with CRT
    P = [fmpz(2),fpmz((p-1)/2)]
    i = 0
    C = []
    #C = push([],crs_fmpz([new_logs[q],Solmod2[indx_q] for q in FB_2],P)) TODO
    for q in FB_2
        i+=1
        C = push(C,crs_fmpz([new_logs[q],Solmod2[i]],P))
    end 
    @debug iszero((F.a).^C - FB_2) ? nothing : (@error "still some logs mod p-1 wrong") 

    #apply on special h 
    randomexp = fmpz(rand(1:i-1))
    while !issmooth(FB_x2,h*(F.a)^randomexp)    # TODO NF-sieve or Q-sieve to find l
        randomexp = fmpz(rand(1:i-1))
    end  
    dict_factors2 = Hecke.factor(FB_x2,h*(F.a)^randomexp)
    log_h = -randomexp + sum([exp*new_logs(prime) for (prime,exp) in dict_factors2])
    @debug (g^log_h == h) ? nothing : (@error "culation of log_h failed")
    return log_h
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
    return newFB_x[1:(l-sum(mask[1:l]))],FactorBase(newFB_x[1:(l-sum(mask[1:l]))]),logdict
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
RELMat,FB,FB_x = Sieve(B, sieve_params(200087,0.02,1.1))
h = B.a^(123)
t = disc_logs(B,h)

println(t)
#find q_logs mod p-1/2
modul = fmpz((200087-1)/2)
g_temp1 = F.a^2
v = wiedemann(RELMat,modul)