using Hecke,Nemo,Revise,Profile,Markdown
include("Magma_sieve.jl"),include("wiedemann.jl")
revise()
ENV["JULIA_DEBUG"] = "" # enable debugging = "all" , disable:  = ""

function check_logdict(F,D,Q)
    for q in Q 
        F.a^D[q] == q || return false 
    end 
    return true 
end 

function cryptoprime(N)
    #return a Prime p with N digits. s.t (p-1)/2 is prime
    p = rand(fmpz(10)^(N-1):fmpz(10)^N)
    while true
        p = next_prime(p+1)
        !isprime(div(p-1,2)) || return p
    end 
end 

@doc Markdown.doc"""
    FB_logs(F::FField,storage=false) -> Tuple{Dict{fmpz, fmpz}, Vector{fmpz_mod}, FactorBase{fmpz}}
Compute a  `Factorbase` and a Dict of its `discrete logarithms` using a Indexcalculus algorithm.
"""
function FB_logs(F::FField,storage=false)
    #for F FField find FB,FB_logs,FB_array
    p = length(F.K)
    modulus_ = fmpz((p-1)/2)
    two = fmpz(2)
    RR = ResidueRing(ZZ,modulus_)
    c,u,v = gcdx(two,modulus_)
    @debug c == 1 || (@error "FB_LOGS: 2 ,(p-1)/2 not coprime")

    #Sieve relations:
    SP = sieve_params(p,0.02,1.15)
    RELMat,FB,FB_x,l = Sieve(F,SP)

    #get kernel mod p-1 / 2 
    RELMat = change_base_ring(RR,RELMat)
    #n,m = size(RELMat);
    #dim,kern = kernel(Matrix(RELMat)) #TODO wiedemann CRT here
    #@debug dim == 1 || (@warn "dim(ker(A)mod(p-1)/2) > 1")
    @label retour
    kern = wiedemann(RELMat,modulus_,storage)
    #TODO exeption if too many loops... / probably inf running time if kernel trivial
    @debug iszero(kern) ? (@info "FB_LOGS: trivial found trivial kernel") : (@info "FB_LOGS: succeded wiedemann")
    !iszero(kern) || @goto retour
    #return inv(v[1]).*v
    #reconstruct mod p (note this works here if (p-1)/2 prime) Only 2 checks necesarry.
    #return kern,u,v
    kern = inv(kern[1]).*kern

    Q,L = Array{fmpz}([]),Array{fmpz}([])
    for i in 1:l
        temp = lift(kern[i])*fmpz(2)*u
        test1 = temp%(p-1)
        test2 = (temp + v*modulus_)%(p-1)
        q_temp = FB_x[i]
        if F.a^test1 == q_temp
            push!(Q,q_temp)
            push!(L,fmpz(test1))
        elseif F.a^test2 == FB_x[i]
            push!(Q,q_temp)
            push!(L,fmpz(test2))
        end 
    end 
    #Indx = Dict(zip(Q,[i  for i=1:length(Q)]))
    Logdict = Dict(zip(Q,L))
    @debug begin 
        check_logdict(F,Logdict,Q) ? (@info "FB_LOGS: Log_dict correct") : (@error "FB_LOGS: Log_dict incorrect")
        length(Logdict) == l ? (@info "FB_LOGS: all FB logs found") : (@warn "FB_LOGS: at least $(length(Logdict)-l) logarithms not found")
    end 
    check_logdict(F,Logdict,Q) ? (@info "FB_LOGS: Log_dict correct") : (@error "FB_LOGS: Log_dict incorrect")
    length(Logdict) == l ? (@info "FB_LOGS: all FB logs found") : (@warn "FB_LOGS: at least $(length(Logdict)-l) logarithms not found") 
    return Logdict,kern,FactorBase(Q)
end
@doc Markdown.doc"""
    Indiv_Log(F,FB,FB_logs,h) -> Tuple{Dict{fmpz, fmpz}, Vector{fmpz_mod}, FactorBase{fmpz}}
Compute the discrete logarithm $log_F.a(h)$ i.e. compute an `x` s.t. `F.a^x = h` given a Factorbase FB with corresponding logarithms in FB_logs using a #TODO sieve.
"""
function Indiv_Log(F,FB,FB_logs,h)
    #return log_a(h) i.e x s.t a^x = h
    p = length(F.K)
    g = F.a
    randomexp = fmpz(rand(1:p-1))
    while !issmooth(FB,lift(h*g^randomexp))    # TODO NF-sieve or Q-sieve to find l
        randomexp = fmpz(rand(1:p-1))
    end  
    factorization = Hecke.factor(FB,lift(h*(F.a)^randomexp))

    log_h = -randomexp + sum([exp*FB_logs[prime] for (prime,exp) in factorization])
    @debug (F.a^lift(log_h) == h) ? nothing : (@error "calculation of log_h failed")
    return log_h
end 

#=
TODO list so far:

>> Implement good logger to quick_overview/present performance/correctness
>> debug/improve with @code_warntype 
>> optimize with  @profile
>> Implement block wiedemann for faster performanve 
>> implement sieve for faster finding l for individual logs

#REMARKS
 do we really need to save FB_x ? maybe only index 
 use original FB to recover logs.
 length(FB.base) for length

#Further notes:
>> flog,clog,iroot
>> badge smoothneth test -> glatttest (factorisieren)

=#

p = cryptoprime(20)
TESTFIELD = FField(GF(p),primitive_elem(GF(p),true))
@time FB_logs(TESTFIELD,false)
@time FB_logs(TESTFIELD,true)
#~5 minutes. for p = cryptoprime(20) 
#~34,5 minutes. for p = 3088833293915623767369443
Profile.print(format=:flat)

function check_logdict_after(F,D)
    for (q,d) in D
        F.a^D[q] == q || return false 
    end 
    return true 
end
check_logdict_after(TESTFIELD,a)
Profile.print(:flat)



@time disc_log_bs_gs(primitive_elem(GF(p),true), GF(p)(449), fmpz(64768854538896063587))
#start ca. 01:00
#abbruch ca. 01:15 
#Preallovation of Memory... 
