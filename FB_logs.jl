##########################################################################################################################################
#
Fachpraktikum = UInt32(0x00000002001) #v0x00.001.001
@info "Version-Fachpraktikum", Fachpraktikum
#
##########################################################################################################################################
using Hecke,Nemo,Revise,Profile,Markdown,TimerOutputs,Random
const to = TimerOutput()#show(to)
include("Magma_sieve.jl"),include("wiedemann_var.jl"),include("preprocessing.jl")
revise()
ENV["JULIA_DEBUG"] = "all" # enable debugging = "all" , disable:  = ""
##########################################################################################################################################
@doc Markdown.doc"""
    FB_logs(F::FField,storage=false) -> Tuple{Dict{fmpz, fmpz}, Vector{fmpz_mod}, FactorBase{fmpz}}
Compute a  `Factorbase` and a Dict of its `discrete logarithms` using a Indexcalculus algorithm.
"""
function FB_logs_new(F,prepro=false,prepro2=false,SP=sieve_params(p,0.02,1.1)::Sparam) #TODO this function for an FField
    ##########################################################################################################################################
    @timeit to "FB_logs_total" begin
    #for F FField find FB,FB_logs,FB_array
    p = characteristic(F.K)
    if p > typemax(Int)
        T = fmpz
    else
        T = Int64
        p = Int64(p)::T
    end 
    modulus_ = div((p-1),2)
    two = T(2)
    RR = ResidueRing(ZZ,modulus_)
    c,u,v = gcdx(two,modulus_)
    c == 1 || (@error "FB_LOGS: 2 ,(p-1)/2 not coprime")
    ##########################################################################################################################################
    #Sieve relations:
    @timeit to "Sieve:" RELMat,_,FB_x,l = Sieve(F,SP)
    ##########################################################################################################################################
    #Preprocessing
    @timeit to "PrePro:" begin
        if prepro 
            n,m = size(RELMat)
            TA = transpose(A)
            @debug @info "FB_LOGS: size of RELMAT before PrePro§; size(A) = $n x $m"
            RELMat,TA = sp_preprocessing_1(RELMat,TA,l)
            (n,m) = size(RELMat)
            @debug @info "FB_LOGS: size of RELMAT after PrePro1§; size(A) = $n x $m"
        end
        if prepro2
            n,m = size(RELMat)
            RELMat,_ = sp_preprocessing_2(RELMat,TA,l)
            (n,m) = size(RELMat)
            @debug @info "FB_LOGS: size of RELMAT after PrePro2§; size(A) = $n x $m"
        end
    end
    n,m = size(RELMat)
    ##########################################################################################################################################
    # get kernel
    RELMat = change_base_ring(RR,RELMat)
    cnt = 0
    @label retour
    @timeit to "Wiedemann" kern = wiedemann_var(RELMat)#,storage)
    cnt+=1
    cnt < 5 || return Dict{fmpz, fmpz}([]),Vector{fmpz_mod}([]),FactorBase(fmpz[])
    #TODO exeption if too many loops... / probably inf running time if kernel trivial
    @debug iszero(kern) ? (@info "FB_LOGS: trivial, found trivial kernel") : (@info "FB_LOGS: succeded wiedemann")
    !iszero(kern) || @goto retour
    kern = inv(kern[1]).*kern #norm kernelvec
    ##########################################################################################################################################
    # recon FB_logs mod p  mod p (note this works here if (p-1)/2 prime) Only 2 checks necesarry.
    @timeit to "Logdict" begin
    Q,L = Array{fmpz}([]),Array{fmpz}([])
    for i in 1:l
        temp = lift(kern[i])*two*u
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
    end
    ##########################################################################################################################################
    Logdict = Dict(zip(Q,L))
    @debug begin 
        check_logdict(F,Logdict,Q) ? (@info "FB_LOGS: Log_dict correct") : (@error "FB_LOGS: Log_dict incorrect")
        length(Logdict) == l ? (@info "FB_LOGS: all FB logs found") : (@warn "FB_LOGS: at least $(length(Logdict)-l) logarithms not found")
    end 
    ##########################################################################################################################################
    # check all Logs ans return 
    check_logdict(F,Logdict,Q) ? (@info "FB_LOGS: Log_dict correct") : (@error "FB_LOGS: Log_dict incorrect")
    length(Logdict) == l ? (@info "FB_LOGS: all FB logs found") : (@warn "FB_LOGS: at least $(length(Logdict)-l) logarithms not found") 
    return Logdict,kern,FactorBase(Q)
    end
end
##########################################################################################################################################
# Aux functions
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
    Indiv_Log(F,FB,FB_logs,h) -> Tuple{Dict{fmpz, fmpz}, Vector{fmpz_mod}, FactorBase{fmpz}}
Compute the discrete logarithm $log_F.a(h)$ i.e. compute an `x` s.t. `F.a^x = h` given a Factorbase FB with corresponding logarithms in FB_logs using a #TODO sieve.
"""
function Indiv_Log(F,FB,FB_logs,h) #TODO work on here.
    #return log_a(h) i.e x s.t a^x = h
    p = length(F.K)
    g = F.a
    randomexp = fmpz(rand(1:p-1))
    while !issmooth(FB,fmpz(lift(h*g^randomexp)))    # TODO NF-sieve or Q-sieve to find l
        randomexp = fmpz(rand(1:p-1))
    end  
    factorization = Hecke.factor(FB,lift(h*(F.a)^randomexp))

    log_h = -randomexp + sum([exp*FB_logs[prime] for (prime,exp) in factorization])
    @debug (F.a^lift(log_h) == h) ? nothing : (@error "calculation of log_h failed")
    return log_h
end 

function check_logdict_after(F,D)
    for (q,d) in D
        F.a^D[q] == q || return false 
    end 
    return true 
end
##########################################################################################################################################
# Testfields
ENV["JULIA_DEBUG"] = "" 
p = magma_p = fmpz(100000000000000000763)
p = cryptoprime(5)
TESTFIELD = BigFField(GF(p),primitive_elem(GF(p),true))

@time A,FB,FBx,l = Sieve(TESTFIELD,sieve_params(p,0.02,1.1))
N = div(p-1,2)
RR = ResidueRing(ZZ,N)
A = change_base_ring(RR,A)
reset_timer!(to)
@time FB_logs, _,FB = FB_logs_new(TESTFIELD,true,true)
show(to)
@time Indiv_Log(TESTFIELD,FB,FB_logs,-121231123)

check_logdict_after(TESTFIELD,D[1])
#=
check_logdict_after(TESTFIELD,Q[1])

ENV["JULIA_DEBUG"] = ""
p = cryptoprime(10)
TESTFIELD = BigFField(GF(p),primitive_elem(GF(p),true))

@profile (@debug @info "ok")
@time FB_logs(TESTFIELD,false)

@profile Sieve(TESTFIELD,SP)

Profile.clear() 
@profile Sieve(TESTFIELD,SP)
Profile.print(format= :flat, sortedby = :count, C = !true)
#Profile.print(sortedby = :count, C = !true)

#include("C:\\Users\\immos\\Documents\\Git\\Fachpraktikum\\IndexCalculus\\FB_logs.jl")

RR = ResidueRing(ZZ,fmpz(2))
a = RR(2)

typeof(a.data)
=#

#fmpz fast SMAt mul! 

##########################################################################################################################################sh
