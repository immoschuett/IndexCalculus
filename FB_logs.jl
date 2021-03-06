##########################################################################################################################################
#
Fachpraktikum = UInt32(0x00000002001) #v0x00.002.001
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
function FB_logs(F,prepro=(true,true,5)::Tuple{Bool, Bool, Int64},SP=sieve_params(p,0.0,0.1,1.1,(Float64,1.0))::Sparam) #TODO this function for an FField
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
        RELMat = change_base_ring(RR,RELMat)
        if prepro[1]
            n,m = size(RELMat)
            TA = transpose(RELMat)
            dens = density(RELMat)
            @info "FB_LOGS: size of RELMAT before PrePro??; size(A) = $n x $m with density $dens"
            RELMat,TA = sp_preprocessing(RELMat,TA,l,prepro[3],prepro[2])
            (n,m) = size(RELMat)
            dens = density(RELMat)
            @info "FB_LOGS: size of RELMAT after PrePro??; size(A) = $n x $m with density $dens"
        end
    end
    n,m = size(RELMat)
    ##########################################################################################################################################
    # get kernel
    cnt = 0
    @label retour
    @timeit to "Wiedemann:" kern = wiedemann_var(RELMat)#,storage)
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
        length(Logdict) == l ? (@info "FB_LOGS: all FB logs found") : (@warn "FB_LOGS:  $(l-length(Logdict)) logarithms not found")
    end 
    ##########################################################################################################################################
    # check all Logs ans return 
    check_logdict(F,Logdict,Q) ? (@info "FB_LOGS: Log_dict correct") : (@error "FB_LOGS: Log_dict incorrect")
    length(Logdict) == l ? (@info "FB_LOGS: all FB logs found") : (@warn "FB_LOGS: $(l-length(Logdict)) logarithms not found") 
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