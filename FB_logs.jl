using Hecke,Nemo
include("Magma_sieve.jl"),include("wiedemann.jl")
function FBlogs(F::FField)
    #for F FField find FB,FB_logs,FB_array
    p = length(F.K)
    g = F.a
    SP = sieve_params(p,0.02,1.3)
    RELMat,FB,FB_x,l = Sieve(F,SP)
    RELMat = change_base_ring(ResidueRing(ZZ,fmpz((p-1))),RELMat)
    n,m = size(RELMat)
    println(n,m)
    dim,kern = kernel(Matrix(RELMat)) #TODO wiedemann CRT here
    o = kern[:,1]
    mask = [F.a^lift(o[i]) for i = 1:length(o)] .== FB_x
    l = min(l,(sum(mask))) #TODO other notation here
    FB = FB_x[mask][1:l]
    o = Array(o)[mask][1:l]
    return Dict(zip(FB,o)),kern,FactorBase(FB)
end

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
    @debug (F.a^log_h == h) ? nothing : (@error "calculation of log_h failed")
    return log_h
end 


B = FField(GF(200087),primitive_elem(GF(200087),true))
g = primitive_elem(GF(200087),true)
FB_logs,kern,FB = FBlogs(B)
kern
FB_logs
Indiv_Log(B,FB,FB_logs,B.a^1234)


function FBlogs_new(F::FField)
    #for F FField find FB,FB_logs,FB_array
    p = length(F.K)
    modulus_ = fmpz((p-1)/2)
    two = fmpz(2)
    RR = ResidueRing(ZZ,modulus_)
    c,u,v = gcdx(two,modulus_)
    @debug c == 1 || (@error "2 ,(p-1)/2 not coprime")

    #Sieve relations:
    SP = sieve_params(p,0.02,1.3)
    RELMat,FB,FB_x,l = Sieve(F,SP)

    #get kernel mod p-1 / 2 
    RELMat = change_base_ring(RR,RELMat)
    n,m = size(RELMat)
    #dim,kern = kernel(Matrix(RELMat)) #TODO wiedemann CRT here
    #@debug dim == 1 || (@warn "dim(ker(A)mod(p-1)/2) > 1")
    @label retour
    kern = wiedemann(RELMat,modulus_)

    @debug iszero(kern) ? (@info "trivial found trivial kernel") : (@info "succeded wiedemann")
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
    println(typeof(Q))
    Logdict = Dict(zip(Q,L))
    @debug check_logdict(F,Logdict,Q) ? (@info "Log_dict correct") : (@error "Log_dict incorrect")
    @debug length(Logdict) ==l ? (@info "all FB logs found") : (@warn "at least " print(length(Logdict)-l) " not found")
    return Logdict,kern,FactorBase(Q)
end

function check_logdict(F,D,Q)
    for q in Q 
        F.a^D[q] == q || return false 
    end 
    return true 
end 

dict,kern,FB = FBlogs_new(B)
typeof(kern[6])
inv(kern[6])
crs_fmpz([fmpz(25931),fmpz(0)],[fmpz(100043),fmpz(2)])

RR = ResidueRing(ZZ,200086)
RR(lift(kern[6]) * fmpz(2) * fmpz(-50021))
