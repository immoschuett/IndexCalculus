include("wiedemann_split.jl")
function FBlogs(F::FField)
    #for F FField find FB,FB_logs,FB_array
    p = length(F.K)
    g = F.a
    SP = sieve_params(p,0.02,1.3)
    RELMat,FB,FB_x,l = Sieve(F,SP)
    RELMat = change_base_ring(ResidueRing(ZZ,fmpz((p-1))),RELMat)
    n,m = size(RELMat)
    println(n,m)
    dim,kern = kernel(Matrix(RELMat))
    o = kern[:,1]
    mask = [F.a^lift(o[i]) for i = 1:length(o)] .== FB_x
    l = min(l,(sum(mask)))
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
Dict_logs
Indiv_Log(B,FB,FB_logs,B.a^1234)