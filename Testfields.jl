##########################################################################################################################################
# Testfields
##########################################################################################################################################
using Hecke,Plots,LaTeXStrings
# TEST for Timing, Paramter adjustment and vice-versa.
function timetest(size)
    time = 0
    for i = 1:10
        p =  cryptoprime(size)
        while p > typemax(Int64)
            p =  cryptoprime(size)
        end 
        parameter = r_sieve_params(P,p)
        TESTFIELD = BigFField(GF(p),primitive_elem(GF(p),true))
        for i = 1:10
            reset_timer!(to)
            FB_logs_new(TESTFIELD,(true,true,5),parameter)
            time += TimerOutputs.time(to["FB_logs_total"])/(10*10*1000*1000)
        end 
    end
    return time # time is defailt in nanosec -> milliseconds
end  
function sizetest(size)
    n = 0
    n2 = 0
    prepro = (true,true,5)
    for i = 1:10
        p =  cryptoprime(size)
        while p > typemax(Int64)
            p =  cryptoprime(size)
        end 
        parameter = r_sieve_params(P,p)
        TESTFIELD = BigFField(GF(p),primitive_elem(GF(p),true))
        for i = 1:10
            A,_,_,l = Sieve(TESTFIELD,parameter)
            n += nrows(A)
            RR = ResidueRing(ZZ,p)
            A = change_base_ring(RR,A)
            TA = transpose(A)
            B,_ = sp_preprocessing(A,TA,l,prepro[3],prepro[2])
            n2 += nrows(B)
        end 
    end
    return (n/100,n2/100) # time is defailt in nanosec -> milliseconds
end 
function testrun1()
    p = fmpz(536485931783382683)
    TESTFIELD = BigFField(GF(p),primitive_elem(GF(p),true))
    upper_q = 7000
    lower_q = 1000
    steps = 100
    L = []
    c_last = 5
    for q in upper_q:-steps:lower_q
        params = Sieve(TESTFIELD,sieve_params_experimental(q,c_last,1.01)) #only works with return sieveparams in SIEVE
        c_last = params.climit
        push!(L,params)
    end 
    return L
end 
function testrun2(L)
    p = fmpz(536485931783382683)
    TESTFIELD = BigFField(GF(p),primitive_elem(GF(p),true))
    Timedata=[]
    for param in L
        reset_timer!(to)
        for i = 1:10
           FB_logs_new(TESTFIELD,(true,true,5),param)
        end 
        push!(Timedata,[param, TimerOutputs.time(to["FB_logs_total"])/10000000 , TimerOutputs.time(to["FB_logs_total"]["Sieve:"])/10000000 , TimerOutputs.time(to["FB_logs_total"]["Wiedemann:"])/10000000 , TimerOutputs.time(to["FB_logs_total"]["PrePro:"])/10000000 ])
    end 
    return Timedata
end 
#=
ENV["JULIA_DEBUG"] = "" 
#p = magma_p = fmpz(100000000000000000763)
p =  cryptoprime(25)
TESTFIELD = BigFField(GF(p),primitive_elem(GF(p),true))
t = sieve_params(p,0.001,0.02,1.01,(Float64,1.0))
t = r_sieve_params(P,p)
reset_timer!(to)
FB_logs(TESTFIELD,(true,true,5),t)
show(to)
Sieve(TESTFIELD,t)
TimerOutputs.time(to["FB_logs_total"])/(1000*1000)
A,_,_,_ = Sieve(TESTFIELD,t)
[sizetest(i) for i in 5:19]
[(i,log2(10^i)) for i = 4:18 ]
# wieviele Zahlen gibt es, die glatt sind ? 
L = testrun1()
LL = copy(L)
DATA = testrun2(L)
D = [i[k] for i in DATA, k in 2:5]
L0 = [i for i in 7000:-100:1000]
labels = ["Gesamtzeit" "Siebzeit" "Kernberechnungszeit" "Preprocessing"]
plot(L0,D,
     seriestype = :scatter,
     discrete_values = true, 
     framestyle = :box, 
     #ticks = 1000:100:4000,
     markers = :xcross,
     #color = 4, 
     #title = latexstring("Im  (\\log_{2}(x))" ), 
     label= labels,
     ylabel = "average time in "*latexstring("\\mu s"),
     xlabel = latexstring("B (= B_F \\log(B_F))")
     )


=#


##########################################################################################################################################
# Example:
#=
ENV["JULIA_DEBUG"] = "" 
p = magma_p = fmpz(100000000000000000763)
p = cryptoprime(15)
TESTFIELD = BigFField(GF(p),primitive_elem(GF(p),true))

FB_logs(TESTFIELD)
Sieve(TESTFIELD,r_sieve_params(P,p))
=#