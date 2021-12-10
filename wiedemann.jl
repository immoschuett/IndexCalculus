using Hecke,Revise
ENV["JULIA_DEBUG"] = "all" # enable debugging
revise()
include("Magma_sieve.jl")

function wiedemann(A,N) # A in Z/NZ ^ n*m
	RR = ResidueRing(ZZ,N)
	#TODO reduce mod N = p-1
	#for now assume N prime (=> Z/NZ field)

	(n,m) = size(A)
	A = change_base_ring(RR, A)
	@debug (rank(Matrix(A)) == (m-1)) ? nothing : (@warn("rank A small"),println(rank(Matrix(A))," != m-1 = ",m-1))
	TA = transpose(A) #later generate random sparse matrix over ZZ

	r = [RR(i) for i in rand(Int8,m)] # later generate random vector over ZZ / sampler ?
	c = [RR(i) for i in rand(Int8,m)]
	randlin = transpose([RR(i) for i in rand(Int8,m)])
	origc = deepcopy(c)
	y = mul(TA,mul(A,r))
	# solve A^tAx = y2 => x -y in kernel(A^tA) to avoid finding zero vec

	#Wiedemann
	#TODO store all to save horner sheme
	#TODO store some and use a little horner sheme
	seq = [randlin*c]
	for i = 1:2*n-1
		c = mul(TA,(mul(A,c))) # generate sequence
		push!(seq,randlin*c)
	end
	done,f = Hecke_berlekamp_massey(seq)
	@debug !iszero(f(transpose(Matrix(A))*Matrix(A))) ? (@warn "something wrong with berlekamp-massey return") : nothing
	@debug !done ? (@warn "modulus N is not prime, TODO: still catch some gcds") : nothing

	delta =0
	if typeof(f) == fmpz
		@error "something wrong with polynom f"
	end
	while iszero(evaluate(f,0))
		delta+=1
		f = divexact(f,gen(parent(f)))
	end
	@debug delta>1 ? (println(delta), @warn "something wrong with delta delta") : nothing
	constpartoff = evaluate(f,0)
	a = -inv(constpartoff)
	reducedf = divexact(f-constpartoff,gen(parent(f)))
	compared = a*reducedf(transpose(Matrix(A))*Matrix(A))*y
	v = mult(a,horner_evaluate(reducedf,TA,A,y))
	@debug !(y == transpose(Matrix(A))*(Matrix(A)*compared)) ? (@error "compared wrong") : nothing
	@debug !(y == (mul(TA,(mul(A,v))))) ? (@error "not Ax = y") : nothing
	@debug !(iszero(mul(TA,(mul(A,v-r))))) ? (@error "not a kernelvec of TA A") : nothing
	@debug !(iszero(mul(A,v-r))) ? (@error "not a kernelvec of A") : nothing
	return v-r
end

# we use horner sheme to use the sparsity of A
function horner_evaluate(f,TA,A,c)
	#return f(A^t *A)*c
	C = collect(coefficients(f))
	n = length(C)
	s = mult(C[end],mul(TA,mul(A,c)))+mult(C[end-1],c)
	for i = n-2:-1:1   #WARNING a_0 in papers, but a_1 in julia
		#s = A^t * A * s + fi * c  inloop
		s = mul(TA,mul(A,s)) + mult(C[i],c)
	end
	@debug !(f(transpose(Matrix(A))*Matrix(A))*c == s) ? (@error "error in horner sheme implementation real value") : nothing
	return s
end

function Hecke_berlekamp_massey(L)#::Vector{fmpz})
	# from Hecke\U6dsX\src\Sparse\Matrix.jl
	 RR = parent(L[1])
	 M = modulus(RR)
	 Ry,x = PolynomialRing(RR, "x", cached = false) ## Ring over TZZ
     R_s = ZZ
     lg = length(L)
     L = [R_s(L[lg-i]) for i in 0:lg-1]
     #Y = gen(Ry)
     g = Ry(L)
     if iszero(g)
       return true, g
     end
     f = x^lg
	 rems = gcd(lift(leading_coefficient(g)),M)
	 if rems != 1
		 return (false,rems)
	 end
     N = R_s(inv(leading_coefficient(g))); g1 = g*N
     v0 = Ry(); v1 = Ry(1)
     while lg <= 2*degree(g1)
       q,r = divrem(f,g1)
       if r==0
         N = R_s(1)
       else
		 rems = gcd(lift(leading_coefficient(r)),M)
  		 if rems != 1
  			 return (false,rems)
  		 end
         N = R_s(inv(leading_coefficient(r)))
       end
       v = (v0-q*v1)*N
       v0 = v1; v1 = v; f = g1; g1= r*N
     end
     return true, divexact(v1, leading_coefficient(v1))
end

function mult(b, V)
  W = deepcopy(V)
  for i=1:length(W)
    W[i]*=b
  end
  return W
end