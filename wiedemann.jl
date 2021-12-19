using Hecke,Revise
revise()
function wiedemann(A,N,storage=false) # A in Z/NZ ^ n*m

	RR = ResidueRing(ZZ,N)
	#TODO reduce mod N = p-1
	#for now assume N prime (=> Z/NZ field)
	(n,m) = size(A);
	A = change_base_ring(RR, A)
	@debug (rank(Matrix(A)) == (m-1)) ? nothing : (@warn("WIEDEMANN: rank A small"),println(rank(Matrix(A))," != m-1 = ",m-1))
	TA = transpose(A)

	r = [RR(i) for i in rand(Int8,m)] # later generate random vector over ZZ / sampler ?
	c = [RR(i) for i in rand(Int8,m)]
	randlin = transpose([RR(i) for i in rand(Int8,m)])

	y = mul(TA,mul(A,r))
	# solve A^tAx = y2 => x -y in kernel(A^tA) to avoid finding zero vec

	#Wiedemann
	#TODO store all to save horner sheme
	#TODO store some and use a little horner sheme
	if storage 
		#store complete sequence
		M = zeros(RR,m,2*n)					 #preallocation in store
		k = @view M[:,1]
		k = c
		for i = 2:2*n
			M_last = @view M[:,i-1]			 
			M_i = @view M[:,i]				 #reduce allocations here
			M_i .= mul(TA,(mul(A,M_last)))   #generate sequence
		end
		done,f = Hecke_berlekamp_massey(M[rand(2:m-1),:])
	else
		seq = zeros(RR,2*n)
		seq0 = @view seq[1]
		seq0 .= randlin*c 					#preallocate here
		for i = 2:2*n
			seqi = @view seq[i]
			c = mul(TA,(mul(A,c))) # generate sequence
			seqi .= randlin*c
		end
		done,f = Hecke_berlekamp_massey(seq)
	end 

	
	@debug begin
		degr = degree(f)
		@info "WIEDEMANN: deg f = $degr where size(A^t*A) = $m"
		typeof(f) != fmpz || (@warn "ERLEKAMP_MASSEY: f may be constant polynom")
		done || (@warn "ERLEKAMP_MASSEY: modulus N is not prime, TODO: still catch some gcds")
		iszero(f(transpose(Matrix(A))*Matrix(A))) ? (@info "BERLEKAMP_MASSEY: valid return") : (@error "BERLEKAMP_MASSEY: unexpected return")
	end

	delta =0
	while iszero(evaluate(f,0)) #TODO collect coeffs:
		delta+=1
		f = divexact(f,gen(parent(f)))
	end

	@debug delta<2 || (@warn "WIEDEMANN: first nonzero coeff of f is a_$delta")
	constpartoff = evaluate(f,0)
	a = -inv(constpartoff)
	reducedf = divexact(f-constpartoff,gen(parent(f)))
	if storage 
		#TODO 
		coeff_vec = collect(coefficients(reducedf))
		v = mult(a,M[:,1:length(coeff_vec)]*coeff_vec)
	else 
		v = mult(a,horner_evaluate(reducedf,TA,A,y))
	end 
	@debug begin 
		y == (mul(TA,(mul(A,v)))) ? (@info "WIEDEMANN: Ax = y",true) : (@error "WIEDEMANN Ax = y",false)
		iszero(mul(A,v-r)) ? (@info "WIEDEMANN: A(v-r) = 0",true) : (@error "WIEDEMANN: A(v-r) = 0",false) 
	end 
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
	@debug f(transpose(Matrix(A))*Matrix(A))*c == s ? (@info "HORNER: f(A^t*A)c = s",true) : (@error  "HORNER: f(A^t*A)c = s",false)
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
	 #rems = gcd(lift(leading_coefficient(g)),M)
	 #if rems != 1
	 #	 return (false,rems)
	 #end
     N = R_s(inv(leading_coefficient(g))); g1 = g*N
     v0 = Ry(); v1 = Ry(1)
     while lg <= 2*degree(g1)
       q,r = divrem(f,g1)
       if r==0
         N = R_s(1)
       else
		 #rems = gcd(lift(leading_coefficient(r)),M)
  		 #if rems != 1
  		 #	 return (false,rems)
  		 #end
         N = R_s(inv(leading_coefficient(r)))
       end
       v = (v0-q*v1)*N
       v0 = v1; v1 = v; f = g1; g1= r*N
     end
     return true, divexact(v1, leading_coefficient(v1))
end
#TODO using .*
function mult(b, V)
  W = deepcopy(V)
  for i=1:length(W)
    W[i]*=b
  end
  return W
end

#TODO log degree f, size of m 
# change horner & save sequence in Matrix.
#TODO check savesequence implementation (not yet any running test)
