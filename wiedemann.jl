using Hecke,Revise
revise()

function wiedemann(A,N=103) # A in Z/NZ ^ n*m
	RR = ResidueRing(ZZ,N)
	#TODO reduce mod N = p-1
	# assume N prime ( Z/NZ field)
	# c in R^m  random ?
	(n,m) = size(A)
	A = change_base_ring(RR, A)
	TA = transpose(A) #later generate random sparse matrix over ZZ
	r_random = rand(Int8,m)
	c_random = rand(Int8,m)
	r = [RR(i) for i in r_random] # later generate random vector over ZZ / sampler ?
	c = [RR(i) for i in c_random]
	origc = deepcopy(c)
	y = mul(TA,(mul(A,r)))
	# solve A^tAx = y  => x -y in kernel(A^tA)
	#Wiedemann
	seq = fmpz_mat((zeros(Int,10,2*n)))
	#seq = convert(Array{fmpz},(zeros(Int,10,2*n)))#store 10 entrys of the vector sequence
	#TODO store all to save horner sheme
	#TODO store some and use a little horner sheme
	seq = change_base_ring(RR, seq)
	#seq = convert(Array{nmod},seq)#store 10 entrys of the vector sequence
	for i = 1:2*n-1
		c = mul(TA,(mul(A,c))) # generate sequence
		seq[:,i+1] = c[1:10]
	end
	#return (seq[2,:])
	done,f = Hecke_berlekamp_massey(seq[2,:])
	@debug @assert done "modulus N is not prime, TODO: still catch some gcds"
	v = horner_evaluate(f,TA,A,origc)
	@debug @assert iszero(mul(TA,(mul(A,v-y)))) "not a kernel vec"
	#return mul(TA,(mul(A,v)))
	return v-y

end

# we use horner sheme to use the sparsity of A
function horner_evaluate(f::nmod_poly,TA,A,c)  where T
	C = collect(coefficients(f))
	n = length(C)
	s = mult(C[end],mul((TA),mul(A,c)))
	for i = n-1:-1:1
		#println(s)
		#TA * A * s + fi * c
		s = mul((TA),mul(A,s)) + mult(C[i],c)
	end
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

function mult(b::nmod, V::Vector{nmod})
  for i=1:length(V)
    V[i]*=b
  end
  return V
end

##
wiedemann(A)
