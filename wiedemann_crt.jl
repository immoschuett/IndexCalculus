function horner_evaluate_var_crt(f,TA,A::SMat{nmod},c) 
	# THIS IS STILL UNUSED. TODO 
    #return f(A^t *A)*c
	T = elem_type(base_ring(A))
	(n,m) = size(A)
	storing_n = Vector{T}(undef,n)
    s =  Vector{T}(undef,m)
	C = collect(coefficients(f))
	n = length(C)
	s =  Hecke.mul!(s,TA, Hecke.mul!(storing_n,A,c)).*C[end]+c.*C[end-1]
	for i = n-2:-1:1
		#s = A^t * A * s + fi * c  inloop
		s =  Hecke.mul!(s,TA, Hecke.mul!(storing_n,A,s)) + c.*C[i]
	end
	@debug iszero(f(transpose(Matrix(A))*Matrix(A))*c - s) ? (@info "HORNER: f(A^t*A)c = s",true) : (@error  "HORNER: f(A^t*A)c = s",false)
	return s
end
function wiedemann_var_crt(A::SMat{fmpz_mod}) #N::fmpz || N::Int64
    ##########################################################################################################################################
	RR = base_ring(A)
	N = modulus(RR)
	T = fmpz
	A = change_base_ring(ZZ, A)::SMat{T}
	(n,m) = nrows(A),ncols(A);
    ##########################################################################################################################################
	@debug begin 
		r = rank(Matrix(A))
		r == (m-1) ? nothing : @warn("WIEDEMANN: rank A small with rank(A) =$r != m-1 = $(m-1)")
	end
    #TODO outsource to PreProcessing
	TA = transpose(A)
	##########################################################################################################################################
	# Prealloc +Randomchoice
    r = fmpz.(rand(RR, m))
	c = fmpz.(rand(RR, m))
    randlin = rand_srow(min(m,10),m,min(10,N),ZZ)
	seq = Vector{T}(undef, 2*n)
	storing_n = zeros(T,n)#Vector{T}(undef,n)
    storing_m = zeros(T,m)#Vector{T}(undef,m)
	z = zero(T)
	##########################################################################################################################################

    #Wiedemann sequence
    # solve A^tAx = y2 => x -y in kernel(A^tA) to avoid finding zero vec
	##########################################################################################################################################
	M = max(maximum(A),-minimum(A))
	B = n*(m*M*N)^2
	tr = nthreads()
	t = root(B,tr)
	p_0 = next_prime(t)
	P = [p_0]
	while prod(P) < B
		push!(P,next_prime(P[end]))
	end 
	t = length(P)
	@debug t == nthreads() 
	##########################################################################################################################################
	ACSC = sparse(A)
	TACSC = sparse(TA)
	y = mul(TA, mul(A,r))
	seq[1] = dot(randlin,c) #randlin*c 	
	Sol = zeros(fmpz,m,t)
	Sol_prealloc = [storing_m for i=1:t] #TODO 
	for i = 2:2*n  #Wiedemann sequence
		for i = 1:t
			Sol[:,i] = LinearAlgebra.mul!(Sol[:,i],TACSC, LinearAlgebra.mul!(storing_n,ACSC,c)).%P[i]
		#c = multact!(c,TA,(mulact!(storing,A,c,zero!(RR),st)),zero(RR),st) # generate sequence
		end 
		c = [ crt(Sol[i,:],P) for i in 1:m ].%N
		seq[i] = dot(randlin,c)  #eleminates
	end
    ##########################################################################################################################################
	seq = RR.(seq)
	done,f = Hecke_berlekamp_massey(seq)
	A = change_base_ring(RR,A)
	TA = change_base_ring(RR,TA)
	y = RR.(y)
	@debug begin
		degr = degree(f)
		@info "WIEDEMANN: deg f = $degr where size(A^t*A) = $m"
		typeof(f) != fmpz || (@warn "ERLEKAMP_MASSEY: f may be constant polynom")
		done || (@warn "ERLEKAMP_MASSEY: modulus N is not prime, TODO: still catch some gcds")
		iszero(f(Matrix(TA)*Matrix(A))) ? (@info "BERLEKAMP_MASSEY: valid return") : (@error "BERLEKAMP_MASSEY: unexpected return")
		#note that second case appears only for debugging storage.(since sequence beginss with A*r)
	end
	delta =0
	while iszero(coeff(f,0)) #TODO collect coeffs:
		delta+=1
		f = divexact(f,gen(parent(f)))
	end
	@debug delta<2 || (@warn "WIEDEMANN: first nonzero coeff of f is a_$delta")
	constpartoff = coeff(f,0)
	a = -inv(constpartoff)
	reducedf = divexact(f-constpartoff,gen(parent(f)))
    ##########################################################################################################################################
    #f(TA*A)'c
    v = horner_evaluate_var(reducedf,TA,A,y).*a
	@debug begin 
		iszero(mul(A,v-r)) ? (@info "WIEDEMANN: A(v-r) = 0",true) : (@error "WIEDEMANN: A(v-r) = 0",false) 
	end 
    ##########################################################################################################################################
	return (v-r)
end

####
# some multiplications:
# somehow this does not work inplace, neither correct if c is inital zero???
function multi!(c::Vector{fmpz}, A::SMat{fmpz}, b::Vector{fmpz})
    t = fmpz()
    for (i, r) in enumerate(A)
      dot!(c[i],r,b,t)
    end
    @debug mul(A,b) == c || @error "MULT for fmpz fatal error"
    return c
end
#dot fmpz
function dot!(s::fmpz, sr::SRow{fmpz}, a::Vector{fmpz},t::fmpz)
    zero!(s)
    zero!(t)
    for (i,v) = sr
      Hecke.mul!(t, v, a[i])
      Hecke.add!(s, s, t)
    end
    return s
end
function multi_ex!(c::Vector{fmpz_mod}, A::SMat{fmpz_mod}, b::Vector{fmpz_mod}) where T
    t = fmpz()
    for (i, r) in enumerate(A)
      c[i] = dot_experimental!(c[i],r,b,t)
    end
    return c
end
function dot_experimental!(s::fmpz_mod, sr::SRow{fmpz_mod}, a::Vector{fmpz_mod},t::fmpz)
    m = modulus(parent(s))
    zero!(s.data)
    zero!(t)
    for (i,v) = sr
      Hecke.mul!(t, v.data, a[i].data)
      Hecke.add!(s.data, s.data, t)
    end
    mod!(s.data, s.data, m)
    return s
end
##dot fmpz_mod
function dot!(s::fmpz_mod, sr::SRow{fmpz_mod}, a::Vector{fmpz_mod})
    m = modulus(parent(s))
    zero!(s.data)
    t = fmpz()
    for (i,v) = sr
      Hecke.mul!(t, v.data, a[i].data)
      Hecke.add!(s.data, s.data, t)
    end
    mod!(s.data, s.data, m)
    return s
end
# some more in external test documents.