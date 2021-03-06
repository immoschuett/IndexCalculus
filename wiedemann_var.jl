function wiedemann_var(A::SMat{nmod}) #N::fmpz || N::Int64
    ##########################################################################################################################################
	RR = base_ring(A)
	N = modulus(RR)
	T= elem_type(RR)
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
    r = rand(RR, m)
	c = rand(RR, m)
    randlin = rand_srow(min(m,10),m,min(10,N),RR)
	seq = Vector{T}(undef, 2*n)
	storing_n = Vector{T}(undef,n)
    storing_m =  Vector{T}(undef,m)
	z = zero(RR)
	##########################################################################################################################################

    #Wiedemann sequence
    # solve A^tAx = A^tAr = y  => x -r in kernel(A^tA) to avoid finding zero vec
	y = Hecke.mul!(storing_m,TA, Hecke.mul!(storing_n,A,r))
	seq[1] = dot(randlin,c,zero!(z)) #randlin*c 		
	for i = 2:2*n  #Wiedemann sequence
        c =  Hecke.mul!(c,TA, Hecke.mul!(storing_n,A,c))
		#c = multact!(c,TA,(mulact!(storing,A,c,zero!(RR),st)),zero(RR),st) # generate sequence
		seq[i] = dot(randlin,c)  #eleminates
	end
    ##########################################################################################################################################
	done,f = Hecke_berlekamp_massey(seq)
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
function wiedemann_var(A::SMat{fmpz_mod}) #N::fmpz || N::Int64
    ##########################################################################################################################################
	RR = base_ring(A)
	N = modulus(RR)
	T= elem_type(RR)
	A = change_base_ring(RR, A)::SMat{T}
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
    r = rand(RR, m)
	c = rand(RR, m)
    randlin = rand_srow(min(m,10),m,min(10,N),RR)
	seq = Vector{T}(undef, 2*n)
	storing_n = zeros(RR,n)#Vector{T}(undef,n)
    storing_m = zeros(RR,m)#Vector{T}(undef,m)
	z = zero(RR)
	##########################################################################################################################################

    #Wiedemann sequence
    # solve A^tAx = y2 => x -y in kernel(A^tA) to avoid finding zero vec
	y = multi!(storing_m,TA, multi!(storing_n,A,r))
	seq[1] = dot(randlin,c,zero!(z)) #randlin*c 		
	for i = 2:2*n  #Wiedemann sequence
        c =  multi!(c,TA, multi!(storing_n,A,c))
		#c = multact!(c,TA,(mulact!(storing,A,c,zero!(RR),st)),zero(RR),st) # generate sequence
		seq[i] = dot(randlin,c)  #eleminates
	end
    ##########################################################################################################################################
	done,f = Hecke_berlekamp_massey(seq)
	@debug begin
		degr = degree(f)
		@info "WIEDEMANN: deg f = $degr where size(A^t*A) = $m"
		typeof(f) != fmpz || (@warn "ERLEKAMP_MASSEY: f may be constant polynom")
		done || (@warn "BERLEKAMP_MASSEY: modulus N is not prime, TODO: still catch some gcds")
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
function horner_evaluate_var(f,TA,A::SMat{nmod},c)
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
function horner_evaluate_var(f,TA,A::SMat{fmpz_mod},c)
    #return f(A^t *A)*c
	R = base_ring(A)
	#T = elem_type(R)
	(n,m) = size(A)
	storing_n = zeros(R,n)#Vector{T}(undef,n)
    s =  zeros(R,m)#Vector{T}(undef,m)
	C = collect(coefficients(f))
	n = length(C)
	s =  multi!(s,TA, multi!(storing_n,A,c)).*C[end]+c.*C[end-1]
	for i = n-2:-1:1
		#s = A^t * A * s + fi * c  inloop
		s =  multi!(s,TA, multi!(storing_n,A,s)) + c.*C[i]
	end
	@debug iszero(f(transpose(Matrix(A))*Matrix(A))*c - s) ? (@info "HORNER: f(A^t*A)c = s",true) : (@error  "HORNER: f(A^t*A)c = s",false)
	return s
end
function rand_srow(l,n,b,R)
    #generate fmpz sparse_row, indx not greater than n limited by n
    #l values not greater than b
    val =  rand(1:b,l)
    pos = randperm!(Vector{Int}(undef, n))[1:l]
    return sparse_row(R,pos,val)
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
function multi!(c::Vector{fmpz_mod}, A::SMat{fmpz_mod}, b::Vector{fmpz_mod}) where T
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

