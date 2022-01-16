function wiedemann_var(A::SMat{fmpz},N) #N::fmpz || N::Int64
    ##########################################################################################################################################
	RR = ResidueRing(ZZ,N)
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
	storing_n = Vector{T}(undef,n)
    storing_m =  Vector{T}(undef,m)
	z = zero(RR)
	##########################################################################################################################################

    #Wiedemann sequence
    # solve A^tAx = y2 => x -y in kernel(A^tA) to avoid finding zero vec
	y = mul!(storing_m,TA,mul!(storing_n,A,r))
	seq[1] = dot(randlin,c,zero!(z)) #randlin*c 		
	for i = 2:2*n  #Wiedemann sequence
        c = mul!(c,TA,mul!(storing_n,A,c))
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
function horner_evaluate_var(f,TA,A,c)
    #return f(A^t *A)*c
	T = elem_type(base_ring(A))
	(n,m) = size(A)
	storing_n = Vector{T}(undef,n)
    s =  Vector{T}(undef,m)
	C = collect(coefficients(f))
	n = length(C)
	s = mul!(s,TA,mul!(storing_n,A,c)).*C[end]+c.*C[end-1]
	for i = n-2:-1:1
		#s = A^t * A * s + fi * c  inloop
		s = mul!(s,TA,mul!(storing_n,A,s)) + c.*C[i]
	end
	@debug iszero(f(transpose(Matrix(A))*Matrix(A))*c - s) ? (@info "HORNER: f(A^t*A)c = s",true) : (@error  "HORNER: f(A^t*A)c = s",false)
	return s
end