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
	#A = mod_sym!(A,N)
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
	M = 2'maximum(abs,A)
	println(M)
	B = n*(m*M*N)^2
	B = 2*n*m^2*N^4
	tr = nthreads()
	t = root(B,tr)  #WARNING size here
	p_0 = next_prime(t)
	P = [p_0]
	while prod(P) < B
		push!(P,next_prime(P[end]))
	end 
	t = length(P)
	E_p = crt_env(P)
	P2 = Int.(P) 
	@debug t == nthreads() 
	##########################################################################################################################################
	ACSC = Int.(sparse(A))
	TACSC = Int.(sparse(TA))
	y = mul(TA, mul(A,r))
	seq[1] = dot(randlin,c) #randlin*c 	
	Sol = zeros(Int64,t,m)
	Sol_2 = zeros(Int64,t,m)
	#Sol_3 = zeros(fmpz,n,t)
	Sol_3 = zeros(Int64,t,n)
	for i = 2:2*n  #Wiedemann sequence
		#c = [c c c c]
		c0 = deepcopy(c)
		#Sol = LinearAlgebra.mul!(Sol,TACSC, LinearAlgebra.mul!(Sol_3,ACSC,c))
		#Sol = transpose(Sol).%P
		#Sol = transpose(Sol)
		c = transpose([c c c c])  # TODO to t times h or c-cat 
		#inplace_conv2int(Sol_2,Sol).%P2
		Sol = Int.(c.%P2)
		Sol_2 = LinearAlgebra.mul!(Sol_2,LinearAlgebra.mul!(Sol_3,Sol,TACSC),ACSC) #TODO with Int64 inplace modular here
		#Sol = Sol.%P2
		#c = multact!(c,TA,(mulact!(storing,A,c,zero!(RR),st)),zero(RR),st) # generate sequence
		Sol = fmpz.(Sol_2)
		#Sol = inplace_conv2fmpz(Sol,Sol_2)
		c = [ crt(Sol[:,i],E_p) for i in 1:m ].%N    # do this inplace here
		#@debug c == mul(TA,mul(A,c0)) ? nothing : (@error "SEQ_wrong at $i" ,break)
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

function mul_ptr!(c::Vector{fmpz}, A::SMat{fmpz}, b::Vector{fmpz})
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

#a = fpmz(4)
#a.d 
#typeof(a.d)
# all boundet by 2^62 -1 if positiv (Int64 negative will result in 1 at 2nd bit) as negative Int64 architecture
#mum prepro... ()  bounded by 2^ * ?  possible.
using Hecke
a = fmpz(100)
while a.d == a 
	a *=2
end 
c = 1000000000000000
a = fmpz(4599999999999999999)
a+=c
b = fmpz(-4599999999999999999)
b.d

bitstring(2^62-1 + 2^62 + 2^63 +1 )
2^65
bitstring(2^62-1)

@time a = fmpz(2^62-1)
@time a.d
@time a+=1
@time a.d * a.d

a = fmpz(100)
@time a*a

@time 100 * 100

@time mul!(a.d,a.d,a.d)
@time a.d *= 2^10

a
c = fmpz

BigInt(2)^1000
a = fmpz(100)
b = fmpz(10)
c = fmpz()
function mul_ptr(c,a,b)
	c.d = a*b.d 
	return c 
end 

function mul_int(c,a,b)
	c = a*b 
	return c 
end 

@time  for i = 1:100000 mul_ptr(c,a,b) end 

@time for i =  1:100000 1231231234%41 end

function inplace_conv2fmpz(v::Vector{fmpz},w::Vector{Int64}) # if small
	for i = 1:length(w)
		v[i].d = w[i]
	end 
	return v
end 
function inplace_conv2fmpz(v::Matrix{fmpz},w::Matrix{Int64},size) # if small
	for i = 1:size[1]
		for j = 1:size[2]
			v[i,j].d = w[i,j]
			println(v[i,j])
		end
	end 
	return v
end 
function inplace_conv2int(v::Vector{Int64},w::Vector{fmpz})  # if small
	for i = 1:length(w)
		v[i] = w[i].d
	end 
	return v
end 
function inplace_conv2int(v::Matrix{Int64},w::Matrix{fmpz})  # if small
	for i = 1:length(w)
		v[i] = w[i].d
	end 
	return v
end 
function inplace_mod_conv2int(v::Vector{Int64},w::Vector{fmpz},modulus::fmpz,z::fmpz) # if small
	for i = 1:length(w)
		zero!(z)
		v[i] = mod!(z,w[i],modulus).d
	end 
	return v
end 
function inplace_modvec_conv2intmatrix(v::Matrix{Int64},w::Vector{fmpz},modulus::Vector{fmpz},z::fmpz) # if small
	for j = 1:length(modulus)
		for i = 1:length(w)
			zero!(z)
			@assert mod!(z,w[i],modulus[j]).d == mod!(z,w[i],modulus[j])
			v[i,j] = mod!(z,w[i],modulus[j]).d
		end 
	end 
	return v
end 



@time inplace_modvec_conv2intmatrix(TestMat,W,g,z)


@time inplace_conv2int(V,W)

V
TestMat = zeros(Int,3,3)
g = [fmpz(3),fmpz(2),fmpz(10)]
z = fmpz()


inplace_mod_conv2int(V,W,g,z)

W = [fmpz(123) fmpz(122) fmpz(18); fmpz(2) fmpz(312) fmpz(1)]
V = zeros(Int64,3)
D = [1, 232, -4]
Z = zeros(fmpz,3)
@time K = inplace_conv2fmpz(Z,D)

U  = zeros(Int64,4,4)
V = zeros(fmpz,4,4)
W = [1 2 3 4; 1 4 55 4; 1 4 555 5; 12 3 4 1]
W = fmpz.(W)
inplace_conv2int(U,W)


V = inplace_conv2fmpz(V,W)
[2, 3]
A = [1 3;  3 4].%[2 ,3]

V
for i in W println(i) end

W


function inplace_conv2fmpz(v::Matrix{fmpz},w::Matrix{Int64}) # if small
	for i = 1:length(W)
		v[i].d = w[i]
		println(v[i],w[i])
	end
 	return v
end 