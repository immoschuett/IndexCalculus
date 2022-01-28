using Hecke, Profile
include("prepro_aux_functions.jl")

ENV["JULIA_DEBUG"] = ""

function sp_preprocessing_1(A, l) #where l denotes the length of the original factor base
    sp_unique(A)
    TA = transpose(A)
    n,m = A.r, TA.r
    #PART 1: Eliminating one entry cols and corresponding rows + deleting rows with <2 entries and zero cols in the right part
    #TODO: nnz. jeweils auch in TA und A ändern 
    done = false
    while !done
        done = true
        for j = l+1:m
            if length(TA[j])==1 
                done = false
                i = TA[j].pos[1]
                delete_col(A, TA, j)
                empty!(TA[j].pos); empty!(TA[j].values)
                TA.nnz-=1
                if A[i].pos != []           
                    delete_col(TA,A,i)
                    A.nnz-=length(A[i])
                    empty!(A[i].pos); empty!(A[i].values)
                end
            end
        end
    end
    A = delete_zero_rows(A)
    TA = transpose(A)
    A = transpose(delete_zero_rows(TA,l+1))
    #TODO: A.cols anpassen
    return A, TA
end

#@profile sp_preprocessing_1(A, l)

#TODO: implement further steps of structured Gauss and test efficiency
function sp_preprocessing_2(A, TA, l)
    n,m = A.r, TA.r
    done = false
    while !done 
        for j=l+1:m
            if length(TA[j]) == 2
                done = false
                a, b = TA[j].pos #indices of rows
                v, w  = TA[j].values #A[a,j], A[b,j]
                p = length(A[a]); q = length(A[b]) #number of entries in rows
                if p  > q        #add A to B -> B = B - w/v *A
                    add_scaled_col_trans!(TA, A, a, b, -divexact(w,v))
                    delete_col(TA, A, a)
                    add_scaled_row!(A, a, b, -divexact(w,v))
                    empty!(A[a].pos); empty!(A[a].values)
                    A.nnz-=length(A[a]) 
                    @debug transpose(A)==TA ? nothing : (print(d,f))
                else             #add B to A -> A = A - v/w*B
                    add_scaled_col_trans!(TA, A, b, a, -divexact(v,w))
                    delete_col(TA, A, b)
                    add_scaled_row!(A, b, a, -divexact(v,w))
                    empty!(A[b].pos); empty!(A[b].values)
                    A.nnz-=length(A[b])
                    @debug transpose(A)==TA ? nothing : (print(d,f))
                end 
            end
        end
    end
    A = delete_zero_rows(A)
    TA = transpose(A)
    A = transpose(delete_zero_rows(TA,l+1))
    return A, TA
end


function sp_preprocessing_3(A, TA, l) #without comparison
    n,m = A.r, TA.r
    done = false
    for counter = 1:10
    #while !done
        for j=1:m
            if length(TA[j])==3
                done = false
                (p,u),(q,v),(r,w) = sort([(TA[j].pos[i],TA[j].values[i]) for i=1:3], by=x->length(A[x[1]]), rev=true)
                add_scaled_col_trans!(TA, A, p, q, -divexact(v,u)) #add P to Q -> Q = Q - v/u *P
                add_scaled_col_trans!(TA, A, p, r, -divexact(w,u)) #add P to R -> R = R - w/u *P
                delete_col(TA, A, p)
                add_scaled_row!(A, p, q, -divexact(v,u))
                add_scaled_row!(A, p, r, -divexact(w,u))
                empty!(A[p].pos); empty!(A[p].values)
                A.nnz-=length(A[p]) 
            end
        end 
    end
    A = delete_zero_rows(A)
    TA = transpose(A)
    A = transpose(delete_zero_rows(TA,l+1))      
    return A, TA
end


 function sp_preprocessing_cases(A, l)
    TA = transpose(A)
    n,m = A.r, TA.r
    done = false
    d = 0
    while !done 
        done = true
        for j = l+1:m
            if length(TA[j]) == 1
                done = false
                i = TA[j].pos[1]
                delete_col(A, TA, j)
                empty!(TA[j].pos); empty!(TA[j].values)
                TA.nnz-=1
                if A[i].pos != []           
                    delete_col(TA,A,i)
                    A.nnz-=length(A[i])
                    empty!(A[i].pos); empty!(A[i].values)
                end
            elseif length(TA[j]) == 2
                done = false
                a, b = TA[j].pos
                v, w  = TA[j].values
                add_scaled_col_trans!(TA, A, a, b, -divexact(w,v))
                delete_col(TA, A, a)
                add_scaled_row!(A, a, b, -divexact(w,v))
                empty!(A[a].pos); empty!(A[a].values)
                A.nnz-=length(A[a]) 
                @assert transpose(A)==TA
            end
        end
    end
    A = delete_zero_rows(A)
    TA = transpose(A)
    A = transpose(delete_zero_rows(TA,l+1))
    return A, TA
end



#TODO: schauen, wie Funktionen sinnnvoll zusammengesetzt werden
#Entweder nach Fällen zweimal durchlaufen oder direkt beide Fälle je Spalte testen






#Example matrix from Sieve
using Markdown, Nemo
include("Magma_sieve.jl")
include("wiedemann.jl")
#include("FB_logs.jl")
p = cryptoprime(10)
TESTFIELD = BigFField(GF(p),primitive_elem(GF(p),true))
SP = sieve_params(p,0.02,1.1)
RELMat,FB,FBx,l = Sieve(TESTFIELD,SP)
p = length(TESTFIELD.K)
modulus_ = fmpz((p-1)/2)
RR = ResidueRing(ZZ,modulus_)
A = change_base_ring(RR,RELMat)
#preprocessing_cases(A, l)
#TA = transpose(A)
A, TA = sp_preprocessing_1(A, l)
#B, TB = A, TA
A, TA = sp_preprocessing_2(A, TA, l)
A, TA = sp_preprocessing_3(A, TA, l)
A, TA = sp_preprocessing_cases(A, l)

entries(A, TA, l)

A, TA = sp_preprocessing_1(A, l)
entries(A, TA, l) 
A, TA = sp_preprocessing_2(A, TA, l)
entries(A, TA, l)

#@time wiedemann(A, modulus_)
@time 2*3
@time wiedemann(A, modulus_)
@time wiedemann(sp_preprocessing_1(A, l), modulus_)



# mini examples
RR = ResidueRing(ZZ,10)

A = [RR(2) RR(0) RR(0);RR(4) RR(3) RR(0);RR(0) RR(0) RR(0);RR(0) RR(1) RR(0)]
A = sparse_matrix(A)
TA = transpose(A)
#delete_row(A, 1)
delete_rows(A, [1,4])

B = [RR(2) RR(0) RR(0);RR(4) RR(0) RR(0);RR(0) RR(0) RR(0);RR(0) RR(1) RR(0);RR(2) RR(0) RR(0);RR(0) RR(1) RR(0)]
B = sparse_matrix(B)
delete_col(B,transpose(B),1)
delete_cols(B, transpose(B), [1,3])
sp_unique(B)

C = [3 0 2 0 0 0 0 1;0 5 0 0 1 -1 0 -2;4 0 0 0 0 0 0 -1;0 1 1 0 0 0 1 1;0 0 0 0 1 0 -2 1;0 0 0 1 0 0 0 -1;3 0 2 0 0 0 0 1]
C = sparse_matrix(C)
sp_preprocessing(C, 4)
#Einträge in Zeile i löschen über A[i].pos = Int64[]; A[i].values = nmod[]

#ideas:
#column operations in left part to produce new one entry columns
#eliminate columns in right part that are multiples of others