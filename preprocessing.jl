using Hecke, Profile
include("prepro_aux_functions.jl")

ENV["JULIA_DEBUG"] = ""

function sp_preprocessing_1(A, TA=transpose(A), l) #where l denotes the length of the original factor base
    sp_unique(A)
    #TA = transpose(A)
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
                delete_col(TA,A,i)
                A.nnz-=length(A[i])
                empty!(A[i].pos); empty!(A[i].values)
            end
        end
    end
    A = delete_zero_rows(A)
    TA = transpose(A)
    A = transpose(delete_zero_rows(TA,l+1))
    @debug sum([length(A[i]) for i = 1:A.r]) == A.nnz
    #TODO: A.cols anpassen
    return A, TA
end

#@profile sp_preprocessing_1(A, l)

#TODO: implement further steps of structured Gauss and test efficiency
function sp_preprocessing_2(A, TA, l)
    n,m = A.r, TA.r
    done = false
    while !done 
        done = true
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
                    A.nnz-=length(A[a])
                    empty!(A[a].pos); empty!(A[a].values)
                else             #add B to A -> A = A - v/w*B
                    add_scaled_col_trans!(TA, A, b, a, -divexact(v,w))
                    delete_col(TA, A, b)
                    add_scaled_row!(A, b, a, -divexact(v,w))
                    A.nnz-=length(A[b])
                    empty!(A[b].pos); empty!(A[b].values)
                end 
            end
        end
    end
    A = delete_zero_rows(A)
    TA = transpose(A)
    A = transpose(delete_zero_rows(TA,l+1))
    @debug sum([length(A[i]) for i = 1:A.r]) == A.nnz
    return A, TA
end


function sp_preprocessing_3(A, TA, l) #without comparison
    n,m = A.r, TA.r
    done = false
    while !done
        done = true
        for j = 1:m
            if length(TA[j]) == 3
                done = false
                (p,u),(q,v),(r,w) = sort([(TA[j].pos[i],TA[j].values[i]) for i=1:3], by=x->length(A[x[1]]), rev=true)
                add_scaled_col_trans!(TA, A, p, q, -divexact(v,u)) #add P to Q -> Q = Q - v/u *P
                add_scaled_col_trans!(TA, A, p, r, -divexact(w,u)) #add P to R -> R = R - w/u *P
                delete_col(TA, A, p)
                add_scaled_row!(A, p, q, -divexact(v,u))
                add_scaled_row!(A, p, r, -divexact(w,u))
                A.nnz-=length(A[p]) 
                empty!(A[p].pos); empty!(A[p].values)
            end
        end 
    end
    A = delete_zero_rows(A)
    TA = transpose(A)
    A = transpose(delete_zero_rows(TA,l+1))    
    @debug sum([length(A[i]) for i = 1:A.r]) == A.nnz  
    return A, TA
end


function sp_preprocessing_4(A, TA, l) #without comparison
    n,m = A.r, TA.r
    done = false
    while !done
        done = true
        for j = 1:m
            if length(TA[j]) == 4
                done = false
                (p,u),(q,v),(r,w),(s,x) = sort([(TA[j].pos[i],TA[j].values[i]) for i=1:4], by=x->length(A[x[1]]), rev=true)
                add_scaled_col_trans!(TA, A, p, q, -divexact(v,u)) #add P to Q -> Q = Q - v/u *P
                add_scaled_col_trans!(TA, A, p, r, -divexact(w,u)) #add P to R -> R = R - w/u *P
                add_scaled_col_trans!(TA, A, p, s, -divexact(x,u)) #add P to S -> S = S - x/u *P
                delete_col(TA, A, p)
                add_scaled_row!(A, p, q, -divexact(v,u))
                add_scaled_row!(A, p, r, -divexact(w,u))
                add_scaled_row!(A, p, s, -divexact(x,u))
                A.nnz-=length(A[p]) 
                empty!(A[p].pos); empty!(A[p].values)
            end
        end 
    end
    A = delete_zero_rows(A)
    TA = transpose(A)
    A = transpose(delete_zero_rows(TA,l+1)) 
    @debug sum([length(A[i]) for i = 1:A.r]) == A.nnz     
    return A, TA
end


function sp_preprocessing_5(A, TA, l) #without comparison
    n,m = A.r, TA.r
    done = false
    while !done
        done = true
        for j = 1:m
            if length(TA[j]) == 5
                done = false
                (p,u),(q,v),(r,w),(s,x),(t,y) = sort([(TA[j].pos[i],TA[j].values[i]) for i=1:5], by=x->length(A[x[1]]), rev=true)
                add_scaled_col_trans!(TA, A, p, q, -divexact(v,u)) #add P to Q -> Q = Q - v/u *P
                add_scaled_col_trans!(TA, A, p, r, -divexact(w,u)) #add P to R -> R = R - w/u *P
                add_scaled_col_trans!(TA, A, p, s, -divexact(x,u)) #add P to S -> S = S - x/u *P
                add_scaled_col_trans!(TA, A, p, t, -divexact(y,u)) #add P to T -> T = T - y/u *P
                delete_col(TA, A, p)
                add_scaled_row!(A, p, q, -divexact(v,u))
                add_scaled_row!(A, p, r, -divexact(w,u))
                add_scaled_row!(A, p, s, -divexact(x,u))
                add_scaled_row!(A, p, t, -divexact(y,u))
                A.nnz-=length(A[p]) 
                empty!(A[p].pos); empty!(A[p].values)
            end
        end 
    end
    A = delete_zero_rows(A)
    TA = transpose(A)
    A = transpose(delete_zero_rows(TA,l+1))      
    @debug sum([length(A[i]) for i = 1:A.r]) == A.nnz
    return A, TA
end


 function sp_preprocessing_cases(A, l)#doesn't work
    sp_unique(A)
    TA = transpose(A)
    n,m = A.r, TA.r
    done = false
    while !done 
        done = true
        for j = l+1:m
            if length(TA[j])==1 
                done = false
                i = TA[j].pos[1]           
                delete_col(TA,A,i)
                A.nnz-=length(A[i])
                empty!(A[i].pos); empty!(A[i].values)
            elseif length(TA[j]) == 2
                done = false
                a, b = TA[j].pos #indices of rows
                v, w  = TA[j].values #A[a,j], A[b,j]
                p = length(A[a]); q = length(A[b]) #number of entries in rows
                if p  > q        #add A to B -> B = B - w/v *A
                    add_scaled_col_trans!(TA, A, a, b, -divexact(w,v))
                    delete_col(TA, A, a)
                    add_scaled_row!(A, a, b, -divexact(w,v))
                    A.nnz-=length(A[a])
                    empty!(A[a].pos); empty!(A[a].values)  
                else             #add B to A -> A = A - v/w*B
                    add_scaled_col_trans!(TA, A, b, a, -divexact(v,w))
                    delete_col(TA, A, b)
                    add_scaled_row!(A, b, a, -divexact(v,w))
                    A.nnz-=length(A[b])
                    empty!(A[b].pos); empty!(A[b].values)
                end 
            end
        end
    end
    A = delete_zero_rows(A)
    TA = transpose(A)
    A = transpose(delete_zero_rows(TA,l+1))
    return A, TA
end

function sp_preprocessing_0(A, TA, l)
    modulus_ = modulus(base_ring(A))
    for i=1:A.r
        if A[i].pos[end]>l
            e = searchsortedfirst(A[i].pos, l+1)#position of first entry on the right in pos array
            v = A[i].values[e]
            if v == (modulus_ - 2)
                idx_col = A[i].pos[e] #index of col
                for idx_row in TA[idx_col]
                    print(idx_row)
                    if idx_row!=i
                        f = findfirst(isequal(idx_col), A[idx_row].pos) #position of this entry in row A[idx_row]
                        w = A[idx_row].values[f]
                        add_scaled_col_trans(TA, A, i, idx_row, -divexact(w, v))
                        add_scaled_row!(A, i, idx_row, -divexact(w, v))
                    end
                end
            end
        end
    end 
    TA = transpose(A)
    return A, TA         
end 



#TODO: schauen, wie Funktionen sinnnvoll zusammengesetzt werden
#Entweder nach Fällen zweimal durchlaufen oder direkt beide Fälle je Spalte testen






#Example matrix from Sieve
using Markdown, Nemo
include("Magma_sieve.jl")
include("wiedemann.jl")
#include("FB_logs.jl")
p = cryptoprime(5)
TESTFIELD = BigFField(GF(p),primitive_elem(GF(p),true))
SP = sieve_params(p,0.02,1.1)
RELMat,FB,FBx,l = Sieve(TESTFIELD,SP)
p = length(TESTFIELD.K)
modulus_ = fmpz((p-1)/2)
RR = ResidueRing(ZZ,modulus_)
A = change_base_ring(RR,RELMat)
density(A)
A, TA = sp_preprocessing_1(A, l)
density(A)
A, TA = sp_preprocessing_2(A, TA, l)
density(A)
A, TA = sp_preprocessing_3(A, TA, l)
density(A)
A, TA = sp_preprocessing_4(A, TA, l)
density(A)
A,TA = sp_preprocessing_5(A, TA, l)
density(A)
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