using Hecke

ENV["JULIA_DEBUG"] = "all"

function delete_row(A, i) 
    non_zeros = length(A[i].pos)
    deleteat!(A.rows, i)
    A.r-=1
    A.nnz-=non_zeros
    return A
end


function delete_rows(A, I, sorted=true) #elements in I need to be ascending
    if !sorted
        sort(I)
    end
    for i in length(I):-1:1
        delete_row(A, I[i])
    end
    return A
end


function delete_zero_rows(A, s=1) #where s denotes the first column where we wanna start
    for i=A.r:-1:s
        if A[i].pos == []
            deleteat!(A.rows, i); A.r-=1
        end
    end
    return A
end


function delete_small_rows(A, s=1)
    for i=A.r:-1:s
        if length(A[i].pos) < 2 
            deleteat!(A.rows, i); A.r-=1
        end
    end
    return A
end


function delete_col(A, TA, j) #only deletes entries in column j, output same size as input
    for row in TA[j].pos 
        i = findall(x->x==j, A[row].pos)
        deleteat!(A[row].pos, i) ; deleteat!(A[row].values, i)
    end
    A.nnz -=length(TA[j].pos)
    return A
end


function delete_cols(A, TA, J)
    for j in J
        delete_col(A, TA, j)
    end
    return A
end


function sp_unique(A)
    out = []
    seen = Set()
   
    for i = 1:A.r
        if !in(A[i], seen)
            push!(seen, A[i])
            push!(out, A[i])
        else
            A.nnz-=length(A[i].pos)
        end
    end
    
    A.rows = out
    A.r = length(out)
    return A
end


function modify_col_add(A, B, c) #B = c*A + B
end


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

#TODO: implement further steps of structured Gauss and test efficiency
function sp_preprocessing_2(A, TA, l)
    done = false
    while !done
        done = true
        for j=l+1:A.c
            if length(TA[j]) == 2
                done = false
                a, b = TA[j].pos #indices of rows
                v, w = TA[j].values #A[a,j], A[b,j]
                p = length(A[a]); q = length(A[b]) #number of entries in rows
                if p  > q        #add A to B -> B = B-w/v * A
                    A.rows[b] = add_scaled_row(A.rows[a],A.rows[b], -divexact(w,v))
                    A.nnz+=(length(A[b])-q)
                    delete_row(A,a)
                else
                    A.rows[a] = add_scaled_row(A.rows[b],A.rows[a], -divexact(v,w))
                    A.nnz+=(length(A[a])-p)
                    delete_row(A, b)
                end
                TA = transpose(A) #TODO: TA anpassen
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

function entries(A, TA, l)
    one_entry = []
    two_entries = []
    three_entries = []
    one_row = []
    for j in l+1:TA.r
        if length(TA[j]) == 1
            push!(one_entry, j)
        elseif length(TA[j]) == 2
            push!(two_entries, j)
        elseif length(TA[j]) == 3
            push!(three_entries, j)
        end
    end 
    for i in 1:A.r
        if length(A[i]) == 1
            push!(one_row, i)
        end
    end
    return one_entry, two_entries, three_entries, one_row
end





#Example matrix from Sieve
#= include("FB_logs.jl")
p = cryptoprime(10)
TESTFIELD = BigFField(GF(p),primitive_elem(GF(p),true))
SP = sieve_params(p,0.02,1.1)
RELMat,FB,FBx,l = Sieve(TESTFIELD,SP)
p = length(TESTFIELD.K)
modulus_ = fmpz((p-1)/2)
RR = ResidueRing(ZZ,modulus_)
A = change_base_ring(RR,RELMat)
TA = transpose(A)
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
RR = ResidueRing(ZZ,6)

A = [RR(2) RR(0) RR(0);RR(4) RR(0) RR(0);RR(0) RR(0) RR(0);RR(0) RR(1) RR(0)]
A = sparse_matrix(A)
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
=#