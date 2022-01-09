using Hecke


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


function sp_preprocessing(A, l) #where l denotes the length of the original factor base
    sp_unique(A)
    TA = transpose(A)
    n,m = A.r, TA.r

    #PART 1: Eliminating one entry cols and corresponding rows + deleting rows with <2 entries and zero cols in the right part
    done = false
    while !done
        done = true
        for j = l+1:m
            if length(TA[j].pos)==1 
                done = false
                i = TA[j].pos[1]
                delete_col(A, TA, j)
                TA[j].pos = Int64[]; TA[j].values = nmod[]
                TA.nnz-=1
                if A[i].pos != []           
                    delete_col(TA,A,i)
                    A.nnz-=length(A[i].pos)
                    A[i].pos = Int64[]; A[i].values = nmod[]
                end
            end
        end
    end
    A = delete_zero_rows(A)
    TA = transpose(A)
    A = transpose(delete_zero_rows(TA,l+1))
    return A, TA
end
#TODO: implement further steps of structured Gauss and test efficiency


#Example matrix from Sieve
p = cryptoprime(10)
TESTFIELD = FField(GF(p),primitive_elem(GF(p),true))
SP = sieve_params(p,0.02,1.1)
RELMat,FB,FBx,l = Sieve(TESTFIELD,SP)
p = length(TESTFIELD.K)
modulus_ = fmpz((p-1)/2)
RR = ResidueRing(ZZ,modulus_)
A = change_base_ring(RR,RELMat)

#@time wiedemann(A, modulus_)
@time 2*3
@time wiedemann(A, modulus_)
@time wiedemann(sp_preprocessing(A, l)[1], modulus_)



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