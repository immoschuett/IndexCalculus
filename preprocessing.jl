##########################################################################################################################################
# Preprocessing for SMat{fmpz_mod/nmods}
#
include("prepro_aux_functions.jl")
function sp_preprocessing(A::SMat{T}, TA::SMat{T}, l, i=1, zero=false) where T <: Union{fmpz_mod, nmod}
    @assert 1<=i<=5
    if zero
        A, TA = sp_preprocessing_0(A, TA, l)
    end
    A, TA = sp_preprocessing_1(A, TA, l)
    if i == 1
        return A, TA
    else
        A, TA = sp_preprocessing_2(A, TA, l)
        if i == 2
            return A, TA
        else
            A, TA = sp_preprocessing_3(A, TA, l)
            if i == 3
                return A, TA
            else
                A, TA = sp_preprocessing_4(A, TA, l)
                if i == 4
                    return A, TA
                else
                    A, TA = sp_preprocessing_5(A, TA, l)
                    return A, TA
                end
            end
        end
    end
    return A, TA
end
##########################################################################################################################################
# Preprocessing 0-5
function sp_preprocessing_0(A::SMat{T}, TA::SMat{T}, l) where T <: Union{fmpz_mod, nmod}
    modulus_ = modulus(base_ring(A))
    for i=1:A.r
        if A[i].pos[end]>l
            e = searchsortedfirst(A[i].pos, l+1)#position of first entry on the right in pos array
            v = A[i].values[e]
            if v == (modulus_ - 2)
                idx_col = A[i].pos[e] #index of col
                for idx_row in TA[idx_col].pos
                    if idx_row!=i
                        f = findfirst(isequal(idx_col), A[idx_row].pos) #position of this entry in row A[idx_row]
                        w = A[idx_row].values[f]
                        add_scaled_col_trans!(TA, A, i, idx_row, -divexact(w, v))
                        add_scaled_row!(A, i, idx_row, -divexact(w, v))
                    end
                end
            end
        end
    end 
    TA = transpose(A)
    return A, TA         
end 
function sp_preprocessing_1(A::SMat{T}, TA::SMat{T}, l) where T <: Union{fmpz_mod, nmod}
    sp_unique(A)
    #TA = transpose(A)
    n,m = A.r, TA.r
    #PART 1: Eliminating one entry cols and corresponding rows + deleting rows with <2 entries and zero cols in the right part
    #TODO: nnz. jeweils auch in TA und A ??ndern 
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
function sp_preprocessing_2_origin(A::SMat{T}, TA::SMat{T}, l) where T <: Union{fmpz_mod, nmod}
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
function sp_preprocessing_2(A::SMat{T}, TA::SMat{T}, l) where T <: Union{fmpz_mod, nmod}
    n,m = A.r, TA.r
    done = false
    while !done
        done = true
        for j = l+1:m
            if length(TA[j]) == 2
                done = false
                (p,u),(q,v) = sort([(TA[j].pos[i],TA[j].values[i]) for i=1:2], by=x->length(A[x[1]]), rev=true)
                add_scaled_col_trans!(TA, A, p, q, -divexact(v,u)) #add P to Q -> Q = Q - v/u *P
                delete_col(TA, A, p)
                add_scaled_row!(A, p, q, -divexact(v,u))
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
function sp_preprocessing_3(A::SMat{T}, TA::SMat{T}, l) where T <: Union{fmpz_mod, nmod}
    n,m = A.r, TA.r
    done = false
    while !done
        done = true
        for j = l+1:m
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
function sp_preprocessing_4(A::SMat{T}, TA::SMat{T}, l) where T <: Union{fmpz_mod, nmod}
    n,m = A.r, TA.r
    done = false
    while !done
        done = true
        for j = l+1:m
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
function sp_preprocessing_5(A::SMat{T}, TA::SMat{T}, l) where T <: Union{fmpz_mod, nmod}
    n,m = A.r, TA.r
    done = false
    while !done
        done = true
        for j = l+1:m
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
#=
function sp_preprocessing_cases(A, l)#doesn't work
    sp_unique(A)
    TA = transpose(A)
    n,m = A.r, TA.r
    done = false
    while !done 
        done = true
        for j = l+1:m
            print(j,",", length(TA[j]))
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
=#
##########################################################################################################################################
#TODO: implement further steps of structured Gauss and test efficiency
#TODO: schauen, wie Funktionen sinnnvoll zusammengesetzt werden
#Entweder nach F??llen zweimal durchlaufen oder direkt beide F??lle je Spalte testen
##########################################################################################################################################
