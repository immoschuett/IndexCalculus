using Hecke

#implemented functions I couldn't use
function add_scaled_row!(A::SMat{T}, i::Int, j::Int, c::T) where T
    A.nnz = A.nnz - length(A[j])
    A.rows[j] = add_scaled_row(A[i], A[j], c)
    A.nnz = A.nnz + length(A[j])
    return A
  end
  

function add_scaled_col!(A::SMat{T}, i::Int, j::Int, c::T) where T #A.nnz was not adapted
    @assert c != 0
  
    @assert 1 <= i <= ncols(A) && 1 <= j <= ncols(A)  
  
    for r in A.rows
      if i in r.pos
        i_i = findfirst(isequal(i), r.pos) #changed
        val_i = r.values[i_i]
        if j in r.pos
          i_j = findfirst(isequal(j), r.pos) #changed
          val_j = r.values[i_j]
  
          r.values[i_j] += c*r.values[i_i]
        else
          k = searchsortedfirst(r.pos, j)
          insert!(r.pos, k, j)
          insert!(r.values, k, c*r.values[i_i])
          A.nnz+=1  #added
        end
      end
    end
    return A
  end


function add_scaled_col_2!(A, TA, a, b, c) #A[_,b] -> A[_,b] + c*A[_,a]    FEHLER!!!
    pos_new_b = sort(union(TA[a].pos, TA[b].pos))   #TODO: find or write function that unites sorted lists
    pos_b = TA[b].pos
    val_b = TA[b].values
    A.nnz-=length(val_b)
    empty!(A[b].pos); empty!(A[b].values)
    for i in pos_new_b
        push!(A[b].pos, i)
        if i in TA[a].pos && i in pos_b
            push!(A[b].values, (TA[b, i]+c*TA[a, i]))
        elseif i in TA[a].pos
            push!(A[b].values, c*TA[a, i])
        else
            push!(A[b].values, TA[b, i])
        end
    end
    A.nnz+=length(pos_new_b)
end


function add_scaled_col_trans!(A, TA, i, j, c) #A[_j]->c*A[_,i]+A[_j]
    @assert c != 0
    @assert 1 <= i <= TA.r && 1 <= j <= TA.r

    for idx in TA[i].pos #indiziert Zeilen, die Eintrag iin Spalte i haben
        idx_i = findfirst(isequal(i), A[idx].pos) #indiziert i in position Array von A[idx]
        if idx in TA[j].pos
            idx_j = findfirst(isequal(j), A[idx].pos) #indiziert j in position Array von A[idx]
            A[idx].values[idx_j] += c*A[idx].values[idx_i]
        else
            k = searchsortedfirst(A[idx].pos, j)
            insert!(A[idx].pos, k, j)
            insert!(A[idx].values, k, c*A[idx].values[idx_i])
            A.nnz+=1
        end
    end
    return A
end


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
        i = searchsortedfirst(A[row].pos, j)
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



# mini examples
RR = ResidueRing(ZZ,10)

A = [RR(2) RR(0) RR(0);RR(4) RR(3) RR(0);RR(0) RR(0) RR(0);RR(0) RR(1) RR(0)]
A = sparse_matrix(A)
TA = transpose(A)
add_scaled_col_trans!(A, TA, 1, 2, RR(2))