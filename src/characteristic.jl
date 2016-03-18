function lattice_tesselation(m,n)
    g = simple_graph(m*n,is_directed=false)
    for j in 1:n-1
        for i in 1:m-1
            ind = (j-1)*m+i
            add_edge!(g,ind,ind+1)
            add_edge!(g,ind,ind+m+1)
            add_edge!(g,ind,ind+m)
        end
        add_edge!(g,(j-1)*m+m,(j-1)*m+2m)
    end
    for i in 1:m-1
        add_edge!(g,(n-1)*m+i,(n-1)*m+i+1)
    end        
    g
end

function excursion_tesselation(A::AbstractMatrix{Bool})
    m,n = size(A)
    g = lattice_tesselation(m,n)
    ge = adjlist(find(A),is_directed=false)
    for j in 1:n
        for i in 1:m
            ind = (j-1)*m+i
            if A[i,j]
                ns = collect(out_neighbors(ind,g))
                ns = ns[ns.>ind]
                for k in ns
                    if A[k]
                        add_edge!(ge,ind,k)
                    end
                end
            end
        end
    end
    ge
end

function Γ(g)
    A = adjacency_matrix_sparse(g,Int)
    Ns = floor(Integer,trace(A^3)/6)
    Ns - num_edges(g) + num_vertices(g)
end

function Γ(A::AbstractMatrix{Bool})
    Γ(excursion_tesselation(A))
end

function Γ(Z::AbstractMatrix{Float64},u)
    Γ(Z.>u)
end
