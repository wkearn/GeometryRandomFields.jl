randc(m,n) = randn(m,n)+im*randn(m,n)

function grf(m,n,C::StationaryCovarianceFunction)
    tx = collect(0:m-1)
    ty = collect(0:n-1)
    Rjs = zeros(n,m)
    for i in 1:m, j in 1:n
        Rjs[j,i] = C([tx[1],ty[1]],[tx[i],ty[j]])
    end
    Sjs = [Rjs Rjs[:,end] Rjs[:,end:-1:2]]
    Sjs2 = [Sjs ; Sjs[end,:] ; Sjs[end:-1:2,1] Sjs[end:-1:2,end:-1:2]]
    Γ = fft(Sjs2)
    Z = randc(size(Sjs2)...)    
    fft(sqrt(Γ./(4*m*n)).*Z)    
end

function grf2(m,n,C::StationaryCovarianceFunction)
    a = [0,0.0]''
    X = gridpoints(m,n)
    D = tril(reshape(pairwise(Euclidean(),a,X),m,n))
    covdist!(D,C,m,n)
    D = Symmetric(D,:L)
    S = [D D[:,end] D[:,end:-1:2]]
    S2 = [S ; S[end,:] ; S[end:-1:2,1] S[end:-1:2,end:-1:2]]+0im
    Γ = fft!(S2)
    Z = randc(size(S2)...)
    Q = @fastmath sqrt(Γ./(4*m*n)).*Z
    fft!(Q)[1:1024,1:1024]
end

function covdist!(R,D,C::StationaryCovarianceFunction,m,n)
    for j in 1:m
        for i in j:n
            R[i,j] = C(D[i,j])
        end
    end
    R
end

function covdist!(D,C::StationaryCovarianceFunction,m,n)
    for j in 1:m
        for i in j:n
            D[i,j] = C(D[i,j])
        end
    end
    D
end

gridpoints(m,n) = gridpoints!(zeros(Int,2,m*n),m,n)

function gridpoints!(A,m,n)
    for j in 1:n
        for i in 1:m
            A[1,(j-1)*m+i] = j-1
            A[2,(j-1)*m+i] = i-1
        end
    end
    A
end
