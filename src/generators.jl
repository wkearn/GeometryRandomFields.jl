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
