EΓ(x,J) = J^2/(2pi)^(3/2)*x.*exp(-x.^2/2) + 2J/(2pi)*exp(-x.^2/2) + 1/sqrt(2pi)*(1-0.5*(1+erf(x/sqrt(2))))

function Γ(Z::AbstractMatrix{Float64},u)
    m,n = size(Z)
    F,E,P = 0,0,0
    for j in 1:n-1
        for i in 1:m-1
            if Z[i,j]>u
                P += 1
                ns1 = Z[i+1,j]>u
                ns2 = Z[i+1,j+1]>u
                ns3 = Z[i,j+1]>u
                E += ns1+ns2+ns3
                F += (ns1&ns2)+(ns2&ns3)
            end
        end
        if Z[m,j]>u
            P += 1
            E += Z[m,j+1]>u
        end
    end
    for i in 1:m-1
        if Z[i,n]>u
            P += 1
            E += Z[i+1,n]>u
        end
    end
    P += Z[m,n]>u
    F-E+P
end

function Γ(Z::AbstractMatrix{Float64},u::AbstractVector)
    γs = zeros(Int,length(u))
    for i in 1:length(u)
        γs[i] = Γ(Z,u[i])
    end
    γs
end
