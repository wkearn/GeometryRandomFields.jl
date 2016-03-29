using GeometryRandomFields
using Base.Test

gf = GaussianCovarianceFunction(10,1)
@test gf(0) == 1
@test_approx_eq gf(10) exp(-1/2)

mf = MatérnCovarianceFunction(1/2,10,1)
@test mf(0) == 1
@test_approx_eq mf(10) exp(-1)

m,n = 1024,512
M = grf(m,n,gf)
@test size(M) == (m,n)

GeometryRandomFields.acf(real(M))

Z1,Z2 = reim(M)

@test Γ(Z1,10) == 0
