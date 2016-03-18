using GeometryRandomFields
using Base.Test

GRF = GeometryRandomFields

gf = GRF.GaussianCovarianceFunction(10,1)
@test gf(0) == 1
@test_approx_eq gf(10) exp(-1/2)

mf = GRF.MatérnCovarianceFunction(1/2,10,1)
@test mf(0) == 1
@test_approx_eq mf(10) exp(-1)

m,n = 1024,512
M = GRF.grf(m,n,gf)
@test size(M) == (2n,2m)

GRF.acf(real(M))
