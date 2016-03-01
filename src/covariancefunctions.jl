abstract AbstractCovarianceFunction

abstract StationaryCovarianceFunction <: AbstractCovarianceFunction

Base.call(f::StationaryCovarianceFunction,x,y) = f(sqrt(sumabs2(x-y)))

doc"""
A Gaussian covariance function is of the form

C(x_1,x_2) = σ^2 e^{-\|x_1-x_2\|^2/2L^2}

where L is the bandwidth of the Gaussian kernel.
"""
type GaussianCovarianceFunction <: StationaryCovarianceFunction
    L::Real
    σ::Real
end

Base.call(f::GaussianCovarianceFunction,h) = f.σ^2*exp(-h^2/(2*f.L^2))

doc"""
A Matérn covariance function takes the form

C(x_1,x_2) = σ^2 \frac{2^{1-ν}}{Γ(ν)}(\sqrt{2ν}\|x_1-x_2\|/ρ)^ν K_ν(\sqrt{2ν}\|x_1-x_2\|/ρ)

where ν gives the order of the modified Bessel function of the second kind, K_ν, ρ is the bandwidth/correlation length and σ is the total covariance.
"""
type MatérnCovarianceFunction <: StationaryCovarianceFunction
    ν::Real
    ρ::Real
    σ::Real
end

Base.call(f::MatérnCovarianceFunction,h) = (h==0?f.σ^2:f.σ^2*2.0^(1-f.ν)/gamma(f.ν)*(sqrt(2*f.ν)*h/f.ρ)^f.ν*besselk(f.ν,sqrt(2*f.ν)*h/f.ρ))

