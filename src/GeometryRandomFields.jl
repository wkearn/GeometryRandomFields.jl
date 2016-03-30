module GeometryRandomFields

export EΓ,Γ,grf,GaussianCovarianceFunction,MatérnCovarianceFunction

include("covariancefunctions.jl")
include("generators.jl")
include("spectralrepresentation.jl")
include("characteristic.jl")

end # module
