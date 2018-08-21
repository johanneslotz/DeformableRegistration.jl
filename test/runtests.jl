using Test

#configure_logging(min_level=:warn)

@testset "DeformableRegistration" begin

@time @testset "submodules" begin
    include("ImageProcessing.jl")
    include("Transformation.jl")
    include("Interpolation.jl")
    include("Distance.jl")
    include("Regularizer.jl")
end

# test full parametric and nonparametric image registration
@time @testset "example registration cases" begin
    include("ParametricRegistration.jl")
    include("NonparametricRegistration.jl")
end

end
