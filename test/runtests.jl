using Base.Test, Logging
Logging.configure(level = Logging.WARNING) # set to DEBUG to see error tables

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
