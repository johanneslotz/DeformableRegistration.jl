module ImageRegistration

# helping functions
include("helpers/matrixFreeHelpers.jl")
include("helpers/registrationOptions.jl");
include("helpers/objectiveFunctionCreation.jl")

# submodules
include("ImageProcessing.jl")
include("Transformation.jl")
include("Interpolation.jl")
include("Visualization.jl")

include("Distance.jl")
include("Regularizer.jl")
include("Optimization.jl")

include("Examples.jl")

end
