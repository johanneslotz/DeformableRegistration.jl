module DeformableRegistration

# helper functions
include("helpers/matrixFreeHelpers.jl")
include("helpers/registrationOptions.jl");

# submodules
include("Types.jl")
include("Transformation.jl")
include("Interpolation.jl")
include("ImageProcessing.jl")


include("Distance.jl")
include("Regularizer.jl")
include("Optimization.jl")

include("Examples.jl")

try
    include("Visualization.jl")
catch err
    println("caught: ", err)
    println("This is 'normal' during automated testing: ")
end

end
