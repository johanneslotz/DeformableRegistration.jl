module DeformableRegistration

# helping functions
include("Types.jl")
include("helpers/matrixFreeHelpers.jl")
include("helpers/registrationOptions.jl");
#include("helpers/objectiveFunctionCreation.jl")

# submodules
include("ImageProcessing.jl")
include("Transformation.jl")
include("Interpolation.jl")

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
