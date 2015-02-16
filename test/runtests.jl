# method tests
include("ImageProcessing.jl")
include("Transformation.jl")
include("Interpolation.jl")
include("Distance.jl")
include("Regularizer.jl")
include("Visualization.jl")

# test full parametric and nonparametric image registration
include("ParametricRegistration.jl")
include("NonparametricRegistration.jl")
include("CellCenteredRegistration.jl")
