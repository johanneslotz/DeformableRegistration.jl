module Distance


using ImageRegistration
using Images
using ImageRegistration.Transformation
import Logging

export ssdDistance, ssdDistanceMatrixFree, maskedSsdDistance, ngfDistance, estimateNGFEpsilon
export calculateHessianStaggered

include("distances/ssdDistance.jl")
include("distances/maskedSsdDistance.jl")
include("distances/ngfDistance.jl")
include("distances/estimateNGFEpsilon.jl")

end
