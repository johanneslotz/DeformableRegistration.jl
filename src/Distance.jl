module Distance

using Images
using Logging

using ImageRegistration
using ImageRegistration.Transformation
using ImageRegistration.Interpolation
using ImageRegistration.ImageProcessing

export ssdDistance, maskedSsdDistance, ngfDistance, estimateNGFEpsilon
export calculateHessianStaggered

include("distances/ssdDistance.jl")
include("distances/maskedSsdDistance.jl")
include("distances/ngfDistance.jl")
include("distances/estimateNGFEpsilon.jl")

end
