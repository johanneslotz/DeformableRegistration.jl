module Distance

using Images

using DeformableRegistration
using DeformableRegistration.Transformation
using DeformableRegistration.Interpolation
using DeformableRegistration.ImageProcessing

export ssdDistance, maskedSsdDistance, ngfDistance, estimateNGFEpsilon
export calculateHessianStaggered

include("distances/ssdDistance.jl")
include("distances/maskedSsdDistance.jl")
include("distances/ngfDistance.jl")
include("distances/estimateNGFEpsilon.jl")

end
