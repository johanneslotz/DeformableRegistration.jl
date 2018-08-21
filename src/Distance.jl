module Distance

using DeformableRegistration
using DeformableRegistration.Transformation
using DeformableRegistration.Interpolation
using DeformableRegistration.ImageProcessing
using DeformableRegistration.Types
using SparseArrays

using Interpolations

export ssdDistance, maskedSsdDistance, ngfDistance, estimateNGFEpsilon
export calculateHessianStaggered
export ssdDistanceArbitraryGrid


include("distances/ssdDistance.jl")
include("distances/maskedSsdDistance.jl")
include("distances/ngfDistance.jl")
include("distances/estimateNGFEpsilon.jl")

end
