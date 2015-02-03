module ImageRegistration

export regOptions

type regOptions
  doDerivative::Bool
  doHessian::Bool
  edgeParameterR::FloatingPoint
  edgeParameterT::FloatingPoint
  centeredGrid::Int
  useEdgeParameterInNumerator::Bool
  parametricOnly::Bool

  # regularizer weight
  α::FloatingPoint

  maxIterGaussNewton::Int
  maxIterCG::Int


  function regOptions()
    doDerivative = false
    doHessian = false
    edgeParameterR = 0.01
    edgeParameterT = 0.01
    centeredGrid = 0
    useEdgeParameterInNumerator = true
    parametricOnly = false
    α = 1
    maxIterGaussNewton = 10
    maxIterCG=2000
    new(doDerivative,doHessian,edgeParameterR,edgeParameterT,centeredGrid,useEdgeParameterInNumerator, parametricOnly, α, maxIterGaussNewton, maxIterCG)
  end
end

# Submodules
include("Transformation.jl")
include("Visualization.jl")
include("ImageProcessing.jl")

include("Distance.jl")
include("Regularizer.jl")
include("Optimization.jl")

include("Examples.jl")

include("helpers.jl")


end
