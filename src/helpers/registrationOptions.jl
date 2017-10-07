export regOptions

type regOptions

  # overall
  parametricOnly::Bool
  matrixFree::Bool
  levels::Array{Int,1}
  interpolateToReferenceImage::Bool

  # distance
  # ngf parameter
  edgeParameterR::Float64
  edgeParameterT::Float64
  useEdgeParameterInNumerator::Bool

  # regularizer
  regularizerWeight::Float64

  # optimization
  maxIterGaussNewton::Int
  maxIterCG::Int
  stopping::Dict

  # additional options in dict
  additionalOptions::Dict

  # set default options
  function regOptions()
    parametricOnly = false
    matrixFree = false
    levels = [5,4,3]
    interpolateToReferenceImage = false
    edgeParameterR = 0.01
    edgeParameterT = 0.01
    useEdgeParameterInNumerator = true
    regularizerWeight = 1
    maxIterGaussNewton = 30
    maxIterCG = 2000
    stopping = Dict(
        "tolJ" => 1e-3,    # tolerance: change of the objective function
        "tolY" => 1e-2,    # tolerance: change of the variables
        "tolG" => 1e-2,    # tolerance: change of the gradient
        "tolQ" => 1e-4     # tolerance: change of quotient)
        )
    additionalOptions = Dict()
    new(parametricOnly,matrixFree,levels, interpolateToReferenceImage,
        edgeParameterR,edgeParameterT,useEdgeParameterInNumerator,
        regularizerWeight,
        maxIterGaussNewton,maxIterCG, stopping,
        additionalOptions)
  end

end
