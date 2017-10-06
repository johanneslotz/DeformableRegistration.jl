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
    additionalOptions = Dict()
    new(parametricOnly,matrixFree,levels, interpolateToReferenceImage,
        edgeParameterR,edgeParameterT,useEdgeParameterInNumerator,
        regularizerWeight,
        maxIterGaussNewton,maxIterCG,
        additionalOptions)
  end

end
