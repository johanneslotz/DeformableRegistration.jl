# helpers to make addition/multiplication of distance and regularizer possible

function *(factor::Number,Regularizer::(Number,Vector,SparseMatrixCSC{Float64,Int64}))
    return factor.*Regularizer[1], factor.*Regularizer[2], factor.*Regularizer[3]
end

function +(Distance::(Number,Number,Number,(Vector,Vector)),Regularizer::(Number,Vector,SparseMatrixCSC{Float64,Int64}))
    return Distance[1].+Regularizer[1], 0, 0
end

function +(Distance::(Number,Vector,Function,(Vector,Vector)),Regularizer::(Number,Vector,SparseMatrixCSC{Float64,Int64}))
    d2J(grid) = Distance[3](grid) + Regularizer[3]*grid
    return Distance[1].+Regularizer[1], Distance[2].+Regularizer[2], d2J
end

function +(Distance::(Number,Vector,SparseMatrixCSC{Float64,Int64},(Vector,Vector)),Regularizer::(Number,Vector,SparseMatrixCSC{Float64,Int64}))
    return Distance[1].+Regularizer[1], Distance[2].+Regularizer[2], Distance[3].+Regularizer[3]
end

import ImageRegistration.Transformation
function getDefaultOptions()
  options = {"doDerivative" => true,
             "doHessian" => true,
             "edgeParameterR" => 0.01,
             "edgeParameterT" => 0.01,
             "centeredGrid" => ImageRegistration.Transformation.getCellCenteredGrid,
             "useEdgeParameterInNumerator" => true,
             "parametricOnly" => true}
  return options
end
