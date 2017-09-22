function ssdDistance(referenceImage::ImageMeta,templateImage::ImageMeta,
                     transformedGrid::Array{Float64,1};
                     doDerivative::Bool=false,doHessian::Bool=false,options::regOptions=regOptions())

  # get relevant options
  parametricOnly = options.parametricOnly
  matrixFree = options.matrixFree
  centeredGrid = options.centeredGrid

  # check grid
  if(checkStaggered(referenceImage,transformedGrid))
    error("StaggeredGrids are not supported. Please transform to cell centered first.")
  end

  # check centered grid
  if(size(centeredGrid)[1]==1)
    centeredGrid = getCellCenteredGrid(referenceImage)
  end

  # interpolation of the template image at transformed grid points
  transformedImage, dY_transformedImage, dX_transformedImage =
      interpolateImage(templateImage,transformedGrid,doDerivative=true)

  # measure the ssd distance
  N = prod(getSize(referenceImage)) # product of sizes
  pixelSize = prod(getPixelSpacing(referenceImage))
  residual = transformedImage .- referenceImage.data[:]
  functionValue = 0.5 * pixelSize * BLAS.dot(N,residual,1,residual,1)[1]

  # calculate ssd derivatives matrix free?
  if(matrixFree)
    dFunctionValue,d2FunctionValue = ssdDerivativesMatrixFree(dX_transformedImage,dY_transformedImage,residual,centeredGrid,pixelSize,N,parametricOnly,doDerivative,doHessian)
  else
    dFunctionValue,d2FunctionValue = ssdDerivativesMatrixBased(dX_transformedImage,dY_transformedImage,residual,centeredGrid,pixelSize,N,parametricOnly,doDerivative,doHessian)
  end

  return functionValue,dFunctionValue,d2FunctionValue,(dX_transformedImage,dY_transformedImage)

end

function ssdDerivativesMatrixBased(dX_transformedImage::Array{Float64,1},
                                   dY_transformedImage::Array{Float64,1},
                                   residual::Array{Float64,1},
                                   centeredGrid::Array{Float64,1},
                                   pixelSize::Number,N::Number,
                                   parametricOnly::Bool,doDerivative::Bool,doHessian::Bool)

  if(doDerivative)
    dTransformedImage = spdiagm((dX_transformedImage,dY_transformedImage),[0,N])
    if(parametricOnly)
      Q = sparse([centeredGrid[1:N] centeredGrid[N+1:end] ones(N)])
      Q = [Q spzeros(size(Q)[1],size(Q)[2]); spzeros(size(Q)[1],size(Q)[2]) Q]
      dTransformedImage = dTransformedImage * Q
    end
    dFunctionValue = pixelSize .* dTransformedImage' * residual
  else
    dFunctionValue = 0.0
    dTransformedImage = 0.0
  end

  if(doHessian)
    d2FunctionValue = pixelSize .* dTransformedImage' * dTransformedImage
  else
    d2FunctionValue = 0.0
  end

  return dFunctionValue, d2FunctionValue

end

function ssdDerivativesMatrixFree(dX_transformedImage::Array{Float64,1},
                                  dY_transformedImage::Array{Float64,1},
                                  residual::Array{Float64,1},
                                  centeredGrid::Array{Float64,1},
                                  pixelSize::Number,N::Number,
                                  parametricOnly::Bool,doDerivative::Bool,doHessian::Bool)

  if(doDerivative)
    dFunctionValue = BLAS.scal(2*N,pixelSize,dTtransposedMultiplication(N,dX_transformedImage,dY_transformedImage,residual),1)
    if(parametricOnly)
      dFunctionValue = QtransposedMultiplication(N,dFunctionValue,centeredGrid)
    end
  else
    dFunctionValue = 0.0
    dTransformedImage = 0.0
  end

  if(doHessian)
    if(parametricOnly)
      d2FunctionValue = computeParametricHessian(N,pixelSize,dX_transformedImage,dY_transformedImage,centeredGrid)
    else
      d2FunctionValue(x) = hessianFunction(N,pixelSize,dX_transformedImage,dY_transformedImage,x)
    end
  else
    d2FunctionValue = 0.0
  end

  return dFunctionValue,d2FunctionValue

end
