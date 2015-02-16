function ssdDistance(referenceImage::Image,templateImage::Image,
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
  N = prod(getSize(referenceImage))
  pixelSize = prod(getPixelSpacing(referenceImage))
  residual = Array{Float64,N}
  residual = transformedImage .- referenceImage.data[:]
  #functionValue = (0.5 .* pixelSize .* residual' * residual)[1]
  functionValue = 0.5 * pixelSize * BLAS.dot(N,residual,1,residual,1)[1]

  # derivative?
  if(doDerivative)
    if(matrixFree)
      dFunctionValue = ssdDerivativeMF(dX_transformedImage,dY_transformedImage,residual,centeredGrid,pixelSize,N,parametricOnly)
    else
      dFunctionValue,dTransformedImage = ssdDerivativeMB(dX_transformedImage,dY_transformedImage,residual,centeredGrid,pixelSize,N,parametricOnly)
    end
  else
    dFunctionValue = 0.0
    dTransformedImage = 0.0
  end

  # hessian?
  if(doHessian)
    if(matrixFree)
      if(parametricOnly)
        d2FunctionValue = computeParametricHessian(N,pixelSize,dX_transformedImage,dY_transformedImage,centeredGrid)
      else
        d2FunctionValue(x) = hessianFunction(N,pixelSize,dX_transformedImage,dY_transformedImage,x)
      end
    else
      d2FunctionValue = pixelSize .* dTransformedImage' * dTransformedImage
    end
  else
    d2FunctionValue = 0.0
  end

  return functionValue,dFunctionValue,d2FunctionValue,(dX_transformedImage,dY_transformedImage)

end

function ssdDerivativeMB(dX_transformedImage,dY_transformedImage,residual,centeredGrid,pixelSize,N,parametricOnly)
  # rearrange the image derivative
  dTransformedImage = spdiagm((dX_transformedImage,dY_transformedImage),[0,N])
  # parametric case
  if(parametricOnly)
    Q = sparse([centeredGrid[1:N] centeredGrid[N+1:end] ones(N)])
    Q = [Q spzeros(size(Q)[1],size(Q)[2]); spzeros(size(Q)[1],size(Q)[2]) Q]
    dTransformedImage = dTransformedImage * Q
  end
  return pixelSize .* dTransformedImage' * residual, dTransformedImage
end

function ssdDerivativeMF(dX_transformedImage,dY_transformedImage,residual,centeredGrid,pixelSize,N,parametricOnly)
  dFunctionValue = BLAS.scal(2*N,pixelSize,dTtransposedMultiplication(N,dX_transformedImage,dY_transformedImage,residual),1)
  # parametric case
  if(parametricOnly)
    dFunctionValue = QtransposedMultiplication(N,dFunctionValue,centeredGrid)
  end
  return dFunctionValue
end
