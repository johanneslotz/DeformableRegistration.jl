function ssdDistance(referenceImage::regImage,templateImage::regImage,
                     transformedGrid::Array{Float64,1};
                     doDerivative::Bool=false,doHessian::Bool=false,options::regOptions=regOptions(),centeredGrid::Array{Float64,1}=zeros(1))

  # get relevant options
  parametricOnly = options.parametricOnly
  matrixFree = options.matrixFree

  # check centered grid
  if(size(centeredGrid)[1]==1)
    centeredGrid = getCellCenteredGrid(referenceImage).data
  end

  # check grid
  if(checkForOddNumberOfGridPoints(referenceImage.data,transformedGrid))
    @printf("Dimension mismatch: 2* image data: %s (%s) <-> transformed grid: %s", size(referenceImage.data),size(referenceImage.data[:]), size(transformedGrid))
    error("StaggeredGrids are not supported. Please transform to cell centered first.")
  end



  # interpolation of the template image at transformed grid points
  transformedImage, dX_transformedImage, dY_transformedImage =
      interpolateImage(templateImage,transformedGrid,doDerivative=true)

  # measure the ssd distance
  N = prod(size(referenceImage.data)) # product of sizes

  residual = Array(transformedImage .- referenceImage.data)[:]
  prodH = prod(referenceImage.voxelsize)
  functionValue = 0.5 * prodH * residual' * residual
  # calculate ssd derivatives matrix free?
  if(doDerivative)
      if(matrixFree)
          dFunctionValue,d2FunctionValue = ssdDerivativesMatrixFree(dX_transformedImage,dY_transformedImage,residual,centeredGrid,prodH,N,parametricOnly,doDerivative,doHessian)
      else
          dFunctionValue,d2FunctionValue = ssdDerivativesMatrixBased(dX_transformedImage,dY_transformedImage,residual,centeredGrid,prodH,N,parametricOnly,doDerivative,doHessian)
      end
  else
      dFunctionValue,d2FunctionValue = 0,0
  end


  return [functionValue,dFunctionValue,d2FunctionValue,(dX_transformedImage,dY_transformedImage)]

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


include("../helpers/smoothing.jl")

function sampleArrayToGridWithSmoothing(a::scaledArray, newGrid::scaledArray)
    # @assert (a.dimensions[1] > newGrid.dimensions[1] &&
    #     a.dimensions[2] > newGrid.dimensions[2]) "new Grid has to be coarser then current grid"
    width = newGrid.voxelsize[1]/a.voxelsize[1]
    #print(width)
    adata = reshape(a.data, a.dimensions)
    adata = smoothArray(adata,3,width)
    #a_int = interpolateArray((adata[1:2:end-1, 1:2:end-1]+adata[2:2:end, 2:2:end]
    a_int = interpolateArray(adata, a.voxelsize, a.shift, newGrid.data, newGrid.dimensions)[1]

    return a_int
end


function ssdDistanceArbitraryGrid(referenceImage::regImage,templateImage::regImage,
                     transformedGrid::scaledArray;
                     doDerivative::Bool=false, doHessian::Bool=false,
                     options::regOptions=regOptions(),
                     interpolationScheme=BSpline(Cubic(Line())))


  parametricOnly = options.parametricOnly

  targetGrid = getCellCenteredGrid(referenceImage)

  # interpolation of the template image at transformed grid points
  transformedImage, dX_transformedImage, dY_transformedImage =
      interpolateImage(templateImage, transformedGrid, targetGrid,
        doDerivative=doDerivative, interpolationScheme=interpolationScheme)

  # measure the ssd distance
  residual = (transformedImage .- referenceImage.data.data)[:]
  prodh = prod(targetGrid.voxelsize)
  functionValue = 0.5 * prodh * sum(residual.^2)


  coarseCenteredGrid = getCellCenteredGrid(transformedGrid)

  residual  = scaledArray(residual,  targetGrid.dimensions, targetGrid.voxelsize, targetGrid.shift)
  residual = sampleArrayToGridWithSmoothing(residual , coarseCenteredGrid)[:]
  prodH = prod(coarseCenteredGrid.voxelsize)

  if doDerivative
      dX_transformedImage  = scaledArray(dX_transformedImage, targetGrid.dimensions, targetGrid.voxelsize, targetGrid.shift)
      dX_transformedImage  = sampleArrayToGridWithSmoothing(dX_transformedImage , coarseCenteredGrid)[:]
      dY_transformedImage  = scaledArray(dY_transformedImage, targetGrid.dimensions, targetGrid.voxelsize, targetGrid.shift)
      dY_transformedImage  = sampleArrayToGridWithSmoothing(dY_transformedImage , coarseCenteredGrid)[:]
  end

  if(doDerivative)
      N = prod(coarseCenteredGrid.dimensions) # product of sizes
      dFunctionValue,d2FunctionValue = DeformableRegistration.Distance.ssdDerivativesMatrixFree(dX_transformedImage, dY_transformedImage ,residual ,[1.0],prodH,N,false,doDerivative,doHessian)
      return [functionValue, dFunctionValue, d2FunctionValue, (dX_transformedImage, dY_transformedImage)]
  end

  return [functionValue, 0, 0, (0, 0)]

end
