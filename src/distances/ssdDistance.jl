function ssdDistance(referenceImage::Image,templateImage::Image,
                     transformedGrid::Array{Float64,1};
                     doDerivative=false,doHessian=false,
                     parametricOnly=false)

  if(checkStaggered(referenceImage,transformedGrid))
      error("StaggeredGrids are not supported. Please transform to cell centered first.")
  end

  # interpolation of the template image at transformed grid points
  transformedImage, dY_transformedImage, dX_transformedImage =
      linearImageInterpolationAtGridWithDerivative(templateImage,transformedGrid)

  # measure the ssd distance
  N = prod(size(referenceImage))
  h = prod(pixelspacing(referenceImage))
  residual = Array{Float64,N}
  residual = transformedImage[:] .- referenceImage.data[:]
  functionValue = (0.5 .* h .* residual' * residual)[1]

  # check if only the distance is wanted
  if(doDerivative)
      # rearrange the image derivative
      dTransformedImage = spdiagm((dX_transformedImage[:],dY_transformedImage[:]),[0,N])
      # parametric case
      if(parametricOnly)
          centeredGrid = getCellCenteredGrid(referenceImage)
          Q = sparse([centeredGrid[1:N] centeredGrid[N+1:end] ones(N)])
          Q = [Q spzeros(size(Q)[1],size(Q)[2]); spzeros(size(Q)[1],size(Q)[2]) Q]
          dTransformedImage = dTransformedImage * Q
      end
      # first derivative
      dFunctionValue    = h .* dTransformedImage' * residual
  else
      dFunctionValue = 0
      dTransformedImage = 0
  end

  if doHessian
      d2FunctionValue   = h .* dTransformedImage' * dTransformedImage
  else
      d2FunctionValue = 0
  end

  return functionValue,dFunctionValue,d2FunctionValue,(dX_transformedImage[:],dY_transformedImage[:])

end

