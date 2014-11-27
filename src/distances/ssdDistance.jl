function ssdDistance(referenceImage::Image,templateImage::Image,
                     transformedGrid::Array{Float64,1};
                     doDerivative=false,doHessian=false,
                     parametricOnly=false,
                     centeredGrid = getCellCenteredGrid(referenceImage))

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
  residual = transformedImage .- referenceImage.data[:]
  functionValue = (0.5 .* h .* residual' * residual)[1]

  # check if only the distance is wanted
  if(doDerivative)
      # rearrange the image derivative
      dTransformedImage = spdiagm((dX_transformedImage,dY_transformedImage),[0,N])
      # parametric case
      if(parametricOnly)
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

  return functionValue,dFunctionValue,d2FunctionValue,(dX_transformedImage,dY_transformedImage)

end

function ssdDistanceMatrixFree(referenceImage::Image,templateImage::Image,
                               transformedGrid::Array{Float64,1};
                               doDerivative=false,doHessian=false,
                               parametricOnly=false,
                               centeredGrid = getCellCenteredGrid(referenceImage))

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
  residual = transformedImage .- referenceImage.data[:]
  functionValue = 0.5 * h * BLAS.dot(N,residual,1,residual,1)[1] # !

  # check if only the distance is wanted
  if(doDerivative)
      dFunctionValue = BLAS.scal(2*N,h,dTtransposedMultiplication(N,dX_transformedImage,dY_transformedImage,residual),1)
      if(parametricOnly)
          dFunctionValue = QtransposedMultiplication(N,dFunctionValue,centeredGrid)
      end
  else
      dFunctionValue = 0
      dTransformedImage = 0
  end

  if doHessian
      if(parametricOnly)
          d2FunctionValue = computeParametricHessian(N,h,dX_transformedImage,dY_transformedImage,centeredGrid)
      else
          d2FunctionValue(x) = hessianFunction(N,h,dX_transformedImage,dY_transformedImage,x)
      end
  else
      d2FunctionValue = 0
  end

  return functionValue,dFunctionValue,d2FunctionValue,(dX_transformedImage,dY_transformedImage)

end

function dTtransposedMultiplication(N::Int,dX::Array{Float64,1},dY::Array{Float64,1},r::Array{Float64,1})
    result = zeros(2*N)
    @simd for i=1:N
        @inbounds result[i] = dX[i] * r[i];
        @inbounds result[N+i] = dY[i] * r[i];
    end
    return result
end

# Q' * hdTr with hdTr = h * dT' * residual
function QtransposedMultiplication(N::Int,hdTr::Array{Float64,1},centeredGrid::Array{Float64,1})
    result = zeros(6)
    @simd for i=1:N
        @inbounds result[1] += centeredGrid[i]*hdTr[i]
        @inbounds result[2] += centeredGrid[N+i]*hdTr[i]
        @inbounds result[3] += hdTr[i]
        @inbounds result[4] += centeredGrid[i]*hdTr[N+i]
        @inbounds result[5] += centeredGrid[N+i]*hdTr[N+i]
        @inbounds result[6] += hdTr[N+i]
    end
    return result
end

# compute parametric hessian
# result = h*(dT*Q)'*(dT*Q) = h* Q' * dT' dT * Q
function computeParametricHessian(N::Int,h::Float64,dTX::Array{Float64,1},dTY::Array{Float64,1},centeredGrid::Array{Float64,1})
    result = zeros(6,6)
    help = zeros(2*N); help[1:N] = centeredGrid[1:N]
    result[:,1] = QtransposedMultiplication(N,hessianFunction(N,h,dTX,dTY,help),centeredGrid)
    help[1:N] = centeredGrid[N+1:end]
    result[:,2] = QtransposedMultiplication(N,hessianFunction(N,h,dTX,dTY,help),centeredGrid)
    help[1:N] = ones(N)
    result[:,3] = QtransposedMultiplication(N,hessianFunction(N,h,dTX,dTY,help),centeredGrid)
    help = zeros(2*N); help[N+1:end] = centeredGrid[1:N]
    result[:,4] = QtransposedMultiplication(N,hessianFunction(N,h,dTX,dTY,help),centeredGrid)
    help[N+1:end] = centeredGrid[N+1:end]
    result[:,5] = QtransposedMultiplication(N,hessianFunction(N,h,dTX,dTY,help),centeredGrid)
    help[N+1:end] = ones(N)
    result[:,6] = QtransposedMultiplication(N,hessianFunction(N,h,dTX,dTY,help),centeredGrid)
    return result
end

# result = h*H*x, H = h*dT'*dT
function hessianFunction(N::Int,h::Float64,dTX::Array{Float64,1},dTY::Array{Float64,1},x::Array{Float64,1})
    result = zeros(2*N)
    @simd for i=1:N
        @inbounds result[i] = h*dTX[i]*dTX[i]*x[i] + h*dTX[i]*dTY[i]*x[N+i]
        @inbounds result[N+i] = h*dTX[i]*dTY[i]*x[i] + h*dTY[i]*dTY[i]*x[N+i]
    end
    return result
end
