function maskedSsdDistance(referenceImage::Array{Float64,2},templateImage::Array{Float64,2},
                           transformedGrid::Array{Float64,1},
                           mask::Array{Float64,1},
                           o::regOptions)

    # interpolation of the template image at transformed grid points
    transformedImage, dY_transformedImage, dX_transformedImage =
      linearImageInterpolationAtGridWithDerivative(templateImage,transformedGrid)

    dY_transformedImage = mask[:] .* dY_transformedImage[:]
    dX_transformedImage = mask[:] .* dX_transformedImage[:]

    # measure the ssd distance
    N = prod(size(referenceImage))
    h = prod(pixelspacing(referenceImage))
    residual = Array{Float64,N}
    residual = mask[:] .* (transformedImage[:] .- referenceImage.data[:])
	  functionValue = 0.5 .* h .* residual' * residual

    # check if only the distance is wanted
    if(~o.doDerivative)
	      return functionValue
    end

    # rearrange the image derivative
    dTransformedImage = spdiagm((dX_transformedImage[:],dY_transformedImage[:]),[0,N])


    # first derivative
    dFunctionValue    = h .* dTransformedImage' * residual

    # approximate second derivative
    d2FunctionValue   = h .* dTransformedImage' * dTransformedImage

    return functionValue,dFunctionValue,d2FunctionValue
end
