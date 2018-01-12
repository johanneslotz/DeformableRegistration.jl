import DeformableRegistration.Distance.ssdDerivativesMatrixBased

include("../src/helpers/smoothing.jl")

function sampleArrayToGridWithSmoothing(a::scaledArray, newGrid::scaledArray)
    @assert (a.dimensions[1] > newGrid.dimensions[1] &&
        a.dimensions[2] > newGrid.dimensions[2]) "new Grid has to be coarser then current grid"

    #kernelSize = Int(ceil(a.dimensions[1]/newGrid.dimensions[1]))
    width = mean(a.dimensions)/mean(newGrid.dimensions)
    #if mod(kernelSize,2) == 0; kernelSize = kernelSize + 1; end
    #debug("kernelsize = ", kernelSize)
    #debug("width = ", width)
    adata = reshape(a.data, a.dimensions)
    #adata = smoothArray(adata, kernelSize, width)
    if width > 1.5
        print(".")
        a_int = interpolateArray((adata[1:2:end-1, 1:2:end-1]+adata[2:2:end, 2:2:end]
                                +adata[2:2:end, 1:2:end-1]+adata[1:2:end-1, 2:2:end])/4, a.voxelsize*2, a.shift, newGrid.data, newGrid.dimensions)
    else
        print("!")
        #error("should not happen now")
        a_int = interpolateArray(adata, a.voxelsize, a.shift, newGrid.data, newGrid.dimensions)
    end
    return a_int
end

function ssdDistanceArbitraryGrid(referenceImage::regImage,templateImage::regImage,
                     transformedGrid::scaledArray;
                     doDerivative::Bool=false,doHessian::Bool=false,options::regOptions=regOptions(),centeredGrid::scaledArray=scaledArray(zeros(1),(0,),[0,],[0,]))

  parametricOnly = options.parametricOnly

  # check centered grid
  if(size(centeredGrid.data)[1]==1)
    centeredGrid = getCellCenteredGrid(referenceImage)
  end

  doUpAndDownsampling = false
  #width = mean(a.dimensions)/mean(newGrid.dimensions)
  if transformedGrid.dimensions[1] < centeredGrid.dimensions[1] && transformedGrid.dimensions[2] < centeredGrid.dimensions[2]
      doUpAndDownsampling = true
  end

  if doUpAndDownsampling
      transformedGridUpsampled = interpolateDeformationField(transformedGrid-getCellCenteredGrid(transformedGrid), centeredGrid,interpolationScheme=BSpline(Constant())) + centeredGrid #interpolationScheme=BSpline(Constant())
      referenceImageData = referenceImage.data
  else
      transformedGridUpsampled = transformedGrid
      referenceImageData = interpolateImage(referenceImage,getCellCenteredGrid(transformedGrid),doDerivative=false)
  end


  # interpolation of the template image at transformed grid points
  transformedImage, dX_transformedImage, dY_transformedImage =
      interpolateImage(templateImage,transformedGridUpsampled,doDerivative=true)

  # measure the ssd distance
  residual = Array(transformedImage .- referenceImageData)[:]
  coarseCenteredGrid = getCellCenteredGrid(transformedGrid)
  #@time begin
  if doUpAndDownsampling
      debug("... resampling to grid resolution")
      dX_transformedImage  = scaledArray(dX_transformedImage, centeredGrid.dimensions, centeredGrid.voxelsize, centeredGrid.shift)
      dX_transformedImage  = sampleArrayToGridWithSmoothing(dX_transformedImage , coarseCenteredGrid)[:]
      dY_transformedImage  = scaledArray(dY_transformedImage, centeredGrid.dimensions, centeredGrid.voxelsize, centeredGrid.shift)
      dY_transformedImage  = sampleArrayToGridWithSmoothing(dY_transformedImage , coarseCenteredGrid)[:]
      residual  = scaledArray(residual,  centeredGrid.dimensions, centeredGrid.voxelsize, centeredGrid.shift)
      residual = sampleArrayToGridWithSmoothing(residual , coarseCenteredGrid)[:]
  end
  #end

  prodH = prod(transformedGrid.voxelsize)
  functionValue = 0.5 * prodH * residual' * residual
  N = prod(transformedGrid.dimensions) # product of sizes

  if(doDerivative)
      dFunctionValue,d2FunctionValue = DeformableRegistration.Distance.ssdDerivativesMatrixFree(dX_transformedImage ,dY_transformedImage ,residual ,[1.0],prodH,N,false,doDerivative,doHessian)
  else
      dFunctionValue,d2FunctionValue = 0,0
  end

  return [functionValue,dFunctionValue,d2FunctionValue,(dX_transformedImage ,dY_transformedImage )]

end
