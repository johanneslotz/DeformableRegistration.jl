using Images
using Logging

using DeformableRegistration: ImageProcessing, Transformation, Interpolation, Distance, Regularizer, Optimization
using DeformableRegistration.regOptions
using Interpolations
import DeformableRegistration.Interpolation.interpolateArray
include("../src/helpers/objectiveFunctionCreation.jl")

#include("../src/distances/ssdDistance.jl")
import DeformableRegistration.Distance.ssdDerivativesMatrixBased


function smoothArray(a::Array{Float64,1}, kernelSize::Int)
    a_smooth = zeros(size(a))
    halfKernelSize = Int((kernelSize+1)/2)
    for i = 1:halfKernelSize
        a_smooth[i] = sum(a[1:i])/i
        a_smooth[end-i+1] = sum(a[end-i+1:end])/i
    end
    for i=halfKernelSize:(length(a)-halfKernelSize)
        a_smooth[i] = sum(a[(i-halfKernelSize+1):(i+halfKernelSize-1)])/kernelSize
    end
    return a_smooth
end

function smoothArray(a::Array{Float64,2}, kernelSize::Int)
    a_smooth = zeros(size(a))
    for i = 1:size(a,1)
        a_smooth[i,:] = smoothArray(a[i,:], kernelSize)
    end
    for i = 1:size(a,2)
        a_smooth[:,i] = smoothArray(a_smooth[:,i], kernelSize)
    end
    return a_smooth
end

function sampleArrayToGridWithSmoothing(a::scaledArray, newGrid::scaledArray)
    @assert (a.dimensions[1] > newGrid.dimensions[1] &&
        a.dimensions[2] > newGrid.dimensions[2]) "new Grid has to be coarser then current grid"

    kernelSize = Int(round(a.dimensions[1]/newGrid.dimensions[1]))
    if mod(kernelSize,2) == 0; kernelSize = kernelSize + 1; end
    debug("kernelsize = ", kernelSize)
    adata = reshape(a.data, a.dimensions)
    adata = smoothArray(adata, kernelSize)

    return interpolateArray(adata, a.voxelsize, a.shift, newGrid.data, newGrid.dimensions)
end

function ssdDistanceArbitraryGrid(referenceImage::regImage,templateImage::regImage,
                     transformedGrid::scaledArray;
                     doDerivative::Bool=false,doHessian::Bool=false,options::regOptions=regOptions(),centeredGrid::scaledArray=scaledArray(zeros(1),(0,),[0,],[0,]))

  parametricOnly = options.parametricOnly

  # check centered grid
  if(size(centeredGrid.data)[1]==1)
    centeredGrid = getCellCenteredGrid(referenceImage)
  end

  transformedGridUpsampled = interpolateDeformationField(transformedGrid-getCellCenteredGrid(transformedGrid), centeredGrid) + centeredGrid

  # interpolation of the template image at transformed grid points
  transformedImage, dX_transformedImage, dY_transformedImage =
      interpolateImage(templateImage,transformedGridUpsampled,doDerivative=true)

  # measure the ssd distance
  residual = Array(transformedImage .- referenceImage.data)[:]

  if transformedGrid.dimensions[1] < centeredGrid.dimensions[1] && transformedGrid.dimensions[2] < centeredGrid.dimensions[2]
      debug("... resampling to grid resolution")
      dX_transformedImage = scaledArray(dX_transformedImage, centeredGrid.dimensions, centeredGrid.voxelsize, centeredGrid.shift)
      dX_transformedImage = sampleArrayToGridWithSmoothing(dX_transformedImage, getCellCenteredGrid(transformedGrid))[:]
      dY_transformedImage = scaledArray(dY_transformedImage, centeredGrid.dimensions, centeredGrid.voxelsize, centeredGrid.shift)
      dY_transformedImage = sampleArrayToGridWithSmoothing(dY_transformedImage, getCellCenteredGrid(transformedGrid))[:]
      residual = scaledArray(residual,  centeredGrid.dimensions, centeredGrid.voxelsize, centeredGrid.shift)
      residual = sampleArrayToGridWithSmoothing(residual, getCellCenteredGrid(transformedGrid))[:]
  end

  prodH = prod(transformedGrid.voxelsize)
  functionValue = 0.5 * prodH * residual' * residual
  N = prod(transformedGrid.dimensions) # product of sizes

  if(doDerivative)
      dFunctionValue,d2FunctionValue = DeformableRegistration.Distance.ssdDerivativesMatrixFree(dX_transformedImage,dY_transformedImage,residual,[1.0],prodH,N,false,doDerivative,doHessian)
  else
      dFunctionValue,d2FunctionValue = 0,0
  end

  return [functionValue,dFunctionValue,d2FunctionValue,(dX_transformedImage,dY_transformedImage)]

end


function registerNonParametricFixedGridResolution(referenceImage, templateImage, options::regOptions;
                                     affineParameters = [1.0,0,0,0,1.0,0],
                                     measureDistance = ssdDistanceArbitraryGrid,
                                     regularizerOperator=createCurvatureOperatorCentered,
                                     initialDisplacement=scaledArray(zeros(1),(1,),[],[]),
                                     constraint = 0,
                                     interpolationScheme=InterpLinearFast, gradientDescentOnly=false)#BSpline(Linear()))

  doConstraints = constraint != 0

  affineParametersInitial = affineParameters
  options.parametricOnly = false
  referenceGrid = []; deformedGrid = []; imageSize = [];
  R_level1 = restrictResolutionToLevel(referenceImage, options.levels[2])
  T_level1 = restrictResolutionToLevel(templateImage,  options.levels[2])
  gridVoxelSize = R_level1.voxelsize
  regularizerMatrix = regularizerOperator(gridVoxelSize, getSize(R_level1))

  if initialDisplacement.data == zeros(1)
      centeredGrid = getCellCenteredGrid(R_level1);
      referenceGrid = transformGridAffine(centeredGrid,affineParameters)
      deformedGrid = referenceGrid
  else
      referenceGrid = transformGridAffine(centeredGrid,affineParameters)
      deformedGrid = interpolateDeformationField(initialDisplacement, referenceGrid, interpolationScheme=interpolationScheme) + referenceGrid
  end


  # start multilevel registration
  for level = options.levels

  	R = restrictResolutionToLevel(referenceImage,level)
  	T = restrictResolutionToLevel(templateImage,level)
    centeredGrid = getCellCenteredGrid(R)

  	imageSize = getSize(R)
  	Logging.info("level ",level,": [",size(R.data)[1],"]x[",size(R.data)[2],"]")

  	# define objective function
  	Jfunc(grid;doDerivative=false,doHessian=false) =
        measureDistance(R, T, scaledArray(grid, size(R_level1.data), R_level1.voxelsize, R_level1.shift), doDerivative=doDerivative,
            doHessian=doHessian, options=options, centeredGrid=centeredGrid)[1:3] +
            options.regularizerWeight * regularizer(grid-referenceGrid.data,regularizerMatrix)
    fValues(grid) = [measureDistance(R, T, scaledArray(grid, size(R_level1.data), R_level1.voxelsize, R_level1.shift), doDerivative=false, doHessian=false, options=options, centeredGrid=centeredGrid)[1:3],
        options.regularizerWeight , regularizer(grid-referenceGrid.data,regularizerMatrix)]
  	## gauss newton method
    if doConstraints
        initialDisplacementCoarse = interpolateDeformationField(initialDisplacement, referenceGrid, interpolationScheme=interpolationScheme) + referenceGrid
        c(x) = constraint(x, initialDisplacementCoarse)
        deformedGrid = opt(Jfunc, deformedGrid.data, referenceGrid.data, options, constraint= c, printFunction = fValues, gradientDescentOnly=gradientDescentOnly)
    else
        deformedGrid = opt(Jfunc, deformedGrid.data, referenceGrid.data, options, printFunction = fValues, gradientDescentOnly=gradientDescentOnly)
    end
    deformedGrid = scaledArray(deformedGrid, size(R_level1.data), R_level1.voxelsize, R_level1.shift)
  end

  # return deformation field
  displacementField = deformedGrid-referenceGrid
  if options.interpolateToReferenceImage
      referenceGrid = transformGridAffine(getCellCenteredGrid(referenceImage),affineParameters)
      return interpolateDeformationField(displacementField, referenceGrid, interpolationScheme=interpolationScheme)
  else
      return displacementField
  end

end
