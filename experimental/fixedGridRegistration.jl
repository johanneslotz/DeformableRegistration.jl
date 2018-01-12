using Images
using Logging

using DeformableRegistration: ImageProcessing, Transformation, Interpolation, Distance, Regularizer, Optimization
using DeformableRegistration.regOptions
using Interpolations
import DeformableRegistration.Interpolation.interpolateArray
include("../src/helpers/objectiveFunctionCreation.jl")
include("ssdDistanceCoarseGrid.jl")

function registerNonParametricFixedGridResolution(referenceImage, templateImage, options::regOptions;
                                     affineParameters = [1.0,0,0,0,1.0,0],
                                     measureDistance = ssdDistanceArbitraryGrid,
                                     regularizerOperator=createCurvatureOperatorCentered,
                                     initialDisplacement=scaledArray(zeros(1),(1,),[],[]),
                                     constraint = 0,
                                     interpolationScheme=InterpLinearFast, gradientDescentOnly=false, gridLevel=4)#BSpline(Linear()))

  doConstraints = constraint != 0

  affineParametersInitial = affineParameters
  options.parametricOnly = false
  referenceGrid = []; deformedGrid = []; imageSize = [];
  R_level1 = restrictResolutionToLevel(referenceImage, gridLevel)
  T_level1 = restrictResolutionToLevel(templateImage,  gridLevel)
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
