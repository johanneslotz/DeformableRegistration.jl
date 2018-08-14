module Examples

using Images
using MicroLogging

using DeformableRegistration: ImageProcessing, Transformation, Interpolation, Distance, Regularizer, Optimization
using DeformableRegistration.regOptions
using DeformableRegistration.Types

using Interpolations

export registerImagesParametric, registerImagesNonParametric

function registerImagesParametric(referenceImage,templateImage, options::regOptions;
                                  affineParameters = [1.0,0,0,0,1.0,0],
                                  measureDistance = ssdDistance)

  affineParametersInitial = affineParameters
  options.parametricOnly = true

  # start multilevel registration
  for level in options.levels

    # define image level
    R = restrictResolutionToLevel(referenceImage,level)
    centeredGrid = getCellCenteredGrid(R);
    T = restrictResolutionToLevel(templateImage,level)

    # output
    @info "level ",level,": [",size(R.data)[1],"]x[",size(R.data)[2],"]"

    # define objective function
    Jfunc(param;doDerivative=false,doHessian=false) =
            measureDistance(R,T,transformGridAffine(centeredGrid,param).data,
                doDerivative=doDerivative,doHessian=doHessian,options=options, centeredGrid=centeredGrid.data)

    # gauss newton method
    affineParameters = optimizeGaussNewton(Jfunc,affineParameters,affineParametersInitial,options)
  end

  return affineParameters
end

function registerImagesNonParametric(referenceImage, templateImage, options::regOptions;
                                     affineParameters = [1.0,0,0,0,1.0,0],
                                     measureDistance = ssdDistance,
                                     regularizerOperator=createCurvatureOperatorCentered,
                                     initialDisplacement=scaledArray(zeros(1),(1,),[],[]))

  affineParametersInitial = affineParameters
  options.parametricOnly = false
  referenceGrid = []; deformedGrid = []; imageSize = [];

  # start multilevel registration
  for level = options.levels

  	R = restrictResolutionToLevel(referenceImage,level)
  	centeredGrid = getCellCenteredGrid(R);
  	T = restrictResolutionToLevel(templateImage,level)

    # define initial Grid
  	affineParametersInitialGrid = transformGridAffine(centeredGrid,affineParametersInitial)
  	if level==options.levels[1]
        if initialDisplacement.data == zeros(1)
  	        # initial deformedGrid
            referenceGrid = transformGridAffine(centeredGrid,affineParameters)
            deformedGrid = referenceGrid
        else
            referenceGrid = transformGridAffine(centeredGrid,affineParameters)
    	    deformedGrid = interpolateDeformationField(initialDisplacement, referenceGrid, interpolationScheme=InterpLinearFast) + referenceGrid
        end
  	else
  	      #deformation field interpolation
  	      displacementField = deformedGrid - referenceGrid
  	      referenceGrid = transformGridAffine(centeredGrid,affineParameters)
  	      deformedGrid = interpolateDeformationField(displacementField, referenceGrid, interpolationScheme=InterpLinearFast) + referenceGrid
  	end

  	imageSize = getSize(R)
  	@info "level ",level,": [",size(R.data)[1],"]x[",size(R.data)[2],"]"

  	# define objective function
  	regularizerMatrix = regularizerOperator(R.voxelsize, imageSize)
  	Jfunc(grid;doDerivative=false,doHessian=false) =
        measureDistance(R, T, grid, doDerivative=doDerivative,
            doHessian=doHessian, options=options, centeredGrid=centeredGrid.data)[1:3] +
            options.regularizerWeight * regularizer(grid-referenceGrid.data,regularizerMatrix)

  	## gauss newton method
  	deformedGrid = optimizeGaussNewton(Jfunc, deformedGrid.data, affineParametersInitialGrid.data, options)
    deformedGrid = scaledArray(deformedGrid, size(R.data), R.voxelsize, R.shift)
  end

  # return deformation field
  displacementField = deformedGrid-referenceGrid
  if options.interpolateToReferenceImage
      referenceGrid = transformGridAffine(getCellCenteredGrid(referenceImage),affineParameters)
      return interpolateDeformationField(displacementField, referenceGrid, interpolationScheme=BSpline(Linear()))
  else
      return displacementField
  end

end

end
