module Examples

using Images
using Logging

using ImageRegistration
using ImageRegistration.ImageProcessing
using ImageRegistration.Transformation
using ImageRegistration.Interpolation
using ImageRegistration.Distance
using ImageRegistration.Regularizer
using ImageRegistration.Optimization

export registerImagesParametric, registerImagesNonparametric

function registerImagesParametric(referenceImage,templateImage, options::regOptions;
                                  affineParameters = [1.0,0,0,0,1.0,0],
                                  measureDistance = ssdDistance)

  affineParametersInitial = affineParameters
	options.parametricOnly = true

  # start multilevel registration
  for level in options.levels

    # define image level
    R = restrictResolutionToLevel(referenceImage,level)
    centeredGrid = getCellCenteredGrid(R); options.centeredGrid = centeredGrid;
    T = restrictResolutionToLevel(templateImage,level)

    # output
    Logging.info("level ",level,": [",size(R)[1],"]x[",size(R)[2],"]")

    # define objective function
    Jfunc(param;doDerivative=false,doHessian=false) = measureDistance(R,T,transformGridAffine(centeredGrid,param),doDerivative=doDerivative,doHessian=doHessian,options=options)

    # gauss newton method
    affineParameters = optimizeGaussNewton(Jfunc,affineParameters,affineParametersInitial,options)

  end

  return affineParameters

end

function registerImagesNonparametric(referenceImage,templateImage,options::regOptions;
                                     affineParameters = [1.0,0,0,0,1.0,0],
                                     measureDistance = ssdDistance)

  affineParametersInitial = affineParameters
	options.parametricOnly = false

  # define variables
  referenceGrid = []; deformedGrid = []; imageSize = []; spatialDomain = []

  # start multilevel registration
  for level in options.levels

    # define image level
    R = restrictResolutionToLevel(referenceImage,level)
    centeredGrid = getCellCenteredGrid(R); options.centeredGrid = centeredGrid;
    spatialDomain = getSpatialDomain(R)
    T = restrictResolutionToLevel(templateImage,level)

    # define initial Grid
    affineParametersInitialGrid = transformGridAffine(centeredGrid,affineParametersInitial)

    if(level==options.levels[1])
      # initial deformedGrid
      referenceGrid = transformGridAffine(centeredGrid,affineParameters)
      deformedGrid = referenceGrid
    else
      #deformation field interpolation
      deformationField = stg2cen(deformedGrid-referenceGrid,imageSize)
      referenceGrid = transformGridAffine(centeredGrid,affineParameters)
      deformedGrid = interpolateDeformationField(deformationField,imageSize,spatialDomain,referenceGrid) + referenceGrid
    end

    # update imageSize
    imageSize = getSize(R)

    # centered to staggered grid
    referenceGrid = cen2stg(referenceGrid,imageSize)
    deformedGrid = cen2stg(deformedGrid,imageSize)
    affineParametersInitialGrid = cen2stg(affineParametersInitialGrid, imageSize)

    # output
    Logging.debug("level ",level,": [",size(R)[1],"]x[",size(R)[2],"]")

    # define objective function
    regularizerMatrix = createElasticOperatorStaggered(getPixelSpacing(R),imageSize)
    Jfunc(grid;doDerivative=false,doHessian=false) = cen2stg(measureDistance(R,T,stg2cen(grid,imageSize),doDerivative=doDerivative,doHessian=doHessian,options=options),R) + options.regularizerWeight * regularizer(grid-referenceGrid,regularizerMatrix)

    # gauss newton method
    deformedGrid = optimizeGaussNewton(Jfunc,deformedGrid,affineParametersInitialGrid,options)

  end

  # return deformation field
  deformationField = stg2cen(deformedGrid-referenceGrid,imageSize)
  referenceGrid = transformGridAffine(getCellCenteredGrid(referenceImage),affineParameters)
  return interpolateDeformationField(deformationField,imageSize,spatialDomain,referenceGrid)

end

end

