module Examples

using ImageRegistration.ImageProcessing
using ImageRegistration.Transformation
using ImageRegistration.Distance
using ImageRegistration.Regularizer
using ImageRegistration.Optimization
using Images
import Logging
using ImageRegistration

export registerImagesParametric, registerImagesNonparametric

function registerImagesParametric(referenceImage,templateImage, options::regOptions;
                                  levels = [5,4,3],
                                  measureDistance = ssdDistance)

    # initial guess affine parameters
    affineParameters = [1.0,0,0,0,1.0,0]

    # start multilevel registration
    for level in levels

        # define image level
        R = restrictResolutionToLevel(referenceImage,level)
        centeredGrid = getCellCenteredGrid(R)
        T = restrictResolutionToLevel(templateImage,level)

        # output
        Logging.info("level ",level,": [",size(R)[1],"]x[",size(R)[2],"]")

        # define objective function
        if measureDistance == ngfDistance
          JfuncWithDerivative(param) = measureDistance(R,T,options,transformGridAffine(centeredGrid,param))
          Jfunc(param) = measureDistance(R,T,options,transformGridAffine(centeredGrid,param))
        else
          options.parametricOnly = true
          optionsNoD = deepcopy(options)
          optionsD = deepcopy(options)
          optionsNoD.doDerivative = false
          optionsNoD.doHessian = false
          optionsD.doDerivative = true
          optionsD.doHessian = true
          JfuncWithDerivative(param) = measureDistance(R,T,transformGridAffine(centeredGrid,param),optionsD)
          Jfunc(param) = measureDistance(R,T,transformGridAffine(centeredGrid,param),optionsNoD)
        end


        # gauss newton method
        affineParameters = optimizeGaussNewton(Jfunc,JfuncWithDerivative,affineParameters,options)

    end

    return affineParameters

end

function registerImagesNonparametric(referenceImage,templateImage,options::regOptions;
                                     affineParameters = [1.0,0,0,0,1.0,0],
                                     levels = [4,3,2],
                                     measureDistance = ssdDistance)

    # define variables
    referenceGrid = []; deformedGrid = []; imageSize = []; spatialDomain = []

    # start multilevel registration
    for level in levels

        # define image level
        R = restrictResolutionToLevel(referenceImage,level)
        centeredGrid = getCellCenteredGrid(R)
        spatialDomain = R["spatialdomain"]
        T = restrictResolutionToLevel(templateImage,level)

        if(level==levels[1])
            # initial deformedGrid
            referenceGrid = transformGridAffine(centeredGrid,affineParameters)
            deformedGrid = referenceGrid
        else
            #deformation field interpolation
            deformationField = stg2cen(deformedGrid-referenceGrid,imageSize)
            referenceGrid = transformGridAffine(centeredGrid,affineParameters)
            deformedGrid = interpolateDeformationFieldAtGrid(deformationField,imageSize,spatialDomain,referenceGrid) + referenceGrid
        end

        # update imageSize
        imageSize = [size(R)[1],size(R)[2]]

        # centered to staggered grid
        referenceGrid = cen2stg(referenceGrid,imageSize)
        deformedGrid = cen2stg(deformedGrid,imageSize)

        # output
        Logging.info("level ",level,": [",size(R)[1],"]x[",size(R)[2],"]")

        # define objective function
        regulizerMatrix = createElasticOperatorStaggered(pixelspacing(R),imageSize)
        if (measureDistance == ngfDistance)
          optionsNoDerivative = ImageRegistration.getDefaultOptions()
          optionsNoDerivative["doDerivative"] = false
          optionsWithDerivative = ImageRegistration.getDefaultOptions()
          Jfunc(grid) = measureDistance(R,T,optionsNoDerivative,stg2cen(grid,imageSize)) + options.α * regularizer(grid-referenceGrid,regulizerMatrix)
          JfuncWithDerivative(grid) = cen2stg(measureDistance(R,T,optionsWithDerivative,stg2cen(grid,imageSize)),R) +
                                      options.α * regularizer(grid-referenceGrid,regulizerMatrix,doDerivative=true,doHessian=true)
        else
          optionsNoD = deepcopy(options)
          optionsD = deepcopy(options)
          optionsNoD.doDerivative = false
          optionsNoD.doHessian = false
          optionsD.doDerivative = true
          optionsD.doHessian = true
          Jfunc(grid) = measureDistance(R,T,stg2cen(grid,imageSize),optionsNoD) + options.α * regularizer(grid-referenceGrid,regulizerMatrix)
          JfuncWithDerivative(grid) = cen2stg(measureDistance(R,T,stg2cen(grid,imageSize),optionsD),R) +
                                      options.α * regularizer(grid-referenceGrid,regulizerMatrix,doDerivative=true,doHessian=true)
        end

        # gauss newton method
        deformedGrid = optimizeGaussNewton(Jfunc,JfuncWithDerivative,deformedGrid,options)

    end

    # return deformation field
    deformationField = stg2cen(deformedGrid-referenceGrid,imageSize)
    referenceGrid = transformGridAffine(getCellCenteredGrid(referenceImage),affineParameters)
    return interpolateDeformationFieldAtGrid(deformationField,imageSize,spatialDomain,referenceGrid)

end

end
