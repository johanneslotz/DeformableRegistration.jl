module Examples

using ..ImageProcessing
using ..Transformation
using ..Distance
using ..Regularizer
using ..Optimization
using Images

export registerImagesParametric, registerImagesNonparametric

function registerImagesParametric(referenceImage,templateImage;
                                  levels = [5,4,3],
                                  measureDistance = ssdDistance,
                                  maxIterGaussNewton = 10,
                                  output=false)

    # initial guess affine parameters
    affineParameters = [1.0,0,0,0,1.0,0]

    # start multilevel registration
    for level in levels

        # define image level
        R = restrictResolutionToLevel(referenceImage,level)
        centeredGrid = getCellCenteredGrid(R)
        T = restrictResolutionToLevel(templateImage,level)

        # output
        if(output)
          println("level ",level,": [",size(R)[1],"]x[",size(R)[2],"]")
        end

        # define objective function
        Jfunc(param) = measureDistance(R,T,transformGridAffine(centeredGrid,param),centeredGrid=centeredGrid)
        JfuncWithDerivative(param) = measureDistance(R,T,transformGridAffine(centeredGrid,param),parametricOnly=true,doDerivative=true,doHessian=true,centeredGrid=centeredGrid)

        # gauss newton method
        affineParameters = optimizeGaussNewton(Jfunc,JfuncWithDerivative,affineParameters,output=output)

    end

    return affineParameters

end

function registerImagesNonparametric(referenceImage,templateImage;
                                     affineParameters = [1.0,0,0,0,1.0,0],
                                     levels = [4,3,2],
                                     measureDistance = ssdDistance,
                                     alpha = 1,
                                     maxIterGaussNewton = 10,
                                     output=false)

    # define variables
    referenceGrid = []; deformedGrid = []; imageSize = []; spatialDomain = []

    # start multilevel registration
    for level in levels

        # define image level
        R = restrictResolutionToLevel(referenceImage,level)
        centeredGrid = getCellCenteredGrid(R)
        spatialDomain = R.properties["spatialdomain"]
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
        if(output)
          println("level ",level,": [",size(R)[1],"]x[",size(R)[2],"]")
        end

        # define objective function
        regulizerMatrix = createElasticOperatorStaggered(pixelspacing(R),imageSize)
        Jfunc(grid) = cen2stg(measureDistance(R,T,stg2cen(grid,imageSize)),R) +
                      alpha * regularizer(grid-referenceGrid,regulizerMatrix)
        JfuncWithDerivative(grid) = cen2stg(measureDistance(R,T,stg2cen(grid,imageSize),doDerivative=true,doHessian=true),R) +
                                    alpha * regularizer(grid-referenceGrid,regulizerMatrix,doDerivative=true,doHessian=true)

        # gauss newton method
        deformedGrid = optimizeGaussNewton(Jfunc,JfuncWithDerivative,deformedGrid,output=output)

    end

    # return deformation field
    deformationField =  stg2cen(deformedGrid-referenceGrid,imageSize)
    referenceGrid = transformGridAffine(getCellCenteredGrid(referenceImage),affineParameters)
    return interpolateDeformationFieldAtGrid(deformationField,imageSize,spatialDomain,referenceGrid)

end

end
