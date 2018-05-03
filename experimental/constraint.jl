using Images
# using Logging

using DeformableRegistration: ImageProcessing, Transformation, Interpolation, Distance, Regularizer, Optimization
using DeformableRegistration.regOptions
using Interpolations
include("../src/helpers/objectiveFunctionCreation.jl")
include("../experimental/aspin.jl")

function matchInitialDisplacementAtIndex(x::Array{Float64,1}, initialDisplacement::scaledArray, side; β = 0.01)
    index = spzeros(initialDisplacement.dimensions...)
    if (side == "left") || (side == 1)
        index[:,end-1:end] = 1
    elseif (side == "right") || (side == 2)
        index[:,1:2] = 1
    elseif (side == "boundary")
        index[[1, end],:] = 1
        index[:,[1, end]] = 1
    else
        throw("side needs to be left or right")
    end
    index = vcat(index[:], index[:])
    f = β * 0.5 * index .* (x - initialDisplacement.data).^2
    dF = β * spdiagm(index .* (x - initialDisplacement.data),0)
    return [f, dF]
end

# function lagrangian(x, λ, μ, c::Function)
#     m = length(x)
#     F = - λ'*c(x)[1]
#     dF = - c(x)[2]*λ
#     d2F = spdiagm(-λ, 0)
#     return [F, dF, d2F]
# end

#
#
# function logBarrier(displacement::Array{Float64,1}, constraint::Array{Float64,1}, constraintIndex::Array{Int64,1})
#         d2functionValue = spzeros(length(displacement),length(displacement))
#         constraintNorm = norm(constraintIndex.*(displacement-constraint))^2
#         functionValue   = log.( constraintNorm )
#         dfunctionValue  = 1./( constraintNorm ) * constraintIndex.*(displacement-constraint)
#         println(functionValue, " c <-> dc ", norm(constraintIndex.*(displacement-constraint)))
#         return [functionValue, dfunctionValue, d2functionValue]
# end

function registerNonParametricFixedGridResolution(referenceImage, templateImage, options::regOptions;
                                     affineParameters = [1.0,0,0,0,1.0,0],
                                     measureDistance = ssdDistance,
                                     regularizerOperator=createCurvatureOperatorCentered,
                                     initialDisplacement=scaledArray(zeros(1),(1,),[],[]),
                                     constraint = 0,
                                     interpolationScheme=InterpLinearFast, gradientDescentOnly=false)#BSpline(Linear()))

  doConstraints = constraint != 0

  affineParametersInitial = affineParameters
  options.parametricOnly = false
  referenceGrid = []; deformedGrid = []; imageSize = [];
  R_level1 = restrictResolutionToLevel(referenceImage, options.levels[1])
  T_level1 = restrictResolutionToLevel(templateImage,  options.levels[1])

  regularizerMatrix = regularizerOperator(R_level1.voxelsize, getSize(R_level1))

  if initialDisplacement.data == zeros(1)
      centeredGrid = getCellCenteredGrid(R);
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

  	imageSize = getSize(R)
  	info("level ",level,": [",size(R.data)[1],"]x[",size(R.data)[2],"]")

  	# define objective function
  	Jfunc(grid;doDerivative=false,doHessian=false) =
        measureDistance(R, T, grid, doDerivative=doDerivative,
            doHessian=doHessian, options=options, centeredGrid=centeredGrid.data)[1:3] +
            options.regularizerWeight * regularizer(grid-referenceGrid.data,regularizerMatrix)
    fValues(grid) = [measureDistance(R, T, grid, doDerivative=false, doHessian=false, options=options, centeredGrid=centeredGrid.data)[1:3],
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


function registerNonParametricConstraint(referenceImage, templateImage, options::regOptions;
                                     affineParameters = [1.0,0,0,0,1.0,0],
                                     measureDistance = ssdDistance,
                                     regularizerOperator=createCurvatureOperatorCentered,
                                     initialDisplacement=scaledArray(zeros(1),(1,),[],[]),
                                     constraint = 0, subspaces = [],
                                     interpolationScheme=InterpLinearFast, gradientDescentOnly=false)#BSpline(Linear()))

  doConstraints = constraint != 0
  doASPIN = size(subspaces,1) > 1

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
    	    deformedGrid = interpolateDeformationField(initialDisplacement, referenceGrid, interpolationScheme=interpolationScheme) + referenceGrid
        end
  	else
  	      #deformation field interpolation
  	      displacementField = deformedGrid - referenceGrid
  	      referenceGrid = transformGridAffine(centeredGrid,affineParameters)
  	      deformedGrid = interpolateDeformationField(displacementField, referenceGrid, interpolationScheme=interpolationScheme) + referenceGrid
  	end

  	imageSize = getSize(R)
  	Logging.info("level ",level,": [",size(R.data)[1],"]x[",size(R.data)[2],"]")

  	# define objective function
  	regularizerMatrix = regularizerOperator(R.voxelsize, imageSize)

  	Jfunc(grid;doDerivative=false,doHessian=false) =
        measureDistance(R, T, grid, doDerivative=doDerivative,
            doHessian=doHessian, options=options, centeredGrid=centeredGrid.data)[1:3] +
            options.regularizerWeight * regularizer(grid-referenceGrid.data,regularizerMatrix)
    fValues(grid) = [measureDistance(R, T, grid, doDerivative=false, doHessian=false, options=options, centeredGrid=centeredGrid.data)[1:3],
        options.regularizerWeight , regularizer(grid-referenceGrid.data,regularizerMatrix)]
  	## gauss newton method
    if doASPIN
        deformedGrid = opt(Jfunc, deformedGrid.data, referenceGrid.data, options, printFunction = fValues, gradientDescentOnly=gradientDescentOnly, subspaces=subspaces)
    elseif doConstraints
        initialDisplacementCoarse = interpolateDeformationField(initialDisplacement, referenceGrid, interpolationScheme=interpolationScheme) + referenceGrid
        c(x) = constraint(x, initialDisplacementCoarse)
        deformedGrid = opt(Jfunc, deformedGrid.data, referenceGrid.data, options, constraint= c, printFunction = fValues)
    else
        deformedGrid = opt(Jfunc, deformedGrid.data, referenceGrid.data, options, printFunction = fValues)
    end
    deformedGrid = scaledArray(deformedGrid, size(R.data), R.voxelsize, R.shift)
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

function opt(Jfunc::Function,  # objective Function
                             y::Array{Float64,1}, yReference::Array{Float64,1}, options::regOptions;
                             constraint::Function = x -> [0, 0, 0], printFunction=x->[], gradientDescentOnly=false, subspaces=[]) #no constraint by default

    if size(subspaces,1) > 1
        return optimizeGaussNewtonASPIN(Jfunc, subspaces, y, yReference, options)
    else
        return optimizeGaussNewtonAugmentedLagrangian(Jfunc,  y, yReference, options,
          constraint=constraint, printFunction=printFunction, gradientDescentOnly=gradientDescentOnly)
    end

end
