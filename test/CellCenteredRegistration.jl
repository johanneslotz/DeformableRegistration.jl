using ImageRegistration
using ImageRegistration.ImageProcessing
using ImageRegistration.Examples
using ImageRegistration.Distance
using ImageRegistration.Transformation
using ImageRegistration.Interpolation
using ImageRegistration.Visualization
using ImageRegistration.Regularizer
using ImageRegistration.Optimization
using Base.Test

# create test images
data = zeros(120,120); data[31:90,21:60] = 1
templateImage = createImage(data) 
data[41:80,41:90] = 1
referenceImage = createImage(data)

displacement = 0
imageSize = 0
spatialDomain = 0

options = regOptions()
options.regularizerWeight = 1
options.levels = [4,3,2]

# start multilevel registration
for level in options.levels

  R = restrictResolutionToLevel(referenceImage,level)
  spatialDomain = getSpatialDomain(R)
  T = restrictResolutionToLevel(templateImage,level)

  identityGrid = getCellCenteredGrid(R)
  displacementInitial = 0.0 .* getCellCenteredGrid(R)

  if(level==options.levels[1])
    displacement = 0.0 .* getCellCenteredGrid(R)
    displacementDim = getSize(R)
  end

  displacement = interpolateDeformationField(displacement,displacementDim,spatialDomain,identityGrid)
  displacementDim = getSize(R)

  # update imageSize
  imageSize = getSize(R)

  # define objective function
  regulizerMatrix = createDiffusiveOperatorCentered(getPixelSpacing(R),imageSize)
  JntObjtvFctn(displacement;doDerivative=false,doHessian=false) =
    ssdDistance(R,T,displacement+identityGrid,doDerivative=doDerivative,doHessian=doHessian,options=options) +
    options.regularizerWeight * regularizer(displacement,regulizerMatrix)

  # gauss newton method
  displacement = optimizeGaussNewton(JntObjtvFctn,displacement,displacementInitial, options)

end

# upsample deformation to reference image
displacement = interpolateDeformationField(displacement,imageSize,spatialDomain,getCellCenteredGrid(referenceImage))

# add this visualization if needed
using PyPlot; pygui(true); close("all")
figure()
visualizeResults(referenceImage,templateImage,deformationField=displacement,numberOfGridLinesX=20,numberOfGridLinesY=20)

@test_approx_eq_eps sum(displacement) -189332.07 1e-2
