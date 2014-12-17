using ImageRegistration
using ImageRegistration.ImageProcessing
using ImageRegistration.Examples
using ImageRegistration.Distance
using ImageRegistration.Transformation
using ImageRegistration.Visualization
using ImageRegistration.Regularizer
using ImageRegistration.Optimization
using Base.Test

# create test images
data = zeros(120,120); data[31:90,21:60] = 1
templateImage = createImage(data) 
data[41:80,41:90] = 1
referenceImage = createImage(data)

levels = [4,3,2]
α = 1

displacement = 0
imageSize = 0
spatialDomain = 0

# start multilevel registration
for level in levels

  R = restrictResolutionToLevel(referenceImage,level)
  spatialDomain = R.properties["spatialdomain"]
  T = restrictResolutionToLevel(templateImage,level)

  identityGrid = getCellCenteredGrid(R)

  if(level==levels[1])
    displacement = 0.0 .* getCellCenteredGrid(R)
    displacementDim = [size(R)[1],size(R)[2]]
  end

  displacement = interpolateDeformationFieldAtGrid(displacement,displacementDim,spatialDomain,identityGrid)
  displacementDim = [size(R)[1],size(R)[2]]

  # update imageSize
  imageSize = [size(R)[1],size(R)[2]]

  # define objective function
  regulizerMatrix = createDiffusiveOperatorCentered(R.properties["pixelspacing"],imageSize)

  JntObjtvFctn(displacement) =
    ssdDistance(R,T,displacement+identityGrid) +
    α * regularizer(displacement,regulizerMatrix)
  JntObjtvFctnWDerivative(displacement) =
    ssdDistance(R,T,displacement+identityGrid,doDerivative=true,doHessian=true) +
    α * regularizer(displacement,regulizerMatrix,doDerivative=true,doHessian=true)

  # gauss newton method
  displacement = optimizeGaussNewton(JntObjtvFctn,JntObjtvFctnWDerivative,displacement,output=false)

end

# upsample deformation to reference image
displacement = interpolateDeformationFieldAtGrid(displacement,imageSize,spatialDomain,getCellCenteredGrid(referenceImage))

# add this visualization if needed
#using PyPlot; pygui(true); close("all")
#figure()
#visualizeResults(referenceImage,templateImage,deformationField=displacement,numberOfGridLinesX=20,numberOfGridLinesY=20)

@test_approx_eq_eps sum(displacement) -189175.77072168788 1e-1


