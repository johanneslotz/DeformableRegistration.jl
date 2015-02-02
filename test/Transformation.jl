using ImageRegistration.Transformation
using ImageRegistration.ImageProcessing
using Base.Test

### GRID

# test getCellCenteredGrid
spatialDomain = [-0.1,5.2,-1.1,3.3]
gridSize=[5,8]
centeredGrid = getCellCenteredGrid(spatialDomain,gridSize)
pixelSpacing = [(spatialDomain[2]-spatialDomain[1])/gridSize[1],(spatialDomain[4]-spatialDomain[3])/gridSize[2]]
# test first point
# x-direction
@test centeredGrid[1] == (spatialDomain[3] + pixelSpacing[2]/2)
# y-direction
@test centeredGrid[prod(gridSize)+1] == (spatialDomain[1]+pixelSpacing[1]/2)

# test staggeredGrid and stg2cen
staggeredGrid = getStaggeredGrid(spatialDomain,gridSize)
centeredGridFromStaggered = stg2cen(staggeredGrid,gridSize)
@test_approx_eq centeredGrid centeredGridFromStaggered

# test cen2stg
staggeredGridFromCenteredGrid = cen2stg(centeredGrid,gridSize)

# visualize cell centered grid
#using ImageRegistration.Visualization
#using PyPlot
#pygui(true); close("all"); PyPlot.svg(true)
#imgdata = rand(gridSize[1],gridSize[2])
#img = createImage(imgdata,spatialDomain=spatialDomain)
#centeredGrid = getCellCenteredGrid(img)
#figure()
#showImage(img)
#showGrid(spatialDomain,gridSize,centeredGrid,showPoints=true,showIndices=true)

### INTERPOLATION

# test linearImageInterpolationAtGrid, linearImageInterpolationAtGridWithDerivative
testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
img = loadImage(testimage)
centeredGrid = getCellCenteredGrid(img)
affineParameters = [0.5,0.4,50,0,0.5,100]
deformationField = zeros(prod(size(img))*2)
deformationField[prod(size(img))+1:end] = 10*sin(0.01*centeredGrid[1:prod(size(img))])
transformedGrid = transformGridAffine(centeredGrid,affineParameters) + deformationField
transformedImageOnly = linearImageInterpolationAtGrid(img,transformedGrid)
transformedImage,dY_transformedImage,dX_transformedImage = linearImageInterpolationAtGridWithDerivative(img,transformedGrid)
transformedImageRef,dY_transformedImageRef,dX_transformedImageRef = interpolateImageAtGridWithDerivative(img,transformedGrid)
@test_approx_eq transformedImageOnly transformedImageRef
@test_approx_eq transformedImage transformedImageRef
@test_approx_eq dY_transformedImage dY_transformedImageRef
@test_approx_eq dX_transformedImage dX_transformedImageRef

# visualisation
#using PyPlot
#using ImageRegistration.Visualization
#pygui(true); close("all")
#figure(figsize=(10,10));
#subplot(2,2,1)
#showImage(img,name="image and transformed grid")
#showGrid(img.properties["spatialdomain"],size(img),centeredGrid,affineParameters=affineParameters,deformationField=deformationField,gridColor="red")
#subplot(2,2,2)
#showImage(transformedImage,name="image at grid points")
#subplot(2,2,3)
#showImage(dY_transformedImage,name="dy")
#subplot(2,2,4)
#showImage(dX_transformedImage,name="dx")

# check timing
#img = createImage(rand(2000,2000))
#centeredGrid = getCellCenteredGrid(img)
#affineParameters = [0.5,0.4,50,0,0.5,100]
#deformationField = zeros(prod(size(img))*2)
#deformationField[prod(size(img))+1:end] = 10*sin(0.01*centeredGrid[1:prod(size(img))])
#transformedGrid = transformGridAffine(centeredGrid,affineParameters) + deformationField
#tic();
#transformedImage,dY_transformedImage,dX_transformedImage = linearImageInterpolationAtGridWithDerivative(img,transformedGrid)
#timing = toc();
#print("... linearImageInterpolationAtGridWithDerivative took ",timing," seconds (1.0s needed to pass)\n")
#@test timing < 1.0
