using ImageRegistration.Transformation
using ImageRegistration.ImageProcessing
using ImageRegistration.Interpolation
using Base.Test
using Logging
using Grid

# setup logging
Logging.configure(level = Logging.INFO) # set to DEBUG to see error tables

# define check derivative
include("helpers/checkDerivative.jl")

# load reference image
testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"

# check interpolation againt Grid.jl
Logging.info("Interpolation: Checking own linear interpolation against the linear interpolation of the Grid.jl package...")
img = loadImage(testimage)
centeredGrid = getCellCenteredGrid(img)
affineParameters = [0.5,0.4,50,0,0.5,100]
deformationField = zeros(prod(size(img))*2)
deformationField[prod(size(img))+1:end] = 10*sin(0.01*centeredGrid[1:prod(size(img))])
transformedGrid = transformGridAffine(centeredGrid,affineParameters) + deformationField
transformedImageOnly = interpolateImage(img,transformedGrid)
transformedImage,dY_transformedImage,dX_transformedImage = interpolateImage(img,transformedGrid,doDerivative=true)
transformedImageRef,dY_transformedImageRef,dX_transformedImageRef = interpolateImage(img,transformedGrid,doDerivative=true,InterpFunction=InterpLinear)
@test_approx_eq transformedImageOnly transformedImageRef
@test_approx_eq transformedImage transformedImageRef
Logging.info("Interpolation: Transformed images are equal ✔")
@test_approx_eq dY_transformedImage dY_transformedImageRef
@test_approx_eq dX_transformedImage dX_transformedImageRef
Logging.info("Interpolation: Derivatives of the transformed images are equal ✔")

# check image derivative
Logging.info("Interpolation: Checking image derivatives...")
centeredGrid = getCellCenteredGrid(img) + 0.1
Ifunc(x) = interpolateImage(img,x)[:]
transformedImage, dY_transformedImage, dX_transformedImage = interpolateImage(img,centeredGrid,doDerivative=true)
dTransformedImage = spdiagm((dX_transformedImage[:],dY_transformedImage[:]),[0, prod(size(img))])
errlin,errquad = checkDerivative(Ifunc,dTransformedImage,centeredGrid)
@test checkErrorDecay(errquad)
Logging.info("Interpolation: Sufficient error decay of taylor approximation ✔")

# check timing
Logging.info("Interpolation: Comparing timings of linear interpolation...")
Logging.info("Interpolation: Image and grid size: 1024x1024")
img = createImage(rand(1024,1024))
centeredGrid = getCellCenteredGrid(img)
affineParameters = [0.5,0.4,50,0,0.5,100]
deformationField = zeros(prod(size(img))*2)
deformationField[prod(size(img))+1:end] = 10*sin(0.01*centeredGrid[1:prod(size(img))])
transformedGrid = transformGridAffine(centeredGrid,affineParameters) + deformationField
# InterpLinearFast
tic();
for i=1:5
  interpolateImage(img,transformedGrid)
end
timing = toq()/5;
Logging.info("Interpolation: interpolateImage (InterpLinearFast) took ",timing," seconds.")
# InterpLinear (Grid.jl)
tic();
for i=1:5
  interpolateImage(img,transformedGrid,InterpFunction=InterpLinear)
end
timingGrid = toq()/5;
Logging.info("Interpolation: interpolateImage (InterpLinear, Grid.jl) took ",timingGrid," seconds.")
# InterpLinearFast with derivative
tic();
for i=1:5
  interpolateImage(img,transformedGrid,doDerivative=true)
end
timing = toq()/5;
Logging.info("Interpolation: interpolateImage with derivative (InterpLinearFast) took ",timing," seconds.")
# InterpLinear (Grid.jl) with derivative
tic();
for i=1:5
  interpolateImage(img,transformedGrid,InterpFunction=InterpLinear,doDerivative=true)
end
timingGrid = toq()/5;
Logging.info("Interpolation: interpolateImage with derivative (InterpLinear, Grid.jl) took ",timingGrid," seconds.")
@test timingGrid>10*timing
Logging.info("Interpolation: Own linear interpolation with derivative is at least ten times faster ✔")

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
