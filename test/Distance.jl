using ImageRegistration.Distance
using ImageRegistration.ImageProcessing
using ImageRegistration.Transformation
using ImageRegistration.Visualization
using PyPlot
using Base.Test
using Images
# define check derivative
include("checkDerivative.jl")

# load reference image
#testimage = dirname(Base.source_path()) * "/testdata/lena.png"
#refImg = loadImage(testimage)
refImg = createImage(100*rand(50,60),spatialDomain=[0.0 50 0 60][:])
typeof(refImg)

println("Checking interpolation on deformed grid...")
# define a cell centered grid, transform it and create a template image
centeredGrid = getCellCenteredGrid(refImg)
affineParameters = [1.5,0.4,-20,0.1,1.5,-10]
deformationField = zeros(prod(size(refImg))*2)
deformationField[prod(size(refImg))+1:end] = 10*sin(0.01*centeredGrid[1:prod(size(refImg))])
transformedGrid = transformGridAffine(centeredGrid,affineParameters) + deformationField
#temImg = linearImageInterpolationAtGrid(refImg,transformedGrid)
#temImg = createImage(temImg)

println("Checking image derivatives...")
# cheack image derivative
centeredGrid = getCellCenteredGrid(refImg) + 0.1
Ifunc(x) = linearImageInterpolationAtGrid(refImg,x)[:]
transformedImage, dY_transformedImage, dX_transformedImage =
  linearImageInterpolationAtGridWithDerivative(refImg,centeredGrid)
dTransformedImage = spdiagm((dX_transformedImage[:],dY_transformedImage[:]),[0, prod(size(refImg))])
errlin,errquad = checkDerivative(Ifunc,dTransformedImage,centeredGrid)
@test checkErrorDecay(errquad)


println("Checking NGF cell centered derivatives...")
# measure distance and check derivative (nonparametric)
ε = 0.1
centeredGrid = getCellCenteredGrid(refImg) + 0.0001
D,dD,d2D,drc = ngfDistance(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true,edgeParameterR=ε,edgeParameterT=ε)
Dfunc(x) = ngfDistance(refImg,refImg,x,edgeParameterR=ε,edgeParameterT=ε)[1]
errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid, output=false)
@test checkErrorDecay(errquad)


println("Checking SSD cell centered derivatives...")
# measure distance and check derivative (nonparametric)
centeredGrid = getCellCenteredGrid(refImg) + 0.1
D,dD,d2D = ssdDistance(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true)
Dfunc(x) = ssdDistance(refImg,refImg,x)[1]
errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid, output=false)
@test checkErrorDecay(errquad)

# measure distance and check derivative (parametic)
centeredGrid = getCellCenteredGrid(refImg) + 0.001
D,dD,d2D = ssdDistance(refImg,refImg,centeredGrid,parametricOnly=true,doDerivative=true,doHessian=true)
Dfunc(p) = ssdDistance(refImg,refImg,transformGridAffine(centeredGrid,p))[1]
errlin,errquad = checkDerivative(Dfunc,dD',[1.0,0,0,0,1,0])
@test checkErrorDecay(errquad)



println("Checking maskedSSD cell centered derivatives...")
# measure distance and check derivative (nonparametric)
mask = ones(refImg.data)[:]
centeredGrid = getCellCenteredGrid(refImg) + 0.1
D,dD,d2D = maskedSsdDistance(refImg,refImg,centeredGrid,mask,doDerivative=true,doHessian=true)
Dfunc(x) = maskedSsdDistance(refImg,refImg,x,mask)[1]
errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid, output=false)
@test checkErrorDecay(errquad)


println("CURRENTLY NOT Checking SSD staggered derivatives... These derivatives are not working, this is a known bug.")
#staggeredGrid = getStaggeredGrid(refImg) + 0.1
#imageSize = [size(refImg)[1],size(refImg)[2]]
#D,dD,d2D = cen2stg( ssdDistance(refImg,refImg,stg2cen(staggeredGrid,imageSize)) ,refImg)
#Dfunc(x) = cen2stg( ssdDistance(refImg,refImg,stg2cen(x,imageSize), doDerivative=true, doHessian=true) ,refImg)[1]
#errlin,errquad = checkDerivative(Dfunc,dD',staggeredGrid)
#@test checkErrorDecay(errquad)

