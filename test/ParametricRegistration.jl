using ImageRegistration.ImageProcessing
using ImageRegistration.Distance
using ImageRegistration.Transformation
using ImageRegistration.Examples
using ImageRegistration.Visualization
using Base.Test

println("testing parametric registration...")

# load reference image
testimage = dirname(Base.source_path()) * "/testdata/lena.png"
refImg = loadImage(testimage)

# define a cell centered grid, transform it and create a template image
centeredGrid = getCellCenteredGrid(refImg)
affineParameters = [1.5,0.4,-200,0.1,1.5,-100]
#deformationField = zeros(prod(size(refImg))*2)
#deformationField[prod(size(refImg))+1:end] = 10*sin(0.05*centeredGrid[1:prod(size(refImg))])
transformedGrid = transformGridAffine(centeredGrid,affineParameters) #+ deformationField
temImg = linearImageInterpolationAtGrid(refImg,transformedGrid)
temImg = createImage(temImg)

# register images
@time affineParameters = registerImagesParametric(refImg,temImg, measureDistance=ssdDistance)
@test_approx_eq_eps affineParameters [0.67762891033478,-0.18051497507358655,117.88728616297367,-0.04520880359582362,0.6765803426265078,59.969880090069566] 1e-1

# register images with ssd matrix free
@time affineParameters = registerImagesParametric(refImg,temImg,measureDistance=ssdDistanceMatrixFree)
@test_approx_eq_eps affineParameters [0.67762891033478,-0.18051497507358655,117.88728616297367,-0.04520880359582362,0.6765803426265078,59.969880090069566] 1e-1

# visualize
#using PyPlot; pygui(true); close("all")
#figure()
#visualizeResults(refImg,temImg,affineParameters=affineParameters)
