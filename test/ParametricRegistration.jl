using ImageRegistration.ImageProcessing
using ImageRegistration.Distance
using ImageRegistration.Transformation
using ImageRegistration.Interpolation
using ImageRegistration.Examples
using ImageRegistration.Distance
using Base.Test
using Logging

# setup logging
Logging.configure(level=Logging.INFO)
Logging.info("Testing parametric registration...")

# load reference image
testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
refImg = loadImage(testimage)

# define a cell centered grid, transform it and create a template image
centeredGrid = getCellCenteredGrid(refImg)
affineParametersInitial = [1.5,0.4,-200,0.1,1.5,-100]
transformedGrid = transformGridAffine(centeredGrid,affineParametersInitial)
temImg = interpolateImage(refImg,transformedGrid)
temImg = createImage(temImg)

# set options
options = ImageRegistration.regOptions()
options.parametricOnly = true;
options.matrixFree = true;

options.useEdgeParameterInNumerator
# register images
affineParameters = registerImagesParametric(refImg,temImg, options)
finalDistanceSSD = ssdDistance(refImg,temImg,transformGridAffine(centeredGrid,affineParameters))
#finalDistanceNGF = ngfDistance(refImg,temImg,transformGridAffine(centeredGrid,affineParameters),options)
#@test_approx_eq_eps finalDistanceNGF[1] 71513.2 3.0
@test_approx_eq_eps finalDistanceSSD[1] 499.85 5e-1
Logging.info(string("Regression test passed (parametric): ", finalDistanceSSD)

## register images with ngf matrix based
#options.edgeParameterR =  1e-2
# options.edgeParameterT =  options.edgeParameterR
# options.levels= [5,4,3,2,1]
#@time affineParameters = ImageRegistration.Examples.registerImagesParametric(refImg,temImg,options,measureDistance=ngfDistance)
#finalDistanceSSD = ssdDistance(refImg,temImg,transformGridAffine(centeredGrid,affineParameters))
#finalDistanceNGF = ngfDistance(refImg,temImg,transformGridAffine(centeredGrid,affineParameters), doDerivative=false)
#finalDistanceNGF
#affineParameters
# @test_approx_eq_eps finalDistanceNGF[1] 71513.2 1e-1
# @test_approx_eq_eps finalDistanceSSD[1] 499.85 1e-1
# Logging.info("regression test passed (parametric): ", ngfDistance)

# visualize
using ImageRegistration.Visualization
 using PyPlot; pygui(true); close("all")
 figure()
 visualizeResults(refImg,temImg,affineParameters=affineParameters)
