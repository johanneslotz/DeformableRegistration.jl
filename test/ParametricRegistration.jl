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
Logging.info("Regression test passed (parametric): ")

# register images with ngf matrix based
#data = zeros(120,120); data[31:90,21:60] = 1
#for i in [1:10]
#	data = conv2([0.2 0.5 0.2
#				 0.5 1.0 0.5
#				 0.2 0.5 0.2],data)[1:120,1:120]
#end
#temImg = createImage(data) 
#dataT = zeros(120,120); dataT[31:90,31:70] = 1
#for i in [1:10]
#	dataT = conv2([0.2 0.5 0.2
#			 0.5 1.0 0.5
#			 0.2 0.5 0.2],dataT)[1:120,1:120]
#end
#refImg = createImage(dataT)
options.levels = [7,6,5,4]
 options.edgeParameterR = 0.1
 options.edgeParameterT = options.edgeParameterR
 options.maxIterGaussNewton = 100
@time affineParameters = ImageRegistration.Examples.registerImagesParametric(refImg,temImg,options,measureDistance=ngfDistance)
finalDistanceSSD = ssdDistance(refImg,temImg,transformGridAffine(centeredGrid,affineParameters))
finalDistanceNGF = ngfDistance(refImg,temImg,transformGridAffine(centeredGrid,affineParameters), doDerivative=false)
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
