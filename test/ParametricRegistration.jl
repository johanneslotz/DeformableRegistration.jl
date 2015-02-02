using ImageRegistration.ImageProcessing
 using ImageRegistration.Distance
 using ImageRegistration.Transformation
 using ImageRegistration.Examples
 #using ImageRegistration.Visualization
 using ImageRegistration.Distance
 using Base.Test
 import Logging

Logging.configure(level=Logging.INFO)

println("testing parametric registration...")

# load reference image
testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
refImg = loadImage(testimage)

# define a cell centered grid, transform it and create a template image
centeredGrid = getCellCenteredGrid(refImg)
affineParameters = [1.5,0.4,-200,0.1,1.5,-100]
#deformationField = zeros(prod(size(refImg))*2)
#deformationField[prod(size(refImg))+1:end] = 10*sin(0.05*centeredGrid[1:prod(size(refImg))])
transformedGrid = transformGridAffine(centeredGrid,affineParameters) #+ deformationField
temImg = linearImageInterpolationAtGrid(refImg,transformedGrid)
temImg = createImage(temImg)

options = ImageRegistration.getDefaultOptions()

# register images
@time affineParameters = registerImagesParametric(refImg,temImg, ImageRegistration.getDefaultOptions(), measureDistance=ssdDistance)

finalDistanceSSD = ssdDistance(refImg,temImg,transformGridAffine(centeredGrid,affineParameters),doDerivative=true,parametricOnly =true)
finalDistanceNGF = ngfDistance(refImg,temImg,options,transformGridAffine(centeredGrid,affineParameters))
@test_approx_eq_eps finalDistanceNGF[1] 71513.2 3.0
@test_approx_eq_eps finalDistanceSSD[1] 499.85 1e-1
Logging.info("regression test passed (parametric): ", ssdDistance)

# register images with ssd matrix free
@time affineParameters = registerImagesParametric(refImg,temImg,ImageRegistration.getDefaultOptions(),measureDistance=ssdDistanceMatrixFree)
finalDistanceSSD = ssdDistance(refImg,temImg,transformGridAffine(centeredGrid,affineParameters))
finalDistanceNGF = ngfDistance(refImg,temImg,options,transformGridAffine(centeredGrid,affineParameters))

@test_approx_eq_eps finalDistanceNGF[1] 71513.2 3.0
@test_approx_eq_eps finalDistanceSSD[1] 499.85 1e-1
Logging.info("regression test passed (parametric): ", ssdDistanceMatrixFree)



# register images with ngf matrix based
# options = ImageRegistration.getDefaultOptions()
#  options["edgeParameterR"] = 0.001#*Distance.estimateNGFEpsilon(refImg)
#  options["edgeParameterT"] = 0.001#*Distance.estimateNGFEpsilon(temImg)

# @time affineParameters = ImageRegistration.Examples.registerImagesParametric(refImg,temImg,options,measureDistance=ngfDistance,levels=[5,4,3,2])
# finalDistanceSSD = ssdDistance(refImg,temImg,transformGridAffine(centeredGrid,affineParameters))
# finalDistanceNGF = ngfDistance(refImg,temImg,transformGridAffine(centeredGrid,affineParameters),options = options)

# @test_approx_eq_eps finalDistanceNGF[1] 71513.2 1e-1
# @test_approx_eq_eps finalDistanceSSD[1] 499.85 1e-1
# Logging.info("regression test passed (parametric): ", ngfDistance)

# visualize

# using ImageRegistration.Visualization
#  using PyPlot; pygui(true); close("all")
#  figure()
#  visualizeResults(refImg,temImg,affineParameters=affineParameters)


