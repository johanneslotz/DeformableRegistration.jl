using ImageRegistration
using ImageRegistration.Distance
using ImageRegistration.ImageProcessing
using ImageRegistration.Transformation
using ImageRegistration.Interpolation
using ImageRegistration.Visualization
using Base.Test
using PyPlot
using Images
using Logging

# setup logging
Logging.configure(level = Logging.INFO) # set to DEBUG to see error tables

# define check derivative
include("helpers/checkDerivative.jl")

# create reference image
refImg = createImage(100*rand(50,60))

# println("Checking NGF cell centered derivatives...")
# # measure distance and check derivative (nonparametric)
# options = ImageRegistration.getDefaultOptions()
# options["edgeParameterR"] = 0.1
# options["edgeParameterT"] = 0.1
# options["centeredGrid"] = getCellCenteredGrid(refImg) + 0.0001
# options["parametricOnly"] = false

# D,dD,d2D,drc = ngfDistance(refImg,refImg,options, options["centeredGrid"] )
# Dfunc(x) = ngfDistance(refImg,refImg,options,x)[1]
# errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid, output=true)
# @test checkErrorDecay(errquad)

# measure distance and check derivative (parametic)
#centeredGrid = getCellCenteredGrid(refImg) + 0.001
#options["parametricOnly"] = true

#D,dD,d2D = ngfDistance(refImg,refImg,options,centeredGrid)
#Dfunc(p) = ngfDistance(refImg,refImg,options,transformGridAffine(centeredGrid,p))[1]
#errlin,errquad = checkDerivative(Dfunc,dD',[1.0,0,0,0,1,0])
#@test checkErrorDecay(errquad)

Logging.info("Distance: Checking SSD...")
# measure distance and check derivative (nonparametric)
centeredGrid = getCellCenteredGrid(refImg)
centeredGrid = centeredGrid + 0.3 *rand(size(centeredGrid))
D,dD,d2D = ssdDistance(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true)
Dfunc(x) = ssdDistance(refImg,refImg,x)[1]
errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid,doPlot=true)
@test checkErrorDecay(errquad)

Logging.info("Distance: SSD nonparametric derivative ✔")
options = regOptions()
options.matrixFree = true;
D,dD,d2D = ssdDistance(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true,options=options)
Dfunc(x) = ssdDistance(refImg,refImg,x)[1]
errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid)
@test checkErrorDecay(errquad)

Logging.info("Distance: SSD nonparametric matrixfree derivative ✔")
# measure distance and check derivative (parametic)
centeredGrid = getCellCenteredGrid(refImg)
options = ImageRegistration.regOptions()
options.parametricOnly=true
evaluationPoint = [1,0,0,0,1,0]+0.1*rand(6)
D,dD,d2D = ssdDistance(refImg,refImg,transformGridAffine(centeredGrid,evaluationPoint),doDerivative=true,doHessian=true,options=options)
Dfunc(p) = ssdDistance(refImg,refImg,transformGridAffine(centeredGrid,p),options=options)[1]
errlin,errquad = checkDerivative(Dfunc,dD',evaluationPoint)
@test checkErrorDecay(errquad)

Logging.info("Distance: SSD parametric derivative ✔")
options.matrixFree = true;
D,dD,d2D = ssdDistance(refImg,refImg,transformGridAffine(centeredGrid,evaluationPoint),doDerivative=true,doHessian=true,options=options)
Dfunc(p) = ssdDistance(refImg,refImg,transformGridAffine(centeredGrid,p),options=options)[1]
errlin,errquad = checkDerivative(Dfunc,dD',evaluationPoint)
@test checkErrorDecay(errquad)
Logging.info("Distance: SSD parametric matrixfree derivative ✔")

Logging.info("Distance: Checking NGF...")
# measure distance and check derivative (nonparametric)
 centeredGrid = getCellCenteredGrid(refImg)
 centeredGrid = centeredGrid + 0.3 *rand(size(centeredGrid))
 D,dD,d2D = ngfDistance(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true)
 Dfunc(x) = ngfDistance(refImg,refImg,x)[1]
 errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid,doPlot=true)
@test checkErrorDecay(errquad)


#Logging.info("Checking maskedSSD cell centered derivatives...")
# measure distance and check derivative (nonparametric)
#mask = ones(refImg.data)[:]
#evaluationPoint = getCellCenteredGrid(refImg)
#evaluationPoint = evaluationPoint + 0.3 *rand(size(centeredGrid))
#options = ImageRegistration.regOptions()
#options.doDerivative=true
#options.parametricOnly=false
#D,dD,d2D = maskedSsdDistance(refImg,refImg,evaluationPoint,mask,options)

#optNoDeriv = ImageRegistration.regOptions()
#optNoDeriv.doDerivative=false
#optNoDeriv.parametricOnly=false
#Dfunc(x) = maskedSsdDistance(refImg,refImg,x,mask,optNoDeriv)[1]
#errlin,errquad = checkDerivative(Dfunc,dD',evaluationPoint)
#@test checkErrorDecay(errquad)

Logging.warn("CURRENTLY NOT Checking SSD staggered derivatives... These derivatives are not working, this is a known bug.")
#staggeredGrid = getStaggeredGrid(refImg) + 0.1
#imageSize = [size(refImg)[1],size(refImg)[2]]
#D,dD,d2D = cen2stg( ssdDistance(refImg,refImg,stg2cen(staggeredGrid,imageSize)) ,refImg)
#Dfunc(x) = cen2stg( ssdDistance(refImg,refImg,stg2cen(x,imageSize), doDerivative=true, doHessian=true) ,refImg)[1]
#errlin,errquad = checkDerivative(Dfunc,dD',staggeredGrid)
#@test checkErrorDecay(errquad)

testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
refImg = loadImage(testimage)
estimatedEps, edgeImage = estimateNGFEpsilon(refImg+0.01*rand(size(refImg)),cutoffPercent=80)
@test estimatedEps >= 0
close("all")
