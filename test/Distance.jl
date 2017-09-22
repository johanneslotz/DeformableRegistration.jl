using DeformableRegistration
using DeformableRegistration.Distance
using DeformableRegistration.ImageProcessing
using DeformableRegistration.Transformation
using DeformableRegistration.Interpolation
using DeformableRegistration.Visualization
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
options = DeformableRegistration.regOptions()
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


Logging.info("Distance: NGF parametric derivative ...")
 evaluationPoint = [1,0,0,0,1,0.0]+0.1*rand(6)
 centeredGrid = getCellCenteredGrid(refImg)
 options.edgeParameterR = 0.1
 options.edgeParameterT = 0.1
 D,dD,d2D = ngfDistance(refImg,refImg,transformGridAffine(centeredGrid,evaluationPoint),doDerivative=true,doHessian=true,options=options)
 Dfunc(p) = ngfDistance(refImg,refImg,transformGridAffine(centeredGrid,p),options=options)[1]
 errlin,errquad = checkDerivative(Dfunc,dD',evaluationPoint,doPlot=true)
@test checkErrorDecay(errquad)

#Logging.info("Checking maskedSSD cell centered derivatives...")
# measure distance and check derivative (nonparametric)
#mask = ones(refImg.data)[:]
#evaluationPoint = getCellCenteredGrid(refImg)
#evaluationPoint = evaluationPoint + 0.3 *rand(size(centeredGrid))
#options = DeformableRegistration.regOptions()
#options.doDerivative=true
#options.parametricOnly=false
#D,dD,d2D = maskedSsdDistance(refImg,refImg,evaluationPoint,mask,options)

#optNoDeriv = DeformableRegistration.regOptions()
#optNoDeriv.doDerivative=false
#optNoDeriv.parametricOnly=false
#Dfunc(x) = maskedSsdDistance(refImg,refImg,x,mask,optNoDeriv)[1]
#errlin,errquad = checkDerivative(Dfunc,dD',evaluationPoint)
#@test checkErrorDecay(errquad)

testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
refImg = loadImage(testimage)
estimatedEps, edgeImage = estimateNGFEpsilon(refImg+0.01*rand(size(refImg)),cutoffPercent=80)
@test estimatedEps >= 0
close("all")
