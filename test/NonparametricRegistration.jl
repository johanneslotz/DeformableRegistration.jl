using ImageRegistration.ImageProcessing
using ImageRegistration.Examples
using ImageRegistration.Distance
using ImageRegistration.Transformation
using Base.Test

# setup logging
Logging.configure(level=Logging.INFO)
Logging.info("Testing nonparametric registration...")

# create test images
data = zeros(120,120); data[31:90,21:60] = 1
temImg = createImage(data) 
data[41:80,41:90] = 1
refImg = createImage(data)
options = ImageRegistration.regOptions()
options.levels = [4,3,2]
options.matrixFree = true;

# register images nonparametric
deformationField = registerImagesNonparametric(refImg,temImg,options)
ssdvalue = ssdDistance(refImg,temImg,getCellCenteredGrid(refImg)+deformationField)[1]
ngfvalue = ngfDistance(refImg,temImg,getCellCenteredGrid(refImg)+deformationField)[1]
@test_approx_eq_eps ssdvalue 124.077 1e-1
@test_approx_eq_eps ngfvalue 663.522 1e-1
Logging.info("Regression test passed (nonparametric): ", ssdDistance)

ngfDistance(refImg, temImg, getCellCenteredGrid(refImg), options=options, doDerivative=true, doHessian = true)
ssdDistance(refImg, temImg, getCellCenteredGrid(refImg), options=options, doDerivative=true, doHessian = true)

# NGF, smooth data first to generate gradients
data = zeros(120,120); data[31:90,21:60] = 1
dataT = copy(data)
data = conv2([0.2 0.5 0.2
			 0.5 1.0 0.5
			 0.2 0.5 0.2],data)[1:120,1:120]
temImg = createImage(data) 
dataT[41:80,41:90] = 1
dataT = conv2([0.2 0.5 0.2
			 0.5 1.0 0.5
			 0.2 0.5 0.2],dataT)[1:120,1:120]
refImg = createImage(dataT)
options.levels = [4,3,2]
 options.edgeParameterR = 1
 options.edgeParameterT = options.edgeParameterR
 options.regularizerWeight = 0.01
 @time deformationField = registerImagesNonparametric(refImg,temImg,options,measureDistance=ngfDistance)
ssdvalue = ssdDistance(refImg,temImg,getCellCenteredGrid(refImg)+deformationField)[1]
ngfvalue = ngfDistance(refImg,temImg,getCellCenteredGrid(refImg)+deformationField)[1]

@test_approx_eq_eps ssdvalue 496.054 1e-1
@test_approx_eq_eps ngfvalue 444.766 1e-1
# println("SSD value after nonlinear NGF registration:")
# println(ssdvalue)

# visualize results
#using ImageRegistration.Visualization
# using PyPlot; pygui(true); close("all")
# figure()
# visualizeResults(refImg,temImg,deformationField=deformationField)
