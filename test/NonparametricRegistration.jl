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
@test_approx_eq_eps ssdvalue 124.077 1e-1
Logging.info("Regression test passed (nonparametric): ", ssdDistance)

# @time deformationField = registerImagesNonparametric(refImg,temImg,alpha=1,measureDistance=ngfDistance)
# ssdvalue = ssdDistance(refImg,temImg,getCellCenteredGrid(refImg)+deformationField)[1]
# @test_approx_eq_eps ssdvalue 124.077 1e-1
# println("SSD value after nonlinear NGF registration:")
# println(ssdvalue)

# visualize results
#using ImageRegistration.Visualization
#using PyPlot; pygui(true); close("all")
#figure()
#visualizeResults(refImg,temImg,deformationField=deformationField)
