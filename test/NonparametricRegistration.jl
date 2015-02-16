using ImageRegistration.ImageProcessing
using ImageRegistration.Examples
using ImageRegistration.Distance
using ImageRegistration.Transformation
using Base.Test

# create test images
data = zeros(120,120); data[31:90,21:60] = 1
temImg = createImage(data) 
data[41:80,41:90] = 1
refImg = createImage(data)
options = ImageRegistration.regOptions()
options.levels = [4,3,2]

# register images nonparametric
@time deformationField = registerImagesNonparametric(refImg,temImg,options)
ssdvalue = ssdDistance(refImg,temImg,getCellCenteredGrid(refImg)+deformationField)[1]
#@test_approx_eq_eps ssdvalue 124.077 1e-1
println("SSD value after nonlinear registration:")
println(ssdvalue)

# register images nonparametric
options.matrixFree = true;
@time deformationField = registerImagesNonparametric(refImg,temImg,options)
ssdvalue = ssdDistance(refImg,temImg,getCellCenteredGrid(refImg)+deformationField)[1]
@test_approx_eq_eps ssdvalue 124.077 1e-1
println("SSD (matrix-free) value after nonlinear registration:")
println(ssdvalue)

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
