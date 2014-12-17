using ImageRegistration.ImageProcessing
using ImageRegistration.Examples
using ImageRegistration.Distance
using ImageRegistration.Transformation
using ImageRegistration.Visualization
using Base.Test

# create test images
data = zeros(120,120); data[31:90,21:60] = 1
temImg = createImage(data) 
data[41:80,41:90] = 1
refImg = createImage(data)

# register images nonparametric
deformationField = registerImagesNonparametric(refImg,temImg,alpha=0.5)
#@test_approx_eq_eps
ssdvalue = ssdDistance(refImg,temImg,getCellCenteredGrid(refImg)+deformationField)[1]

println("The recursion test is temporarily disabled, SSD value after nonlinear registration:")
println(ssdvalue)

# visualize results
# using PyPlot; pygui(true); close("all")
# figure()
# visualizeResults(refImg,temImg,deformationField=deformationField)
