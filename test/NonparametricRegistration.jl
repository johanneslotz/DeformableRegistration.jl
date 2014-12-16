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

print(refImg["spatialdomain"])
print(refImg["pixelspacing"])
print(size(refImg),"\n")

# register images nonparametric
deformationField = registerImagesNonparametric(refImg,temImg)
#@test_approx_eq_eps ssdDistance(refImg,temImg,getCellCenteredGrid(refImg)+deformationField)[1] 117.80394650617919 1e-0

print(refImg["spatialdomain"])
print(refImg["pixelspacing"])
print(size(refImg))

# visualize results
using PyPlot; pygui(true); close("all")
figure()
visualizeResults(refImg,temImg,deformationField=deformationField)
