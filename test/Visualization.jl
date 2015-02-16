using ImageRegistration.ImageProcessing
using ImageRegistration.Transformation
using ImageRegistration.Visualization
using PyPlot
using Base.Test

# setup logging
Logging.configure(level = Logging.INFO) # set to DEBUG to see error tables
Logging.info("Visualization: Testing visualization...")

# test showImage
testimage = dirname(Base.source_path()) * "/testdata/testimage.png"
img = loadImage(testimage,spatialDomain=[-1.1,2.1,-1.0,2.21])
figure();
showImage(img,name="/testdata/testimage.png")
# check if figure was created
@assert typeof(gcf()) == Figure
Logging.info("Visualization: showImage: Figure was created ✔")
# check y and x limits
@test ylim() == (2.1,-1.1)
@test xlim() == (-1.0,2.21)
Logging.info("Visualization: showImage: x,y limits are alright ✔")

# test showGrid
gridSize = [256,256]
spatialDomain = [-1.1,2.1,-1.0,2.21]
centeredGrid = getCellCenteredGrid(spatialDomain,gridSize)
affineParameters = [0.5,0.1,0,0,0.5,0.5]
deformationField = zeros(prod(gridSize)*2)
deformationField[prod(gridSize)+1:end] = 0.05*sin(5*centeredGrid[1:prod(gridSize)])
figure();
showGrid(spatialDomain,gridSize,centeredGrid,
         affineParameters=affineParameters,
         deformationField=deformationField)
# check if figure was created
@assert typeof(gcf()) == Figure
Logging.info("Visualization: showGrid: Figure was created ✔")
# check y and x limits
@test ylim() == (2.1,-1.1)
@test xlim() == (-1.0,2.21)
Logging.info("Visualization: showGrid: x,y limits are alright ✔")
