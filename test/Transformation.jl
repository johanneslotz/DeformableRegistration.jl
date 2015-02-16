using ImageRegistration.Transformation
using ImageRegistration.ImageProcessing
using Base.Test
using Logging

# setup logging
Logging.configure(level = Logging.INFO) # set to DEBUG to see error tables
Logging.info("Transformation: Testing grid points creation, transformation and grid changes...")

# test getCellCenteredGrid
spatialDomain = [-0.1,5.2,-1.1,3.3]
gridSize=[5,8]
centeredGrid = getCellCenteredGrid(spatialDomain,gridSize)
pixelSpacing = [(spatialDomain[2]-spatialDomain[1])/gridSize[1],(spatialDomain[4]-spatialDomain[3])/gridSize[2]]
# test first point
# x-direction
@test centeredGrid[1] == (spatialDomain[3] + pixelSpacing[2]/2)
# y-direction
@test centeredGrid[prod(gridSize)+1] == (spatialDomain[1]+pixelSpacing[1]/2)
Logging.info("Transformation: getCellCenteredGrid ✔")

# test staggeredGrid and stg2cen
staggeredGrid = getStaggeredGrid(spatialDomain,gridSize)
centeredGridFromStaggered = stg2cen(staggeredGrid,gridSize)
@test_approx_eq centeredGrid centeredGridFromStaggered
Logging.info("Transformation: getStaggeredGrid, stg2cen ✔")

# test cen2stg
staggeredGridFromCenteredGrid = cen2stg(centeredGrid,gridSize)
Logging.info("Transformation: cen2stg ✔")

# visualize cell centered grid
#using ImageRegistration.Visualization
#using PyPlot
#pygui(true); close("all"); PyPlot.svg(true)
#imgdata = rand(gridSize[1],gridSize[2])
#img = createImage(imgdata,spatialDomain=spatialDomain)
#centeredGrid = getCellCenteredGrid(img)
#figure()
#showImage(img)
#showGrid(spatialDomain,gridSize,centeredGrid,showPoints=true,showIndices=true)

