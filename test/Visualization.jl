using DeformableRegistration: ImageProcessing, Transformation, Visualization
 using PyPlot
 using Base.Test
#  using Logging

# setup logging

@testset "Visualization" begin
    testimage = dirname(Base.source_path()) * "/testdata/testimage.png"
    img = loadImage(testimage)
    figure();
    showImage(img,name="/testdata/testimage.png")
    @test typeof(gcf()) == Figure
end

close("all")
#
# # test showGrid
# gridSize = (256,256)
#  voxelsize = [0.1, 0.1]
#  shift = [5,5.0]
#  centeredGrid = getCellCenteredGrid(voxelsize, shift, gridSize)
#  affineParameters = [0.5,0.1,0,0,0.5,0.5]
#  deformationField = zeros(prod(gridSize)*2)
#  deformationField[prod(gridSize)+1:end] = 0.05*sin(5*centeredGrid[1:prod(gridSize)])
#  figure();
#  showGrid(spatialDomain,gridSize,centeredGrid,
#          affineParameters=affineParameters,
#          deformationField=deformationField)
# # check if figure was created
# @assert typeof(gcf()) == Figure
# Logging.info("Visualization: showGrid: Figure was created ✔")
# # check y and x limits
# @test ylim() == (2.1,-1.1)
# @test xlim() == (-1.0,2.21)
# Logging.info("Visualization: showGrid: x,y limits are alright ✔")
