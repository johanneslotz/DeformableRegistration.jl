using DeformableRegistration.Transformation
 using DeformableRegistration.ImageProcessing
 using Test

@testset "Transformation" begin
@testset "grid creation" begin
    # test getCellCenteredGrid
    spatialDomain = [-0.1,5.2,-1.1,3.3]
    gridSize=(5,8)
    voxelsize = [(spatialDomain[2]-spatialDomain[1])/gridSize[1],(spatialDomain[4]-spatialDomain[3])/gridSize[2]]
    shift = spatialDomain[1:2:3]
    centeredGrid = getCellCenteredGrid(voxelsize, shift ,gridSize)

    @test length(centeredGrid.data) == prod(gridSize)*2

    @test centeredGrid.data[1] == (shift[1] + voxelsize[1]/2)
    @test centeredGrid.data[prod(gridSize)+1] == (shift[2]+voxelsize[2]/2)

    nodalGrid = getNodalGrid(voxelsize,shift,gridSize)
    @test length(nodalGrid.data) == (gridSize[1]+1)*(gridSize[2]+1)*2


end
end


# @testset "stg <-> cen conversion" begin
#     gridSize=(5,8)
#     voxelsize = [1.06,0.55]
#     shift = [-0.1,-1.1]
#     centeredGrid = getCellCenteredGrid(voxelsize, shift ,gridSize)
#     staggeredGrid = getStaggeredGrid(voxelsize,shift,gridSize)
#     centeredGridFromStaggered = stg2cen(staggeredGrid,gridSize)
#     @test_approx_eq centeredGrid centeredGridFromStaggered
#     # test cen2stg
#     staggeredGridFromCenteredGrid = cen2stg(centeredGrid,gridSize)
# end

# visualize cell centered grid
#using DeformableRegistration.Visualization
#using PyPlot
#pygui(true); close("all"); PyPlot.svg(true)
#imgdata = rand(gridSize[1],gridSize[2])
#img = createImage(imgdata,spatialDomain=spatialDomain)
#centeredGrid = getCellCenteredGrid(img)
#figure()
#showImage(img)
#showGrid(spatialDomain,gridSize,centeredGrid,showPoints=true,showIndices=true)
