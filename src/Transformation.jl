module Transformation

#using Images
using Interpolations
using DeformableRegistration.Types
# using DeformableRegistration.ImageProcessing


export getCellCenteredGrid, getNodalGrid, stg2cen, cen2stg, checkForOddNumberOfGridPoints, getCellCenteredGridRanges
export transformGridAffine, transformWorldCoordinate


function getCellCenteredGrid(I::regImage)
  return getCellCenteredGrid(I.voxelsize,I.shift,size(I.data))
end

function getCellCenteredGrid(I::scaledArray)
  return getCellCenteredGrid(I.voxelsize,I.shift,I.dimensions)
end

function getCellCenteredGrid(voxelsize::Array{Float64,1},shift::Array{Float64,1}, gridSize::Tuple{Int64,Int64})
    x = Array(shift[1]+voxelsize[1]/2:voxelsize[1]:shift[1]+voxelsize[1]*(gridSize[1]))
    y = Array(shift[2]+voxelsize[2]/2:voxelsize[2]:shift[2]+voxelsize[2]*(gridSize[2]))
    xgrid = repmat(x,1,gridSize[2])[:]
    ygrid = repmat(y',gridSize[1],1)[:]
    return scaledArray(vcat(xgrid,ygrid), gridSize, voxelsize, shift)
end

function getNodalGrid(image::regImage)
  return getNodalGrid(I.voxelsize,I.shift, size(I.data))
end

function getNodalGrid(voxelsize::Array{Float64,1},shift::Array{Float64,1}, gridSize::Tuple{Int64,Int64})
    x = Array(shift[1]:voxelsize[1]:shift[1]+voxelsize[1]*(gridSize[1]))
    y = Array(shift[2]:voxelsize[2]:shift[2]+voxelsize[2]*(gridSize[2]))
    xgrid = repmat(x,1,length(y))[:]
    ygrid = repmat(y',length(x),1)[:]
    return scaledArray(vcat(xgrid,ygrid), gridSize, voxelsize, shift)
end

function checkForOddNumberOfGridPoints(I::Array{Float64,2},grid::Array{Float64,1})
  if( 2*prod(size(I)) == size(grid,1) )
      return false
  else
      return true
  end
end

function getCellCenteredGridRanges(I::regImage, gridSize)
    x = Array(I.shift[1]+I.voxelsize[1]/2:I.voxelsize[1]:I.shift[1]+I.voxelsize[1]*(gridSize[1]))
    y = Array(I.shift[2]+I.voxelsize[2]/2:I.voxelsize[2]:I.shift[2]+I.voxelsize[2]*(gridSize[1]))
    return x,y
end

function transformWorldCoordinate(I::regImage,y::Float64,x::Float64)
    xNew = (x - I.shift[1]) / I.voxelsize[1] + .5
    yNew = (y - I.shift[2]) / I.voxelsize[2] + .5
    return yNew,xNew
end

function transformWorldCoordinate(I::regImage,grid::Array{Float64,1})
    N = size(grid,1)
    newGrid = zeros(N)
    for i=1:N
      newGrid[i] = (grid[i] - I.shift[1]) / I.voxelsize[1] + .5
      newGrid[N+i] = (grid[N+i] - I.shift[2]) / I.voxelsize[2] + .5
    end
    return newGrid
end

function transformGridAffine(grid::scaledArray,affineParameters)
    transformationMatrix = [affineParameters[1] affineParameters[2]; affineParameters[4] affineParameters[5]]
    translationParameters = [affineParameters[3]; affineParameters[6]]
    return transformGridAffine(transformationMatrix, grid, translationParameters)
end

function transformGridAffine(transformationMatrix::Matrix,gridOriginal::scaledArray,translationParameters::Array{Float64,1})
    N = prod(gridOriginal.dimensions)
    grid = zeros(2*N)
    @simd for i=1:N
        @inbounds grid[i]   = transformationMatrix[1,1]*gridOriginal.data[i]+transformationMatrix[1,2]*gridOriginal.data[N+i]+translationParameters[1]
        @inbounds grid[N+i] = transformationMatrix[2,1]*gridOriginal.data[i]+transformationMatrix[2,2]*gridOriginal.data[N+i]+translationParameters[2]
    end
    return scaledArray(grid, gridOriginal.dimensions, gridOriginal.voxelsize, gridOriginal.shift)
end

end
