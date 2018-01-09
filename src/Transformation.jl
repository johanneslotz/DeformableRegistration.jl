module Transformation

using Images
using Interpolations

using DeformableRegistration.ImageProcessing

export getCellCenteredGrid, getNodalGrid, stg2cen, cen2stg, checkForOddNumberOfGridPoints, getCellCenteredGridRanges
export transformGridAffine, transformWorldCoordinate, scaledArray, + , - , *


struct scaledArray
           data::Array{Float64, 1}
           dimensions::Tuple{Vararg{Int64}}
           voxelsize::Array{Float64, 1}
           shift::Array{Float64, 1}
end

import Base.+
import Base.-
import Base.*
import Base.size

function +(a::scaledArray, b::scaledArray)
    @assert a.voxelsize == b.voxelsize "Array's world matrices must match"
    @assert a.shift == b.shift "Array's world matrices must match"
    return(scaledArray(a.data+b.data, a.dimensions, a.voxelsize, a.shift))
end

function +(a::scaledArray, b::Array{Float64,1}) # for derivative check
    return(scaledArray(a.data+b, a.dimensions, a.voxelsize, a.shift))
end

function -(a::scaledArray, b::scaledArray)
    return(a + (-1 * b))
end

function *(s::Number, a::scaledArray)
    b = deepcopy(a)
    b.data[:] = s * b.data[:]
    return b
end

function size(a::scaledArray)
    return size(a.data)
end

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

# # staggered grid to Cellcentered grid
# function stg2cen(staggeredGrid::Array{Float64,1},gridSize::Tuple{Int,Int})
#     N = prod(gridSize)
#     M = gridSize[1]*(gridSize[2]+1)
#     centeredGrid = zeros(2*N)
#     for j=1:gridSize[2]
#         for i=1:gridSize[1]
#             centeredGrid[  (j-1)*gridSize[1]+i] = (staggeredGrid[(j-1)* (gridSize[1]+1)+i]      + staggeredGrid[j*gridSize[1]+i]            ) / 2
#             centeredGrid[N+(j-1)*gridSize[1]+i] = (staggeredGrid[M+(j-1)*(gridSize[1]+1)+i] + staggeredGrid[M+(j-1)*(gridSize[1]+1)+i+1]) / 2
#         end
#     end
#     return centeredGrid
# end
#
# # Cellcentered grid to staggeredGrid
# function cen2stg(centeredGrid::Array{Float64,1},gridSize::Array{Int64,1})
#     NX = (gridSize[2]+1)*gridSize[1]
#     NY = (gridSize[1]+1)*gridSize[2]
#     M = prod(gridSize)
#     staggeredGrid = zeros(NX+NY)
#
#     # staggered grid in horizontal direction
#     # inner part
#     for i=gridSize[1]+1:M
#         staggeredGrid[i] = ( centeredGrid[i-gridSize[1]] + centeredGrid[i] ) / 2
#     end
#     # outer part (edge area)
#     for i=1:gridSize[1]
#         staggeredGrid[i] = centeredGrid[i] / 2
#     end
#     for i=M+1:NX
#         staggeredGrid[i] = centeredGrid[i-gridSize[1]] / 2
#     end
#
#     # staggered grid in vertical direction
#     # inner part
#     for j=1:gridSize[2]
#         staggeredGrid[NX+(j-1)*(gridSize[1]+1)+1] = centeredGrid[M+1+(j-1)*gridSize[1]] / 2
#         for i=2:gridSize[1]
#             staggeredGrid[NX+(j-1)*(gridSize[1]+1)+i] = (centeredGrid[M+i+(j-1)*gridSize[1]-1] + centeredGrid[M+i+(j-1)*gridSize[1]]) / 2
#         end
#         staggeredGrid[NX+j*(gridSize[1]+1)] = centeredGrid[M+j*gridSize[1]] / 2
#     end
#
#     return staggeredGrid
# end
#
# function cen2stg(distanceOutput::Tuple{Number,Array{Float64,1},Function,Tuple{Array{Float64,1},Array{Float64,1}}}, refImg::ImageMeta)
#   D = distanceOutput[1]
#   dD = cen2stg(distanceOutput[2],getSize(refImg))
#   d2D(grid) = cen2stg(distanceOutput[3](stg2cen(grid,getSize(refImg))),getSize(refImg))
#   dTransformedImage = distanceOutput[4]
#   return(D, dD, d2D, dTransformedImage)
# end
#
# function cen2stg(distanceOutput::Tuple{Number,Array{Float64,1},Function}, refImg::ImageMeta)
#   D = distanceOutput[1]
#   dD = cen2stg(distanceOutput[2],getSize(refImg))
#   d2D(grid) = cen2stg(distanceOutput[3](stg2cen(grid,getSize(refImg))),getSize(refImg))
# 	return(D, dD, d2D)
# end

# function cen2stg(distanceOutput::Tuple{Number,Array{Float64,1},SparseMatrixCSC{Float64,Int64},Tuple{Array{Float64,1},Array{Float64,1}}}, refImg::ImageMeta)
#   D = distanceOutput[1]
#   dD = cen2stg(distanceOutput[2],getSize(refImg))
#   d2D(grid) = cen2stg(distanceOutput[3]*stg2cen(grid,getSize(refImg)),getSize(refImg))
#   dTransformedImage = distanceOutput[4]
#   return(D, dD, d2D, dTransformedImage)
# end
#
# function cen2stg(distanceOutput::Tuple{Number,Array{Float64,1},SparseMatrixCSC{Float64,Int64}}, refImg::ImageMeta)
#   D = distanceOutput[1]
#   dD = cen2stg(distanceOutput[2],getSize(refImg))
#   d2D(grid) = cen2stg(distanceOutput[3]*stg2cen(grid,getSize(refImg)),getSize(refImg))
#   return(D, dD, d2D)
# end
#
# function cen2stg(distanceOutput::Tuple{Number,Any,Any,Tuple{Array{Float64,1},Array{Float64,1}}}, refImg::ImageMeta)
#   D = distanceOutput[1]
#   dD = distanceOutput[2]
#   d2D = distanceOutput[3]
#   dTransformedImage = distanceOutput[4]
#   return(D, dD, d2D, dTransformedImage)
# end
#
# function cen2stg(distanceOutput::Tuple{Number,Number,Number}, refImg::ImageMeta)
#   D = distanceOutput[1]
#   dD = distanceOutput[2]
#   d2D = distanceOutput[3]
#   return(D, dD, d2D)
# end
#
function checkForOddNumberOfGridPoints(I::ImageMeta,grid::Array{Float64,1})
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
