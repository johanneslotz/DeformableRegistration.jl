module Transformation

using Images
using Grid

using ImageRegistration.ImageProcessing

export getCellCenteredGrid, getStaggeredGrid, stg2cen, cen2stg, checkStaggered, getCellCenteredGridRanges
export transformGridAffine, transformWorldCoordinate

function getCellCenteredGrid(I::Image)
  return getCellCenteredGrid(getSpatialDomain(I),getSize(I))
end

function getCellCenteredGrid(spatialDomain::Array{Float64,1},gridSize::Array{Int,1})
    pixelSpacing = [(spatialDomain[2]-spatialDomain[1])/gridSize[1],
                    (spatialDomain[4]-spatialDomain[3])/gridSize[2]]
    y = [spatialDomain[1]+pixelSpacing[1]/2:pixelSpacing[1]:spatialDomain[2]]
    x = [spatialDomain[3]+pixelSpacing[2]/2:pixelSpacing[2]:spatialDomain[4]]
    return [repmat(x,1,gridSize[1])'[:],repmat(y,1,gridSize[2])[:]]
end

function getStaggeredGrid(image::Image)
  return getStaggeredGrid(getSpatialDomain(I),getSize(I))
end

function getStaggeredGrid(spatialDomain::Array{Float64,1},gridSize::Array{Int,1})
    pixelSpacing = [(spatialDomain[2]-spatialDomain[1])/gridSize[1],
                    (spatialDomain[4]-spatialDomain[3])/gridSize[2]]
    y = [spatialDomain[1]:pixelSpacing[1]:spatialDomain[2]]
    x = [spatialDomain[3]:pixelSpacing[2]:spatialDomain[4]]
    return [repmat(x,1,gridSize[1])'[:],repmat(y,1,gridSize[2])[:]]
end

function stg2cen(staggeredGrid::Array{Float64,1},gridSize::Array{Int,1})
    N = prod(gridSize)
    M = gridSize[1]*(gridSize[2]+1)
    centeredGrid = zeros(2*N)
    for j=1:gridSize[2]
        for i=1:gridSize[1]
            centeredGrid[  (j-1)*gridSize[1]+i] = (staggeredGrid[(j-1)* gridSize[1]+i]      + staggeredGrid[j*gridSize[1]+i]            ) / 2
            centeredGrid[N+(j-1)*gridSize[1]+i] = (staggeredGrid[M+(j-1)*(gridSize[1]+1)+i] + staggeredGrid[M+(j-1)*(gridSize[1]+1)+i+1]) / 2
        end
    end
    return centeredGrid
end

function cen2stg(centeredGrid::Array{Float64,1},gridSize::Array{Int64,1})
    NX = (gridSize[2]+1)*gridSize[1]
    NY = (gridSize[1]+1)*gridSize[2]
    M = prod(gridSize)
    staggeredGrid = zeros(NX+NY)

    # staggered grid in horizontal direction
    # inner part
    for i=gridSize[1]+1:M
        staggeredGrid[i] = ( centeredGrid[i-gridSize[1]] + centeredGrid[i] ) / 2
    end
    # outer part (edge area)
    for i=1:gridSize[1]
        staggeredGrid[i] = centeredGrid[i] / 2
    end
    for i=M+1:NX
        staggeredGrid[i] = centeredGrid[i-gridSize[1]] / 2
    end

    # staggered grid in vertical direction
    # inner part
    for j=1:gridSize[2]
        staggeredGrid[NX+(j-1)*(gridSize[1]+1)+1] = centeredGrid[M+1+(j-1)*gridSize[1]] / 2
        for i=2:gridSize[1]
            staggeredGrid[NX+(j-1)*(gridSize[1]+1)+i] = (centeredGrid[M+i+(j-1)*gridSize[1]-1] + centeredGrid[M+i+(j-1)*gridSize[1]]) / 2
        end
        staggeredGrid[NX+j*(gridSize[1]+1)] = centeredGrid[M+j*gridSize[1]] / 2
    end

    return staggeredGrid
end

function cen2stg(distanceOutput::(Number,Array{Float64,1},Function,(Array{Float64,1},Array{Float64,1})), refImg::Image)
  D = distanceOutput[1]
  dD = cen2stg(distanceOutput[2],getSize(refImg))
  d2D(grid) = cen2stg(distanceOutput[3](stg2cen(grid,getSize(refImg))),getSize(refImg))
  dTransformedImage = distanceOutput[4]
  return(D, dD, d2D, dTransformedImage)
end

function cen2stg(distanceOutput::(Number,Array{Float64,1},SparseMatrixCSC{Float64,Int64},(Array{Float64,1},Array{Float64,1})), refImg::Image)
  D = distanceOutput[1]
  dD = cen2stg(distanceOutput[2],getSize(refImg))
  d2D(grid) = cen2stg(distanceOutput[3]*stg2cen(grid,getSize(refImg)),getSize(refImg))
  dTransformedImage = distanceOutput[4]
  return(D, dD, d2D, dTransformedImage)
end

function cen2stg(distanceOutput::(Number,Array{Float64,1},Matrix,(Array{Float64,1},Array{Float64,1})), refImg::Image)
  D = distanceOutput[1]
  dD = distanceOutput[2]
  d2D = distanceOutput[3]
  dTransformedImage = distanceOutput[4]
  return(D, dD, d2D, dTransformedImage)
end

function cen2stg(distanceOutput::(Number,Number,Number,(Array{Float64,1},Array{Float64,1})), refImg::Image)
  D = distanceOutput[1]
  dD = distanceOutput[2]
  d2D = distanceOutput[3]
  dTransformedImage = distanceOutput[4]
  return(D, dD, d2D, dTransformedImage)
end

function checkStaggered(I::Image,grid::Array{Float64,1})
  if( 2*prod(size(I)) == size(grid,1) )
      return false
  else
      return true
  end
end

function getCellCenteredGridRanges(I::Image)
    spatialDomain = getSpatialDomain(I)
    pixelSpacing = getPixelSpacing(I)
    y = spatialDomain[1]+pixelSpacing[1]/2:pixelSpacing[1]:spatialDomain[2]
    x = spatialDomain[3]+pixelSpacing[2]/2:pixelSpacing[2]:spatialDomain[4]
    return y,x
end

function transformWorldCoordinate(I::Image,y::Float64,x::Float64)
    spatialDomain = getSpatialDomain(I)
    pixelSpacing = getPixelSpacing(I)
    yNew = (y - spatialDomain[1]) / pixelSpacing[1] + .5
    xNew = (x - spatialDomain[3]) / pixelSpacing[2] + .5
    return yNew,xNew
end

function transformWorldCoordinate(I::Image,grid::Array{Float64,1})
    spatialDomain = getSpatialDomain(I)
    pixelSpacing = getPixelSpacing(I)
    N = size(grid,1)
    newGrid = zeros(N)
    for i=1:N
      newGrid[i] = (grid[i] - spatialDomain[3]) / pixelSpacing[2] + .5
      newGrid[N+i] = (grid[N+i] - spatialDomain[1]) / pixelSpacing[1] + .5
    end
    return newGrid
end

function transformGridAffine(grid,affineParameters)
    transformationMatrix = [affineParameters[1] affineParameters[2]; affineParameters[4] affineParameters[5]]
    translationParameters = [affineParameters[3]; affineParameters[6]]
    return transformGridAffine(transformationMatrix,grid,translationParameters)
end

function transformGridAffine(transformationMatrix::Matrix,gridOrg::Array{Float64,1},translationParameters::Array{Float64,1})
    N = int(size(gridOrg,1)/2)
    grid = zeros(2*N)
    @simd for i=1:N
        @inbounds grid[i]   = transformationMatrix[1,1]*gridOrg[i]+transformationMatrix[1,2]*gridOrg[N+i]+translationParameters[1]
        @inbounds grid[N+i] = transformationMatrix[2,1]*gridOrg[i]+transformationMatrix[2,2]*gridOrg[N+i]+translationParameters[2]
    end
    return grid
end

end
