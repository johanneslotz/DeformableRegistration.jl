module Transformation

using Images
using Grid

export getCellCenteredGrid, getStaggeredGrid, stg2cen, cen2stg, checkStaggered
export transformGridAffine
export interpolateDeformationFieldAtGrid
export interpolateImageAtGrid, interpolateImageAtGridWithDerivative
export linearImageInterpolationAtGrid, linearImageInterpolationAtGridWithDerivative

function getCellCenteredGrid(image::Image)
  return getCellCenteredGrid(image.properties["spatialdomain"],[size(image)[1],size(image)[2]])
end

function getCellCenteredGrid(spatialDomain::Array{Float64,1},gridSize::Array{Int,1})
    pixelSpacing = [(spatialDomain[2]-spatialDomain[1])/gridSize[1],
                    (spatialDomain[4]-spatialDomain[3])/gridSize[2]]
    y = [spatialDomain[1]+pixelSpacing[1]/2:pixelSpacing[1]:spatialDomain[2]]
    x = [spatialDomain[3]+pixelSpacing[2]/2:pixelSpacing[2]:spatialDomain[4]]
    return [repmat(x,1,gridSize[1])'[:],repmat(y,1,gridSize[2])[:]]
end

function getStaggeredGrid(image::Image)
  return getStaggeredGrid(image.properties["spatialdomain"],[size(image)[1],size(image)[2]])
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

function cen2stg(distanceOutput::(Number,Array{Float64,1},Any,(Array{Float64,1},Array{Float64,1})), refImg::Image)
  D = distanceOutput[1]
  dD = distanceOutput[2]
  d2D = distanceOutput[3]
  dTransformedImage = distanceOutput[4]
  dD = cen2stg(dD,[size(refImg)[1], size(refImg)[2]])
  pixelSpacing = prod(refImg.properties["pixelspacing"])
  d2D = calculateHessianStaggered((dTransformedImage), [size(refImg)[1], size(refImg)[2]], pixelSpacing)
  return(D, dD, d2D, dTransformedImage)
end

function cen2stg(distanceOutput::(Number,Any,Any,(Array{Float64,1},Array{Float64,1})), refImg::Image)
  D = distanceOutput[1]
  dD = distanceOutput[2]
  d2D = distanceOutput[3]
  dTransformedImage = distanceOutput[4]
  return(D, dD, d2D, dTransformedImage)
end

function calculateHessianStaggered(dTransformedImage, sizeOfReferenceImage, pixelSpacing)
  drstg = cen2stg( [dTransformedImage[1]; dTransformedImage[2]] , sizeOfReferenceImage)
  mYstg = (sizeOfReferenceImage[1]+1)*sizeOfReferenceImage[2]; mXstg = sizeOfReferenceImage[1]*(sizeOfReferenceImage[2]+1)
  dTransformedImage = spdiagm((drstg[1:mXstg],drstg[mXstg+1:end]),[0,mXstg])
  d2FunctionValue   = pixelSpacing .* dTransformedImage' * dTransformedImage
  return d2FunctionValue
end

function checkStaggered(image,grid)
  if( (2*prod(size(image))) == size(grid,1) )
      return false
  else
      return true
  end
end

function transformGridAffine(grid,affineParameters)
    transformationMatrix = [affineParameters[1] affineParameters[2]; affineParameters[4] affineParameters[5]]
    translationParameters = [affineParameters[3]; affineParameters[6]]
    return transformGridAffine(transformationMatrix,grid,translationParameters)
end

function transformGridAffine(transformationMatrix,gridOrg,translationParameters)
    grid = copy(gridOrg)
    grid=reshape(grid,int(size(grid,1)/2),2)
    for i=1:int(size(grid,1))
        grid[i,:] = (transformationMatrix*grid[i,:]'+translationParameters)'
    end
    return grid[:]
end

function interpolateDeformationFieldAtGrid(deformationField::Array{Float64,2}, spatialDomain::Array{Float64,1}, newPoints::Array{Float64,1})

    # determine number of new points and pixel spacing
    numberOfPoints::Int = int(size(newPoints,1)/2)
    pixelSpacing::Array{Float64,1} = [(spatialDomain[2]-spatialDomain[1])/size(deformationField)[1],
                                      (spatialDomain[4]-spatialDomain[3])/size(deformationField)[2]]

    # define interpolation grid
    deformationFieldInt = InterpGrid(deformationField, BCreflect, InterpLinear)

    # interpolate deformation at new points
    deformationFieldAtNewPoints=Array(Float64, numberOfPoints); y=0.0; x=0.0;
    for i=1:numberOfPoints
        x = (newPoints[i] - spatialDomain[3]) / pixelSpacing[2] + .5
        y = (newPoints[numberOfPoints+i] - spatialDomain[1]) / pixelSpacing[1] + .5
        deformationFieldAtNewPoints[i] = deformationFieldInt[y,x]
    end

    return deformationFieldAtNewPoints

end

function interpolateDeformationFieldAtGrid(deformationField::Array{Float64,1}, gridSize, spatialDomain::Array{Float64,1}, newPoints::Array{Float64,1})

    deformationFieldX = reshape(deformationField[1:prod(gridSize)],gridSize[1],gridSize[2])
    deformationFieldY = reshape(deformationField[prod(gridSize)+1:end],gridSize[1],gridSize[2])

    return [interpolateDeformationFieldAtGrid(deformationFieldX,spatialDomain,newPoints),
            interpolateDeformationFieldAtGrid(deformationFieldY,spatialDomain,newPoints)]

end

function interpolateImageAtGrid(image::Image,transformedGrid::Array{Float64,1};
                                         InterpFunction=InterpLinear)

    # determine number of new points, pixel spacing and spatial domain of the image
    numberOfPoints::Int = int(size(transformedGrid,1)/2)
    pixelSpacing::Array{Float64,1} = pixelspacing(image)
    spatialDomain::Array{Float64,1} = image.properties["spatialdomain"]
    sizeY = size(image)[1]
    sizeX = size(image)[2]

    # define interpolation grid
    imageInt = InterpGrid(image.data, 0.0, InterpFunction)

    # interpolate image at new points
    transformedImage = zeros(numberOfPoints); y=0.0; x=0.0;
    for i=1:numberOfPoints
        x = (transformedGrid[i] - spatialDomain[3]) / pixelSpacing[2] + .5
        y = (transformedGrid[numberOfPoints+i] - spatialDomain[1]) / pixelSpacing[1] + .5
        if ((x <= 0) | (y <= 0) | (x >= sizeX+1) | (y >= sizeY+1))
            continue
        end
        transformedImage[i] = imageInt[y,x]
    end

    return reshape(transformedImage,height(image),width(image))

end

function interpolateImageAtGridWithDerivative(image::Image,transformedGrid::Array{Float64,1};
                                                       InterpFunction=InterpLinear)

    # determine number of new points, pixel spacing and spatial domain of the image
    numberOfPoints::Int = int(size(transformedGrid,1)/2)
    pixelSpacing::Array{Float64,1} = pixelspacing(image)
    spatialDomain::Array{Float64,1} = image.properties["spatialdomain"]
    sizeY = size(image)[1]
    sizeX = size(image)[2]

    # define interpolation grid
    imageInt = InterpGrid(image.data, 0.0, InterpFunction)

    # interpolate image at new points with derivative
    transformedImage = zeros(numberOfPoints); y=0.0; x=0.0;
    dY_transformedImage = zeros(numberOfPoints)
    dX_transformedImage = zeros(numberOfPoints)
    for i=1:numberOfPoints
        x = (transformedGrid[i] - spatialDomain[3]) / pixelSpacing[2] + .5
        y = (transformedGrid[numberOfPoints+i] - spatialDomain[1]) / pixelSpacing[1] + .5
        if ((x <= 0) | (y <= 0) | (x >= sizeX+1) | (y >= sizeY+1))
            continue
        end
        transformedImage[i],(dY_transformedImage[i],dX_transformedImage[i]) = valgrad(imageInt,y,x)
    end

    return reshape(transformedImage,height(image),width(image)),
           reshape(dY_transformedImage,height(image),width(image)),
           reshape(dX_transformedImage,height(image),width(image))
end

function linearImageInterpolationAtGrid(image::Image,transformedGrid::Array{Float64,1})

    # determine number of new points, pixel spacing and spatial domain of the image
    numberOfPoints::Int = int(size(transformedGrid,1)/2)
    pixelSpacing::Array{Float64,1} = pixelspacing(image)
    spatialDomain::Array{Float64,1} = image.properties["spatialdomain"]

    # interpolate image at new points
    transformedImage = zeros(numberOfPoints)
    p::Array{Float64,1} = zeros(4)
    sizeY = size(image)[1]
    sizeX = size(image)[2]
    x::Float64 = 0.0; y::Float64 = 0.0;
    xf::Int = 0; yf::Int = 0;

    for i=1:numberOfPoints

        x = (transformedGrid[i]   - spatialDomain[3]) / pixelSpacing[2] + .5
        y = (transformedGrid[numberOfPoints+i] - spatialDomain[1]) / pixelSpacing[1] + .5

        if ((x <= 0) | (y <= 0) | (x >= sizeX+1) | (y >= sizeY+1))
            continue
        end

        xf = ifloor(x)
        yf = ifloor(y)


        p[1] = ((xf<1)          | (yf<1))        ?   0.0: image.data[yf,xf]
        p[2] = ((xf+1>sizeX)    | (yf<1))        ?   0.0: image.data[yf,xf+1]
        p[3] = ((xf<1)          | (yf+1>sizeY))  ?   0.0: image.data[yf+1,xf]
        p[4] = ((xf+1>sizeX)    | (yf+1>sizeY))  ?   0.0: image.data[yf+1,xf+1]

        x = x - xf
        y = y - yf

        transformedImage[i] = (p[1] * (1-x) + p[2] * x) * (1-y) + (p[3] * (1-x) + p[4] * x) * y

    end

    return reshape(transformedImage,height(image),width(image))

end

function linearImageInterpolationAtGridWithDerivative(image::Image,transformedGrid::Array{Float64,1})

    # determine number of new points, pixel spacing and spatial domain of the image
    numberOfPoints::Int = int(size(transformedGrid,1)/2)
    pixelSpacing::Array{Float64,1} = pixelspacing(image)
    spatialDomain::Array{Float64,1} = image.properties["spatialdomain"]

    # interpolate image at new points
    transformedImage = zeros(numberOfPoints)
    dY_transformedImage = zeros(numberOfPoints)
    dX_transformedImage = zeros(numberOfPoints)
    p = zeros(4)
    sizeY= size(image)[1]
    sizeX = size(image)[2]
    x::Float64 = 0.0; y::Float64 = 0.0;
    xf::Int = 0; yf::Int = 0;

    for i=1:numberOfPoints

        x = (transformedGrid[i]   - spatialDomain[3]) / pixelSpacing[2] + .5
        y = (transformedGrid[numberOfPoints+i] - spatialDomain[1]) / pixelSpacing[1] + .5

        if ((x <= 0) | (y <= 0) | (x >= sizeX+1) | (y >= sizeY+1))
            continue
        end

        xf = ifloor(x)
        yf = ifloor(y)

        p[1] = ((xf<1)          | (yf<1))        ?   0.0: image.data[yf,xf]
        p[2] = ((xf+1>sizeX)    | (yf<1))        ?   0.0: image.data[yf,xf+1]
        p[3] = ((xf<1)          | (yf+1>sizeY))  ?   0.0: image.data[yf+1,xf]
        p[4] = ((xf+1>sizeX)    | (yf+1>sizeY))  ?   0.0: image.data[yf+1,xf+1]

        x = x - xf
        y = y - yf

        transformedImage[i] = (p[1] * (1-x) + p[2] * x) * (1-y) + (p[3] * (1-x) + p[4] * x) * y
        dX_transformedImage[i] = ((p[2] - p[1]) * (1-y) + (p[4] - p[3]) * y) / pixelSpacing[2]
        dY_transformedImage[i] = ((p[3] - p[1]) * (1-x) + (p[4] - p[2]) * x) / pixelSpacing[1]

    end

    return reshape(transformedImage,height(image),width(image)),
           reshape(dY_transformedImage,height(image),width(image)),
           reshape(dX_transformedImage,height(image),width(image))

end

end
