module Interpolation

using Images
using Grid

using ImageRegistration.ImageProcessing

export interpolateDeformationFieldAtGrid
export interpolateImageAtGrid, interpolateImageAtGridWithDerivative
export linearImageInterpolationAtGrid, linearImageInterpolationAtGridWithDerivative

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
    pixelSpacing::Array{Float64,1} = getPixelSpacing(image)
    spatialDomain::Array{Float64,1} = getSpatialDomain(image)
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
    pixelSpacing::Array{Float64,1} = getPixelSpacing(image)
    spatialDomain::Array{Float64,1} = getSpatialDomain(image)
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
    pixelSpacing::Array{Float64,1} = getPixelSpacing(image)
    spatialDomain::Array{Float64,1} = getSpatialDomain(image)

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
    pixelSpacing::Array{Float64,1} = getPixelSpacing(image)
    spatialDomain::Array{Float64,1} = getSpatialDomain(image)

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

    return transformedImage,dY_transformedImage,dX_transformedImage

end

end
