module Interpolation
using Images
using Interpolations
import Interpolations.Throw # see https://github.com/JuliaMath/Interpolations.jl/issues/127


using DeformableRegistration.Transformation
using DeformableRegistration.ImageProcessing

export interpolateImage,InterpLinearFast
export interpolateDeformationField

function interpolateImage(
    I::regImage,transformedGrid::scaledArray;
    doDerivative=false,interpolationScheme=InterpLinearFast)

    return interpolateImage(I, transformedGrid.data, doDerivative=doDerivative,
                            interpolationScheme=interpolationScheme)
end

function interpolateImage(
    I::regImage,transformedGrid::Array{Float64,1};
    doDerivative=false,interpolationScheme=InterpLinearFast)

    return interpolateArray(Array(I.data), I.voxelsize, I.shift, transformedGrid, size(I.data);
                            doDerivative=doDerivative,interpolationScheme=interpolationScheme)
end

function interpolateArray(
    data::Array{Float64,2}, voxelsize::Array{Float64,1},
    shift::Array{Float64,1}, transformedGrid::Array{Float64,1}, newDimensions::Tuple{Vararg{Int64}};
    doDerivative=false, interpolationScheme=InterpLinearFast)

  # use faster linear interpolation as default
  if(interpolationScheme == InterpLinearFast)
    return InterpLinearFast(data, voxelsize, shift, transformedGrid, newDimensions, doDerivative=doDerivative)
  end

  # determine number of new points, pixel spacing and spatial domain of the image
  numberOfPoints::Int = ceil(Int64,size(transformedGrid,1)/2)

  imageIntNoScaling = interpolate(data, interpolationScheme, OnCell())

  xRange = shift[1]+voxelsize[1]/2:voxelsize[1]:shift[1]+voxelsize[1]*(size(data,1))
  yRange = shift[2]+voxelsize[1]/2:voxelsize[2]:shift[2]+voxelsize[2]*(size(data,2))

  imageIntTmp = extrapolate(imageIntNoScaling, 0.0)
  imageInt = scale(imageIntTmp, xRange, yRange)


  transformedImage = zeros(numberOfPoints);
  for i=1:numberOfPoints
      transformedImage[i] = imageInt[transformedGrid[i],transformedGrid[numberOfPoints+i]]
  end

  if(doDerivative)
    imageIntTmp = extrapolate(imageIntNoScaling, Throw())
    imageInt = scale(imageIntTmp, xRange, yRange)

    dY_transformedImage = zeros(numberOfPoints)
    dX_transformedImage = zeros(numberOfPoints)
    for i=1:numberOfPoints
        try # zero if outside boundary (extrapolate(..., 0.0) does not work with gradient)
            dX_transformedImage[i],dY_transformedImage[i] = gradient(imageInt,transformedGrid[i],transformedGrid[i+numberOfPoints])
        catch x
            if isa(x, BoundsError)
                #print(".")
            else
                throw(x)
            end
        end
    end
    return reshape(transformedImage, newDimensions), dX_transformedImage, dY_transformedImage
  end
  return reshape(transformedImage, newDimensions)
end

function InterpLinearFast(data::Array{Float64,2}, voxelsize::Array{Float64,1}, shift::Array{Float64,1},
                            transformedGrid::Array{Float64,1}, newDimensions::Tuple{Vararg{Int64}};doDerivative=false)

  # determine number of new points, pixel spacing and spatial domain of the image
  numberOfPoints::Int = convert(Int64,(size(transformedGrid,1)/2))

  # setup variables
  transformedImage = zeros(numberOfPoints)
  dY_transformedImage = zeros(numberOfPoints)
  dX_transformedImage = zeros(numberOfPoints)
  p = zeros(4)
  sizeY= size(data)[2]
  sizeX = size(data)[1]
  x::Float64 = 0.0; y::Float64 = 0.0;
  xf::Int = 0; yf::Int = 0;

  # interpolate image at new points
  for i=1:numberOfPoints

    x = (transformedGrid[i]   - shift[1]) / voxelsize[1] + .5
    y = (transformedGrid[numberOfPoints+i] - shift[2]) / voxelsize[2] + .5

    if ((x <= 0) | (y <= 0) | (x >= sizeX+1) | (y >= sizeY+1))
      continue
    end

    xf = floor(Int64,x)
    yf = floor(Int64,y)

    p[1] = ((xf<1)          | (yf<1))        ?   0.0: data[xf,yf]
    p[2] = ((xf+1>sizeX)    | (yf<1))        ?   0.0: data[xf+1,yf]
    p[3] = ((xf<1)          | (yf+1>sizeY))  ?   0.0: data[xf,yf+1]
    p[4] = ((xf+1>sizeX)    | (yf+1>sizeY))  ?   0.0: data[xf+1,yf+1]

    x = x - xf
    y = y - yf

    transformedImage[i] = (p[1] * (1-x) + p[2] * x) * (1-y) + (p[3] * (1-x) + p[4] * x) * y
    dX_transformedImage[i] = ((p[3] - p[1]) * (1-x) + (p[4] - p[2]) * x) / voxelsize[1]
    dY_transformedImage[i] = ((p[2] - p[1]) * (1-y) + (p[4] - p[3]) * y) / voxelsize[2]

  end

  if(doDerivative)
    return reshape(transformedImage,newDimensions),dY_transformedImage,dX_transformedImage
  else
    return reshape(transformedImage,newDimensions)
  end
end


function interpolateDeformationField(deformationField::scaledArray, newPoints::scaledArray; interpolationScheme = InterpLinearFast)# BSpline(Linear())
    dims = deformationField.dimensions
    deformationFieldX = reshape(deformationField.data[1:prod(dims)], dims[1], dims[2])
    deformationFieldY = reshape(deformationField.data[prod(dims)+1:end], dims[1], dims[2])

    return scaledArray([
        interpolateArray(deformationFieldX, deformationField.voxelsize, deformationField.shift,
            newPoints.data, (prod(newPoints.dimensions),), interpolationScheme = interpolationScheme, doDerivative=false);
        interpolateArray(deformationFieldY, deformationField.voxelsize, deformationField.shift,
            newPoints.data, (prod(newPoints.dimensions),), interpolationScheme = interpolationScheme, doDerivative=false)
        ],newPoints.dimensions, newPoints.voxelsize, newPoints.shift)


end

end
