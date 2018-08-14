module Interpolation
using Interpolations
import Interpolations.Throw # see https://github.com/JuliaMath/Interpolations.jl/issues/127

# using Logging.debug
using DeformableRegistration.Transformation
using DeformableRegistration.Types


export interpolateImage,InterpLinearFast, interpolateArray
export interpolateDeformationField

function interpolateImage(
    I::regImage,transformedGrid::scaledArray;
    doDerivative=false,interpolationScheme=InterpLinearFast)

    return interpolateArray(Array(I.data), I.voxelsize, I.shift, transformedGrid.data, transformedGrid.dimensions;
                            doDerivative=doDerivative,interpolationScheme=interpolationScheme)
end

function interpolateImage(
    I::regImage,transformedGrid::scaledArray, targetGrid::scaledArray;
    doDerivative=false,interpolationScheme=BSpline(Cubic(Line())))

    return interpolateArray(Array(I.data), I.voxelsize, I.shift, transformedGrid, targetGrid;
                            doDerivative=doDerivative,interpolationScheme=interpolationScheme)
end

function interpolateImage(
    I::regImage,transformedGrid::Array{Float64,1};
    doDerivative=false,interpolationScheme=InterpLinearFast)
    @assert prod(size(I.data.data))==length(transformedGrid)/2
    return interpolateArray(Array(I.data), I.voxelsize, I.shift, transformedGrid, size(I.data);
                            doDerivative=doDerivative,interpolationScheme=interpolationScheme)
end
##
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
  yRange = shift[2]+voxelsize[2]/2:voxelsize[2]:shift[2]+voxelsize[2]*(size(data,2))

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
        if (xRange[1] <= transformedGrid[i] < xRange[end] + voxelsize[1]) &&
            (yRange[1] <= transformedGrid[i+numberOfPoints] < yRange[end] + voxelsize[2])
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
    end
    return reshape(transformedImage, newDimensions), dX_transformedImage, dY_transformedImage
  end
  return reshape(transformedImage, newDimensions),0,0
end
##
function interpolateArray(
    data::Array{Float64,2}, voxelsize::Array{Float64,1},
    shift::Array{Float64,1}, transformedGrid::scaledArray,
    targetGrid::scaledArray;
    doDerivative=false, interpolationScheme=BSpline(Cubic(Line())))

  # determine number of new points, pixel spacing and spatial domain of the image
  imageIntNoScaling = interpolate(data, interpolationScheme, OnCell())

  xRange = shift[1]+voxelsize[1]/2:voxelsize[1]:shift[1]+voxelsize[1]*(size(data,1))
  yRange = shift[2]+voxelsize[2]/2:voxelsize[2]:shift[2]+voxelsize[2]*(size(data,2))

  imageIntTmp = extrapolate(imageIntNoScaling, 0.0)
  imageInt = scale(imageIntTmp, xRange, yRange)

  deformationAtTargetGrid = interpolateDeformationField(transformedGrid-getCellCenteredGrid(transformedGrid), targetGrid, interpolationScheme = interpolationScheme) + getCellCenteredGrid(targetGrid)

  numberOfPoints::Int = prod(deformationAtTargetGrid.dimensions)

  transformedImage = zeros(numberOfPoints);
  for i=1:numberOfPoints
      transformedImage[i] = imageInt[deformationAtTargetGrid.data[i],deformationAtTargetGrid.data[numberOfPoints+i]]
  end
  if(doDerivative)
    imageIntTmp = extrapolate(imageIntNoScaling, Throw())
    imageInt = scale(imageIntTmp, xRange, yRange)

    dY_transformedImage = zeros(numberOfPoints)
    dX_transformedImage = zeros(numberOfPoints)
    for i=1:numberOfPoints
        # first checking range ourselves (much faster then try-catch but result is different (check for . between comma in output))
        if (xRange[1] <= deformationAtTargetGrid.data[i] < xRange[end] + voxelsize[1]) &&
            (yRange[1] <= deformationAtTargetGrid.data[i+numberOfPoints] < yRange[end] + voxelsize[2])
            try # zero if outside boundary (extrapolate(..., 0.0) does not work with gradient)
                dX_transformedImage[i],dY_transformedImage[i] = gradient(imageInt,deformationAtTargetGrid.data[i],deformationAtTargetGrid.data[i+numberOfPoints])
            catch x
                if isa(x, BoundsError)
                    #print(".")
                else
                    throw(x)
                end
            end
        else
            #print(",")
        end

    end

    return reshape(transformedImage, targetGrid.dimensions), dX_transformedImage, dY_transformedImage
  end
  return reshape(transformedImage, targetGrid.dimensions),0,0
end
##

function InterpLinearFast(data::Array{Float64,2}, voxelsize::Array{Float64,1}, shift::Array{Float64,1},
                            transformedGrid::Array{Float64,1}, newDimensions::Tuple{Vararg{Int64}};doDerivative=false)

  # determine number of new points, pixel spacing and spatial domain of the image
  numberOfPoints::Int = convert(Int64,(size(transformedGrid,1)/2))

  # setup variables
  transformedImage = zeros(numberOfPoints)
  dX_transformedImage = zeros(numberOfPoints)
  dY_transformedImage = zeros(numberOfPoints)
  p = zeros(4)
  sizeX = size(data)[1]
  sizeY= size(data)[2]
  x::Float64 = 0.0; y::Float64 = 0.0;
  xf::Int = 0; yf::Int = 0;

  # interpolate image at new points
  for i=1:numberOfPoints

    x = (transformedGrid[i]   - shift[1]) / voxelsize[1] + .5
    y = (transformedGrid[numberOfPoints+i] - shift[2]) / voxelsize[2] + .5

    if ((x <= 0) | (y <= 0) | (x >= sizeX+1) | (y >= sizeY+1))
      continue
    end

    xf = floor(Int64,x);
    yf = floor(Int64,y);

    # p[1] = ((xf<1)          | (yf<1))        ?   0.0: data[xf,yf]
    # p[2] = ((xf+1>sizeX)    | (yf<1))        ?   0.0: data[xf+1,yf]
    # p[3] = ((xf<1)          | (yf+1>sizeY))  ?   0.0: data[xf,yf+1]
    # p[4] = ((xf+1>sizeX)    | (yf+1>sizeY))  ?   0.0: data[xf+1,yf+1]
    boundsX(v) = max(1, min(sizeX, v))
    boundsY(v) = max(1, min(sizeY, v))

    p[1] = data[boundsX(xf), boundsY(yf)]
    p[2] = data[boundsX(xf+1), boundsY(yf)]
    p[3] = data[boundsX(xf), boundsY(yf+1)]
    p[4] = data[boundsX(xf+1), boundsY(yf+1)]

    x = x - xf
    y = y - yf

    transformedImage[i] = (p[1] * (1-x) + p[2] * x) * (1-y) + (p[3] * (1-x) + p[4] * x) * y
    dX_transformedImage[i] = ((p[3] - p[1]) * (1-x) + (p[4] - p[2]) * x) / voxelsize[1]
    dY_transformedImage[i] = ((p[2] - p[1]) * (1-y) + (p[4] - p[3]) * y) / voxelsize[2]

  end

  if(doDerivative)
    return reshape(transformedImage,newDimensions),dY_transformedImage,dX_transformedImage
  else
    return reshape(transformedImage,newDimensions),0,0
  end
end


function interpolateDeformationField(deformationField::scaledArray,
        newPoints::scaledArray;
        interpolationScheme = InterpLinearFast)# BSpline(Linear())

    dims = deformationField.dimensions
    deformationFieldX = reshape(deformationField.data[1:prod(dims)], dims[1], dims[2])
    deformationFieldY = reshape(deformationField.data[prod(dims)+1:end], dims[1], dims[2])

    return scaledArray([
        interpolateArray(deformationFieldX, deformationField.voxelsize, deformationField.shift,
            newPoints.data, (prod(newPoints.dimensions),), interpolationScheme = interpolationScheme, doDerivative=false)[1];
        interpolateArray(deformationFieldY, deformationField.voxelsize, deformationField.shift,
            newPoints.data, (prod(newPoints.dimensions),), interpolationScheme = interpolationScheme, doDerivative=false)[1]
        ],newPoints.dimensions, newPoints.voxelsize, newPoints.shift)


end

end
