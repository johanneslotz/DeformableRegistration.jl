module ImageProcessing

using FileIO
using ColorTypes
using DeformableRegistration.Types
using DeformableRegistration.Transformation
using DeformableRegistration.Interpolation
using MicroLogging
export createImage, loadImage, setImageProperties
export restrictResolutionToLevel
export getSize, getSpatialDomain, getPixelSpacing

include("./helpers/smoothing.jl")
export smoothArray




# load image and convert it into gray image
function loadImage(pathToFile;
                   voxelsize = [1.0,1.0], shift = [0.0, 0.0])

    # load image
    image = load(pathToFile)
    #image = permutedims(image, (2,1))

    image = convert(Array{ColorTypes.Gray},image)

    # convert image data to float (easier access)
    image = convert(Array{Float64,2},image)

    # set image properties

    return regImage(image, voxelsize, shift)

end

function createImage(imageData;
                     voxelsize = [1.0,1.0], shift = [0.0, 0.0])

    # load image
    #image = grayim(copy(imageData))

    return regImage(imageData, voxelsize, shift)

end

function restrictRegImage(im::regImage)
    newImage = deepcopy(im)
    newImage.data[:] = smoothArray(newImage.data, 3, 3.0)
    oldsize = size(im)
    newsize = map(Int,ceil.(oldsize./2))
    newGrid = getCellCenteredGrid(im.voxelsize*2, im.shift, newsize)
    newImage = createImage(interpolateImage(newImage,newGrid)[1], voxelsize = newGrid.voxelsize, shift = newGrid.shift)
    return newImage
end


function restrictResolutionToLevel(im::regImage,level)

    # copy restrict image, there seems to be a bug in the image copy method, take care of the deepcopy here
    #im = deepcopy(image.data)
    maxlevel = 0;
    for l=1:level
        if( (size(im.data,1)>2) && (size(im.data,2)>2) )
            im = restrictRegImage(im)
        end
    end
    return im
end

function getSize(I::regImage)
  return [size(I.data,1), size(I.data,2)]
end


end
