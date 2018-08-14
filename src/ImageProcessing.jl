module ImageProcessing

using Images
using DeformableRegistration.Types

export createImage, loadImage, setImageProperties
export restrictResolutionToLevel
export getSize, getSpatialDomain, getPixelSpacing

include("./helpers/smoothing.jl")
export smoothArray




# load image and convert it into gray image
function loadImage(pathToFile;
                   voxelsize = [1.0,1.0], shift = [0.0, 0.0])

    # load image
    image = Images.load(pathToFile)
    #image = permutedims(image, (2,1))

    # convert image to gray image
    image = Images.Gray.(image)

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


function restrictResolutionToLevel(image::regImage,level)

    # copy restrict image, there seems to be a bug in the image copy method, take care of the deepcopy here
    im = deepcopy(image.data)
    maxlevel = 0;
    for l=1:level
        if( (width(im)>2) && (height(im)>2) )
            im = restrict(im)
        end
    end

    voxelsizeNew = zeros(2)
    for xy = 1:2
        voxelsizeNew[xy] = (image.voxelsize[xy] * size(image.data,xy))/size(im,xy)
    end

    return regImage(im,voxelsizeNew,image.shift)

end

function getSize(I::regImage)
  return [height(I.data), width(I.data)]
end


end
