module ImageProcessing

using Images

export createImage, loadImage, setImageProperties
export restrictResolutionToLevel
export getSize, getSpatialDomain, getPixelSpacing

# load image and convert it into gray image
function loadImage(pathToFile;
                   spatialDomain = [0.0,0.0,0.0,0.0])

    # load image
    image = Images.imread(pathToFile)
    image.data = image.data'

    # convert image to gray image
    image = convert(Image{Gray}, image)

    # convert image data to float (easier access)
    image = grayim(convert(Array{Float64,2},image.data))

    # set image properties
    setImageProperties(image,spatialDomain = spatialDomain)

    return image

end

function createImage(imageData;
                     spatialDomain = [0.0,0.0,0.0,0.0])

    # load image
    image = grayim(copy(imageData))

    # set image properties
    setImageProperties(image,spatialDomain = spatialDomain)

    return image

end

function setImageProperties(image;
                            spatialDomain = [0.0,0.0,0.0,0.0])

    # define spatial order
    # y: vertical axis
    # x: horizontal axis
    image["spatialorder"] = ["y","x"]

    # set spatial domain
    # vertical direction:   spatialDomain[1] to spatialDomain[2]
    # horizontal direction: spatialDomain[3] to spatialDomain[4]
    if(spatialDomain==[0.0,0.0,0.0,0.0])
        spatialDomain[2] = size(image)[1]
        spatialDomain[4] = size(image)[2]
    end
    image["spatialdomain"] = spatialDomain

    return image

end

function restrictResolutionToLevel(image,level)

    # copy restrict image, there seems to be a bug in the image copy method, take care of the deepcopy here
    restrictedImage = Image(copy(image.data), deepcopy(image.properties))
    maxlevel = 0;
    for l=1:level
        if( (width(restrictedImage)>2) && (height(restrictedImage)>2) )
            restrictedImage = restrict(restrictedImage)
            maxlevel = l
        end
    end

    return restrictedImage

end

function getSize(I::Image)
  return [height(I), width(I)]
end

function getSpatialDomain(I::Image)
    return I["spatialdomain"]
end

function getPixelSpacing(I::Image)
    # vertical direction:   getPixelSpacing(image)[1]
    # horizontal direction: getPixelSpacing(image)[2]
    return [(I["spatialdomain"][2]-I["spatialdomain"][1])/height(I),
            (I["spatialdomain"][4]-I["spatialdomain"][3])/width(I)]
end

end

