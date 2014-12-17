module ImageProcessing

using Images
export setImageProperties, loadImage, createImage, getPartOfImage, restrictResolutionToLevel, sizev

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

    # set pixel spacing
    # vertical direction:   pixelspacing(image)[1]
    # horizontal direction: pixelspacing(image)[2]
    image["pixelspacing"] = [(spatialDomain[2]-spatialDomain[1])/height(image),
                             (spatialDomain[4]-spatialDomain[3])/width(image)]

    return image

end

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

function getPartOfImage(image,verticalRange::UnitRange{Int64},horizontalRange::UnitRange{Int64})

    # define new image, copying properties
    partOfImage = similar(image)
    partOfImage.data = image[verticalRange,horizontalRange]

    # set new spatial domain
    spatialdomain = zeros(4)
    spatialdomain[1] = (verticalRange[1]-1)*pixelspacing(image)[1] + image["spatialdomain"][1]
    spatialdomain[2] =   verticalRange[end]*pixelspacing(image)[1] + image["spatialdomain"][1]
    spatialdomain[3] = (horizontalRange[1]-1)*pixelspacing(image)[2] + image["spatialdomain"][3]
    spatialdomain[4] =   horizontalRange[end]*pixelspacing(image)[2] + image["spatialdomain"][3]
    partOfImage["spatialdomain"] = spatialdomain

    # return new image
    return partOfImage

end

function restrictResolutionToLevel(image,level)

    # copy restrict image, there seems to ber a bug in the image copy method, take care of the deepcopy here
    restrictedImage = Image(copy(image.data), deepcopy(image.properties))
    maxlevel = 0;
    for l=1:level
        if( (width(restrictedImage)>2) && (height(restrictedImage)>2) )
            restrictedImage = restrict(restrictedImage)
            maxlevel = l
        end
    end

    # change pixel spacing accordingly
    spatialdomain = image["spatialdomain"]
    restrictedImage["pixelspacing"][1] = (spatialdomain[2]-spatialdomain[1]) / height(restrictedImage)
    restrictedImage["pixelspacing"][2] = (spatialdomain[4]-spatialdomain[3]) / width(restrictedImage)

    return restrictedImage

end

function sizev(I::Image)
  return [size(I.data,1), size(I.data,2) ]
end

end

