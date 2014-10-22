using ImageRegistration.ImageProcessing
using Images
using Base.Test

function checkImageProperties(img::Image)
    # check if image has two dimensions
    assert2d(img)
    # check if image values are between 0 and 1
    @test float64(minfinite(img))>=0.0
    @test float64(maxfinite(img))<=1.0
    # check if colorspace is gray
    @test colorspace(img) == "Gray"
end

# test loadImage
testimage = dirname(Base.source_path()) * "/testdata/testimage.png"
img = loadImage(testimage)
checkImageProperties(img)
# check if spatial domain is the one that was specified
img = loadImage(testimage,spatialDomain = [-1.1,2.1,-1.0,2.21])
@test img.properties["spatialdomain"] == [-1.1,2.1,-1.0,2.21]

# test createImage
imgdata = rand(43,123)
img = createImage(imgdata)
checkImageProperties(img)
@test size(imgdata) == (43,123)
# check if spatial domain is the one that was specified
img = createImage(imgdata,spatialDomain = [-1.1,2.1,-1.0,2.21])
@test img.properties["spatialdomain"] == [-1.1,2.1,-1.0,2.21]

# test setImageProperties
img = createImage(imgdata)
img = setImageProperties(img)
checkImageProperties(img)

# test getPartOfImage
img = loadImage(testimage)
imgpart = getPartOfImage(img,1:161,1:401)
# check if image part has defined size
@test size(imgpart) == (161,401)
# check if spatialdomain is correct
@test (imgpart.properties["spatialdomain"][2]-img.properties["spatialdomain"][1])/height(imgpart) == pixelspacing(imgpart)[1]
@test (imgpart.properties["spatialdomain"][4]-img.properties["spatialdomain"][3])/width(imgpart) == pixelspacing(imgpart)[2]
# check if it returns image data
imgpart = getPartOfImage(img,1:height(img),1:width(img))
@test imgpart == img

# test restrictResolutionToLevel
img = loadImage(testimage)
# check image properties for each level of the created image
for i=1:10
  restrictedImage = restrictResolutionToLevel(img,i)
  checkImageProperties(restrictedImage)
  @test pixelspacing(img) != pixelspacing(restrictedImage)
end

# load image, write image, load the written image and compare it
testimagewrite = dirname(Base.source_path()) * "/testdata/testimage_write.png"
#show(img)
imwrite(img, testimagewrite)
imgwrite = loadImage(testimagewrite)
@test_approx_eq_eps(img,imgwrite,1e-3)
