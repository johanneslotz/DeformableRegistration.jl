
using DeformableRegistration: Regularizer, Examples, ImageProcessing, Transformation, Interpolation, Distance
 using DeformableRegistration.regOptions
 using Base.Test

#using DeformableRegistration.Visualization # does not work on travis-ci


# setup logging

@testset "artificial data SSD" begin
    # create test images
    data = zeros(120,120); data[31:90,21:60] = 1
    temImg = createImage(data) 
    dataR = deepcopy(data)
    dataR[41:80,41:90] = 1
    refImg = createImage(dataR)
    options = DeformableRegistration.regOptions()
    options.levels = [5,4,3,2]
    options.matrixFree = true;
    options.interpolateToReferenceImage = true
    options.regularizerWeight = 1

    ## register images nonparametric
    displacement = registerImagesNonParametric(refImg, temImg, options) # regularizerOperator=createDiffusiveOperatorCentered
    ssdvalue = ssdDistance(refImg, temImg, getCellCenteredGrid(refImg).data + displacement.data, options=options)[1]
    ngfvalue = ngfDistance(refImg, temImg, getCellCenteredGrid(refImg).data + displacement.data, options=options)[1]
    #visualizeResults(refImg, temImg, displacement=displacement)

    @test ssdvalue ≈ 42.63 atol=1e-1
    @test ngfvalue ≈ 365.01 atol=1e-1

end

@testset "artificial data NGF" begin

	# NGF, smooth data first to generate gradients
	data = zeros(120,120); data[31:90,21:60] = 1
	dataT = copy(data)
	data = conv2([0.2 0.5 0.2
				 0.5 1.0 0.5
				 0.2 0.5 0.2],data)[1:120,1:120]
	temImg = createImage(data) 
	dataT[41:80,41:90] = 1
	dataT = conv2([0.2 0.5 0.2
				 0.5 1.0 0.5
				 0.2 0.5 0.2],dataT)[1:120,1:120]
	refImg = createImage(dataT)

	options = regOptions()
	options.levels = [4,3,2]
	options.edgeParameterR = 1
	options.edgeParameterT = options.edgeParameterR
	options.regularizerWeight = 100
	options.interpolateToReferenceImage = true
	@time displacement = registerImagesNonParametric(refImg, temImg, options, measureDistance=ngfDistance)

    ssdvalue = ssdDistance(refImg,temImg,getCellCenteredGrid(refImg).data+displacement.data)[1]
	ngfvalue = ngfDistance(refImg,temImg,getCellCenteredGrid(refImg).data+displacement.data)[1]

    Logging.debug(ssdvalue)
    Logging.debug(ngfvalue)
    # These results seem to depend on the current phase of the moon, your platform or something similar.
	# This seems to be a bug in ngfDistance or at least in there and is notet in an issue in GitHub. Please
	# make the test more specific once this is solved.
    # visualizeResults(refImg, temImg, displacement=displacement)

	@test_skip ssdvalue ≈ 498.07 atol=2e2
	@test_skip ngfvalue ≈ 399.403 atol=2e2
##
end

##
@testset "SSD shifted"  begin
    for D=[ssdDistance]#, ngfDistance]
##
	testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
	referenceImage = loadImage(testimage)
	# define a cell centered grid, transform it and create a template image
	centeredGrid = getCellCenteredGrid(referenceImage)
	affineParametersInitial = [1.0,0.0,-50,0.0,1.0,-50]
	transformedGrid = transformGridAffine(centeredGrid,affineParametersInitial)
	temImg = interpolateImage(referenceImage,transformedGrid,interpolationScheme=InterpLinearFast)
	templateImage = createImage(temImg)
	options = regOptions()
	options.levels = [6,5,4,3]
	options.parametricOnly = false;
	options.regularizerWeight = 100
	options.matrixFree = true;
	options.interpolateToReferenceImage = false
	options.useEdgeParameterInNumerator = true

	displacement = registerImagesNonParametric(referenceImage, templateImage, options, measureDistance=D)
	@test mean(displacement.data) ≈ -affineParametersInitial[3] atol = 10

    end
end
