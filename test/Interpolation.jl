using DeformableRegistration: Transformation, ImageProcessing, Interpolation
 using Base.Test
 using Interpolations: BSpline, Linear, Cubic, Free
 using Logging
 Logging.configure(level=INFO)

@testset "Interpolation" begin
##
@testset "Interpolation: interpolation" begin
    ## load reference image

    testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
    img = loadImage(testimage)
    centeredGrid = getCellCenteredGrid(img)
    affineParameters = [0.5,0.4,50,0,0.5,100]
    deformationField = zeros(prod(size(img.data))*2)
    deformationField[prod(size(img.data))+1:end] = 10*sin.(0.01*centeredGrid.data[1:prod(size(img.data))])
    deformationField = scaledArray(deformationField, size(img.data), [1, 1.0], [0.0, 0])
    transformedGrid = transformGridAffine(centeredGrid,affineParameters) + deformationField
    transformedImageOnly,dx,dy = interpolateImage(img,transformedGrid,doDerivative=true)
    transformedImage,dY_transformedImage,dX_transformedImage = interpolateImage(img,transformedGrid,doDerivative=true)
    transformedImageRef,dY_transformedImageRef,dX_transformedImageRef = interpolateImage(img,transformedGrid,doDerivative=true,interpolationScheme=BSpline(Linear()))

    ## using PyPlot
    #dY_transformedImage = reshape(dY_transformedImage, size(img.data))
    #dY_transformedImageRef = reshape(dY_transformedImageRef, size(img.data))
    #imshow(dY_transformedImage-dY_transformedImageRef)
    @test norm(transformedImageOnly - transformedImageRef) < 1.0
    @test transformedImageOnly[1:10] ≈ transformedImageRef[1:10]

    @test norm(transformedImage - transformedImageRef) < 1.0
    @test transformedImage[1:10] ≈ transformedImageRef[1:10]

    @test norm(dY_transformedImage - dY_transformedImageRef) ≈ 0 atol=20
    @test dY_transformedImage[1:10] ≈ dY_transformedImageRef[1:10]

    @test norm(dX_transformedImage - dX_transformedImageRef) ≈ 0 atol=20
    @test dX_transformedImage[1:10] ≈ dX_transformedImageRef[1:10]
end

##
@testset "Interpolation: check derivatives" begin
    include("helpers/checkDerivative.jl")
    testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
    img = loadImage(testimage)

    # check image derivative
    centeredGrid = getCellCenteredGrid(img).data + 0.1
    Ifunc(x) = interpolateImage(img,x)[:]
    transformedImage, dX_transformedImage, dY_transformedImage = interpolateImage(img,centeredGrid,doDerivative=true)
    dTransformedImage = spdiagm((dX_transformedImage[:],dY_transformedImage[:]),[0, prod(size(img.data))])
    errlin,errquad = checkDerivative(Ifunc,dTransformedImage,centeredGrid)
    @test checkErrorDecay(errquad)
end
##
@testset "Interpolation: check timing of linear interpolation" begin
    using Interpolations: BSpline, Linear
    # check timing
    img = createImage(rand(1024,1024))
    centeredGrid = getCellCenteredGrid(img)
    affineParameters = [0.5,0.4,50,0,0.5,100]
    deformationField = zeros(prod(size(img.data))*2)
    deformationField[prod(size(img.data))+1:end] = 10*sin.(0.01*centeredGrid.data[1:prod(size(img.data))])
    deformationField = scaledArray(deformationField, size(img.data), [1, 1.0], [0.0, 0])
    transformedGrid = transformGridAffine(centeredGrid,affineParameters) + deformationField
    # InterpLinearFast
    tic();
    for i=1:5
      interpolateImage(img,transformedGrid)
    end
    timing = toq()/5;
    Logging.info("Interpolation: interpolateImage (InterpLinearFast) took ",timing," seconds.")
    # InterpLinear (Grid.jl)
    tic();
    for i=1:5
      interpolateImage(img,transformedGrid,interpolationScheme=BSpline(Linear()))
    end
    timingGrid = toq()/5;
    Logging.info("Interpolation: interpolateImage (InterpLinear, Interpolations.jl) took ",timingGrid," seconds.")
    # InterpLinearFast with derivative
    tic();
    for i=1:5
      interpolateImage(img,transformedGrid,doDerivative=true)
    end
    timing = toq()/5;
    Logging.info("Interpolation: interpolateImage with derivative (InterpLinearFast) took ",timing," seconds.")
    # InterpLinear (Grid.jl) with derivative
    tic();
    for i=1:5
      interpolateImage(img,transformedGrid,interpolationScheme=BSpline(Linear()),doDerivative=true)
    end
    timingGrid = toq()/5;
    Logging.info("Interpolation: interpolateImage with derivative (InterpLinear, Interpolations.jl) took ",timingGrid," seconds.")
    @test timingGrid>timing
    #Logging.info("Interpolation: Own linear interpolation with derivative is faster ✔")
end
##
@testset "compare different interpolation types on real data" begin
    testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
    refImg = loadImage(testimage)
    # define a cell centered grid, transform it and create a template image
    centeredGrid = getCellCenteredGrid(refImg)
    affineParametersInitial = [1.2,0.0,-50,0.0,1.2,-0]
    transformedGrid = transformGridAffine(centeredGrid,affineParametersInitial)
    temImg = interpolateImage(refImg,transformedGrid,interpolationScheme=InterpLinearFast)
    temImg = createImage(temImg)
    ##

    temImg2 = interpolateImage(refImg,transformedGrid,interpolationScheme=BSpline(Cubic(Free())))
    temImg2 = createImage(temImg2)


    if Logging.LogLevel == Logging.DEBUG
    using PyPlot
    figure(); PyPlot.imshow(Array(temImg2.data)-temImg.data)
    end

    @test norm(Array(temImg2.data)[:]-temImg.data[:])  < 20
end

@testset "compare different interpolation types on real data" begin
    testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
    refImg = loadImage(testimage)
    # define a cell centered grid, transform it and create a template image
    centeredGrid = getCellCenteredGrid(refImg)
    affineParametersInitial = [1.2,0.0,-50,0.0,1.2,-0]
    transformedGrid = transformGridAffine(centeredGrid,affineParametersInitial)
    temImg = interpolateImage(refImg,transformedGrid,interpolationScheme=InterpLinearFast)
    temImg = createImage(temImg)

    temImg2 = interpolateImage(refImg,transformedGrid,interpolationScheme=BSpline(Linear()))
    temImg2 = createImage(temImg2)

    if Logging.LogLevel == Logging.DEBUG
    #using PyPlot
    #figure(); PyPlot.imshow(Array(temImg2.data)-temImg.data)
    end
    @test norm(Array(temImg2.data)[:]-temImg.data[:])  < 10
end

##
@testset "interpolation with identity grid results in orgiginal data: $scheme" for scheme in (InterpLinearFast, BSpline(Linear()), BSpline(Cubic(Free())))
    data = [1.0 2; 3 4]
    testimage = createImage(data)
    g = getCellCenteredGrid(testimage)
    dataInterpolated = interpolateImage(testimage, g, interpolationScheme=scheme)
    @test data ≈ dataInterpolated
end


end
