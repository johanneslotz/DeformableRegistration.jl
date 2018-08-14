using DeformableRegistration: Transformation, ImageProcessing, Interpolation, Distance
 using Base.Test
 using Interpolations: BSpline, Linear, Cubic, Line, Free
#  using Logging
#  Logging.configure(level=INFO)


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

@testset "Interpolation: interpolation with target grid" begin
    ## load reference image

    testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
    img = loadImage(testimage)
    img = restrictResolutionToLevel(img,1)
    centeredGrid = getCellCenteredGrid(img)
    affineParameters = [0.5,0.4,50,0,0.5,100]
    deformationField = zeros(prod(size(img.data))*2)
    deformationField[prod(size(img.data))+1:end] = 10*sin.(0.01*centeredGrid.data[1:prod(size(img.data))])
    deformationField = scaledArray(deformationField, size(img.data), img.voxelsize, img.shift)
    transformedGrid = transformGridAffine(centeredGrid,affineParameters) + deformationField

    targetGrid = getCellCenteredGrid(deformationField)


    transformedImage,dx,dy = interpolateImage(img,transformedGrid,doDerivative=true, interpolationScheme=BSpline(Linear()))
    tic()
    transformedImageTG,dxTG,dyTG = interpolateImage(img,transformedGrid,targetGrid, doDerivative=true, interpolationScheme=BSpline(Linear()))
    timing = toc()

    # @test timing < 0.3 not working on CI server :(

    @test norm(transformedImage - transformedImageTG) < 1e-10
    @test norm(dx - dxTG) < 1e-10
    @test norm(dy - dyTG) < 1e-10

    newSize = size(img.data)./2
    newSize = (Int(ceil(newSize[1])), Int(ceil(newSize[2])))
    targetGrid = getCellCenteredGrid(img.voxelsize*2, img.shift, newSize)
    @time transformedImageTG,dxTG,dyTG = interpolateImage(img,transformedGrid,targetGrid, doDerivative=true)

    transformedImageCoarse = transformedImage[1:2:end,1:2:end]
    @test norm(transformedImageCoarse[:] - transformedImageTG[:]) < 4
end


@testset "Interpolation: interpolation of part of image with target grid" begin
    ## load reference image

    testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
    img = loadImage(testimage)
    img = restrictResolutionToLevel(img,1)
    centeredGrid = getCellCenteredGrid(img.voxelsize*3, img.shift, (200,200))
    affineParameters = [0.5,0.4,25,0,0.5,50]
    transformedGrid = transformGridAffine(centeredGrid,affineParameters)

    # using PyPlot


    @time transformedImage,dx,dy = interpolateImage(img,transformedGrid,doDerivative=true, interpolationScheme=BSpline(Linear()))
    # subplot(221)
    # imshow(transformedImage)
    # subplot(222)
    # imshow(reshape(dx[1:200*200],transformedGrid.dimensions))

    targetGrid = getCellCenteredGrid(img.voxelsize*6, img.shift, (100,100))
    @time transformedImageTG,dxTG,dyTG = interpolateImage(img,transformedGrid,targetGrid, doDerivative=true)
    # subplot(223)
    # imshow(transformedImageTG)
    # subplot(224)
    # imshow(reshape(dxTG[1:100*100],targetGrid.dimensions))
    @test norm(transformedImage[1:2:200,1:2:200] - transformedImageTG) < 1

    targetGrid = getCellCenteredGrid(img.voxelsize*3, img.shift, (100,100))
    @time transformedImageTG,dxTG,dyTG = interpolateImage(img,transformedGrid,targetGrid, doDerivative=true)
    @test norm(transformedImage[1:100,1:100] - transformedImageTG) < 1


end

@testset "Interpolation: check derivatives with target grid 1" begin
    for interpolationScheme=[BSpline(Linear()), BSpline(Cubic(Line()))]

        include("helpers/checkDerivative.jl")
        testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
        img = loadImage(testimage)

        targetGrid = getCellCenteredGrid(img)
        # check image derivative
        transformedGrid = getCellCenteredGrid(img)
        transformedGrid.data[:] = transformedGrid.data[:] + 0.1
        Ifunc(x) = interpolateImage(img, x, targetGrid, interpolationScheme=interpolationScheme)[1][:]
        transformedImage, dX_transformedImage, dY_transformedImage =
            interpolateImage(img,transformedGrid,targetGrid,doDerivative=true,
                interpolationScheme=interpolationScheme)
        dTransformedImage = spdiagm((dX_transformedImage[:],dY_transformedImage[:]),[0, prod(size(img.data))])
        @time errlin,errquad = checkDerivative(Ifunc,dTransformedImage,transformedGrid)
        @test checkErrorDecay(errquad)
    end
end


##
@testset "Interpolation: check derivatives (all interpolation types)" for
    interpolationScheme=[InterpLinearFast, BSpline(Linear()), BSpline(Cubic(Line()))]
    include("helpers/checkDerivative.jl")
    testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
    img = loadImage(testimage)

    # check image derivative
    centeredGrid = getCellCenteredGrid(img).data + 0.1
    Ifunc(x) = interpolateImage(img,x, interpolationScheme=interpolationScheme)[1][:]
    @time transformedImage, dX_transformedImage, dY_transformedImage =
        interpolateImage(img,centeredGrid,doDerivative=true,
            interpolationScheme=interpolationScheme)
    dTransformedImage = spdiagm((dX_transformedImage[:],dY_transformedImage[:]),[0, prod(size(img.data))])
    @time errlin,errquad = checkDerivative(Ifunc,dTransformedImage,centeredGrid)
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
    info("Interpolation: interpolateImage (InterpLinearFast) took ",timing," seconds.")
    # InterpLinear (Grid.jl)
    tic();
    for i=1:5
      interpolateImage(img,transformedGrid,interpolationScheme=BSpline(Linear()))
    end
    timingGrid = toq()/5;
    info("Interpolation: interpolateImage (InterpLinear, Interpolations.jl) took ",timingGrid," seconds.")
    # InterpLinearFast with derivative
    tic();
    for i=1:5
      interpolateImage(img,transformedGrid,doDerivative=true)
    end
    timing = toq()/5;
    info("Interpolation: interpolateImage with derivative (InterpLinearFast) took ",timing," seconds.")
    # InterpLinear (Grid.jl) with derivative
    tic();
    for i=1:5
      interpolateImage(img,transformedGrid,interpolationScheme=BSpline(Linear()),doDerivative=true)
    end
    timingGrid = toq()/5;
    info("Interpolation: interpolateImage with derivative (InterpLinear, Interpolations.jl) took ",timingGrid," seconds.")
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
    temImg = interpolateImage(refImg,transformedGrid,interpolationScheme=InterpLinearFast)[1]
    temImg = createImage(temImg)
    ##
    targetGrid = getCellCenteredGrid(refImg)
    targetGrid.data[:] = targetGrid.data[:] - 0
    temImg2 = interpolateImage(refImg,transformedGrid,targetGrid,interpolationScheme=BSpline(Cubic(Free())))[1]
    temImg2 = createImage(temImg2)


    # if Logging.LogLevel == Logging.DEBUG
    # using PyPlot
    # figure(); PyPlot.imshow(Array(temImg2.data)-temImg.data)
    # end

    # border pixels differe in both interpolation schemes, most likely due
    # to InterpLinearFast padding the image too far...
    # allowing 2 pixels around the border of the image to be different
    @test norm(Array(temImg2.data)[:]-temImg.data[:])  < 32
    nDifferentPixels = sum(abs.(Array(temImg2.data)[:]-temImg.data[:]) .> 0.01)
    @test nDifferentPixels < maximum(size(refImg.data))*8
end

@testset "compare different interpolation types on real data" begin
    testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
    refImg = loadImage(testimage)
    # define a cell centered grid, transform it and create a template image
    centeredGrid = getCellCenteredGrid(refImg)
    affineParametersInitial = [1.2,0.0,-50,0.0,1.2,-0]
    transformedGrid = transformGridAffine(centeredGrid,affineParametersInitial)
    temImg = interpolateImage(refImg,transformedGrid,interpolationScheme=InterpLinearFast)[1]
    temImg = createImage(temImg)

    temImg2 = interpolateImage(refImg,transformedGrid,interpolationScheme=BSpline(Linear()))[1]
    temImg2 = createImage(temImg2)

    # if Logging.LogLevel == Logging.DEBUG
    # #using PyPlot
    # #figure(); PyPlot.imshow(Array(temImg2.data)-temImg.data)
    # end

    # border pixels differe in both interpolation schemes, most likely due
    # to InterpLinearFast padding the image too far...
    # allowing 2 pixels around the border of the image to be different
    @test norm(Array(temImg2.data)[:]-temImg.data[:])  < 32
    nDifferentPixels = sum(abs.(Array(temImg2.data)[:]-temImg.data[:]) .> 0.01)
    @test nDifferentPixels < maximum(size(refImg.data))*4
end

##
@testset "interpolation with identity grid results in orgiginal data: $scheme" for scheme in (InterpLinearFast, BSpline(Linear()), BSpline(Cubic(Free())))
    data = [1.0 2; 3 4]
    testimage = createImage(data)
    g = getCellCenteredGrid(testimage)
    dataInterpolated = interpolateImage(testimage, g, interpolationScheme=scheme)[1]
    @test data ≈ dataInterpolated
end

@testset "smoothing" begin
    for k=[3,5,7,11]
        k_half = Int(ceil(k/2))
        for x = [ones(10,11)*2.0, rand(63,26)]
            x[1:k_half,:] = 0
            x[:,1:k_half] = 0
            x[end-k_half:end,:] = 0
            x[:,end-k_half:end] = 0

            #x = ones(10,11)*2.0;
            xnew = smoothArray(x,k, Float64(k))
            @test sum(x[:]) ≈ sum(xnew[:])
        end
    end

end


end
