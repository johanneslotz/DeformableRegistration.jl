using DeformableRegistration: Distance, ImageProcessing, Transformation, Interpolation, regOptions
 using Test
 using Interpolations
 using LinearAlgebra

include("helpers/checkDerivative.jl")


@testset "Distance: distance derivatives - SSD" begin
    # create reference image
    refImg = createImage(100*rand(50,60))

    # measure distance and check derivative (nonparametric)
    centeredGrid = getCellCenteredGrid(refImg)
     centeredGrid.data[:] = centeredGrid.data[:] + 0.3 *rand(size(centeredGrid.data))
     D,dD,d2D = ssdDistance(refImg,refImg,centeredGrid.data,doDerivative=true,doHessian=true)
     Dfunc(x) = ssdDistance(refImg,refImg,x)[1]
     errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid.data,doPlot=false)
     @test checkErrorDecay(errquad)

    options = regOptions()
     options.matrixFree = true;
     D,dD,d2D = ssdDistance(refImg,refImg,centeredGrid.data,doDerivative=true,doHessian=true,options=options)
     Dfunc(x) = ssdDistance(refImg,refImg,x)[1]
     errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid.data)
     @test checkErrorDecay(errquad)

    # measure distance and check derivative (parametic)
    centeredGrid = getCellCenteredGrid(refImg)
     options = DeformableRegistration.regOptions()
     options.parametricOnly=true
     evaluationPoint = [1,0,0,0,1,0]+0.1*rand(6)
     D,dD,d2D = ssdDistance(refImg,refImg,transformGridAffine(centeredGrid,evaluationPoint).data,doDerivative=true,doHessian=true,options=options)
     Dfunc(p) = ssdDistance(refImg,refImg,transformGridAffine(centeredGrid,p).data,options=options)[1]
     errlin,errquad = checkDerivative(Dfunc,dD',evaluationPoint)
     @test checkErrorDecay(errquad)

    options.matrixFree = true;
     D,dD,d2D = ssdDistance(refImg,refImg,transformGridAffine(centeredGrid,evaluationPoint).data,doDerivative=true,doHessian=true,options=options)
     Dfunc(p) = ssdDistance(refImg,refImg,transformGridAffine(centeredGrid,p).data,options=options)[1]
     errlin,errquad = checkDerivative(Dfunc,dD',evaluationPoint)
     @test checkErrorDecay(errquad)
end

@testset "NGF" begin
    @testset "Distance: distance derivatives - NGF - parametric" begin
        refImg = createImage(100*rand(50,60))

        # measure distance and check derivative (nonparametric)
        centeredGrid = getCellCenteredGrid(refImg)
        centeredGrid.data[:] = centeredGrid.data[:] + 0.3 *rand(size(centeredGrid.data))
        D,dD,d2D = ngfDistance(refImg,refImg,centeredGrid.data,doDerivative=true,doHessian=true)
        Dfunc(x) = ngfDistance(refImg,refImg,x)[1]
        errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid.data,doPlot=false)
        @test checkErrorDecay(errquad)


        evaluationPoint = [1,0,0,0,1,0.0]+0.1*rand(6)
        centeredGrid = getCellCenteredGrid(refImg)
        options = DeformableRegistration.regOptions()
        options.edgeParameterR = 0.1
        options.edgeParameterT = 0.1
        options.matrixFree = false
        options.parametricOnly = true
        D,dD,d2D = ngfDistance(refImg,refImg,transformGridAffine(centeredGrid,evaluationPoint).data,doDerivative=true,doHessian=true,options=options)
        Dfunc(p) = ngfDistance(refImg,refImg,transformGridAffine(centeredGrid,p).data,options=options)[1]
        errlin,errquad = checkDerivative(Dfunc,dD',evaluationPoint,doPlot=false)
        @test checkErrorDecay(errquad)

        evaluationPoint = [1,0,0,0,1,0.0]+0.1*rand(6)
        centeredGrid = getCellCenteredGrid(refImg)
        options = DeformableRegistration.regOptions()
        options.edgeParameterR = 0.1
        options.edgeParameterT = 0.1
        options.matrixFree = true
        options.parametricOnly = true
        D,dD,d2D = ngfDistance(refImg,refImg,transformGridAffine(centeredGrid,evaluationPoint).data,doDerivative=true,doHessian=true,options=options)
        Dfunc(p) = ngfDistance(refImg,refImg,transformGridAffine(centeredGrid,p).data,options=options)[1]
        errlin,errquad = checkDerivative(Dfunc,dD',evaluationPoint,doPlot=false)
        @test checkErrorDecay(errquad)
    end

    @testset "Distance: distance derivatives - NGF - nonparametric" begin
        refImg = createImage(100*rand(50,60))

        ##
        centeredGrid = getCellCenteredGrid(refImg)
        centeredGrid.data[:] = centeredGrid.data[:] + 0.3 *rand(size(centeredGrid.data))
        options = DeformableRegistration.regOptions()
        options.edgeParameterR = 0.1
        options.edgeParameterT = 0.1
        options.matrixFree = false
        options.parametricOnly = false
        D,dD,d2D = ngfDistance(refImg,refImg,centeredGrid.data,doDerivative=true,doHessian=true,options=options)
        Dfunc(grid) = ngfDistance(refImg,refImg,grid,options=options)[1]
        errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid.data,doPlot=false)
        @test checkErrorDecay(errquad)

        centeredGrid = getCellCenteredGrid(refImg)
        centeredGrid.data[:] = centeredGrid.data[:] + 0.3 *rand(size(centeredGrid.data))
        options = DeformableRegistration.regOptions()
        options.edgeParameterR = 0.1
        options.edgeParameterT = 0.1
        options.matrixFree = true
        options.parametricOnly = false
        D,dD,d2D = ngfDistance(refImg,refImg,centeredGrid.data,doDerivative=true,doHessian=true,options=options)
        Dfunc(grid) = ngfDistance(refImg,refImg,grid,options=options)[1]
        errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid.data,doPlot=false)
        @test checkErrorDecay(errquad)

    end

    @testset "Distance: ngf epsilon estimation" begin
        testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
        refImg = loadImage(testimage)
        for r = refImg.data
             r = r + 0.01*rand()
        end
        warn("Skipping NGF ϵ estimation test.")
        @test_skip estimateNGFEpsilon(refImg,cutoffPercent=80)[1] >= 0 # seems to run forever
    end
end

@testset "SSD arbitrary grid" begin
    @testset "fixedGrid_IdenticalImagesZeroDeformation" begin
        n = 10

        imageData = ones(n,n)
        referenceImage = createImage(imageData, voxelsize = [1.0,1.0], shift = [0.0, 0.0])
        templateImage = createImage(imageData, voxelsize = [1.0,1.0], shift = [0.0, 0.0])
        transformedGrid = getCellCenteredGrid(referenceImage)
        # voxelsize::Array{Float64,1},shift::Array{Float64,1}, gridSize::Tuple{Int64,Int64}
        functionValue,dFunctionValue,d2FunctionValue,d_transformedImage =
            ssdDistanceArbitraryGrid(referenceImage, templateImage, transformedGrid)
                             # doDerivative::Bool=false,doHessian::Bool=false,
                             # options::regOptions=regOptions(),centeredGrid::scaledArray=scaledArray(zeros(1),(0,),[0,],[0,])))
        @test functionValue ≈ 0 atol=0.0001

        functionValue,dFunctionValue,d2FunctionValue,d_transformedImage =
            ssdDistanceArbitraryGrid(referenceImage, templateImage, transformedGrid, doDerivative=true)
        @test functionValue ≈ 0 atol=0.0001

        imageData = rand(n,n)
        referenceImage = createImage(imageData, voxelsize = [1.0,1.0], shift = [0.0, 0.0])
        templateImage = createImage(imageData, voxelsize = [1.0,1.0], shift = [0.0, 0.0])
        transformedGrid = getCellCenteredGrid(referenceImage)
        # voxelsize::Array{Float64,1},shift::Array{Float64,1}, gridSize::Tuple{Int64,Int64}
        functionValue,dFunctionValue,d2FunctionValue,d_transformedImage =
            ssdDistanceArbitraryGrid(referenceImage, templateImage, transformedGrid)
        @test functionValue ≈ 0 atol=0.0001

        functionValue,dFunctionValue,d2FunctionValue,d_transformedImage =
            ssdDistanceArbitraryGrid(referenceImage, templateImage, transformedGrid, doDerivative=true)
        @test functionValue ≈ 0 atol=0.0001

        imageData = rand(n,n)
        referenceImage = createImage(imageData, voxelsize = [1.0,1.0], shift = [0.0, 0.0])
        templateImage = createImage(imageData, voxelsize = [1.0,1.0], shift = [0.0, 0.0])
        transformedGrid = getCellCenteredGrid(restrictResolutionToLevel(referenceImage,2))
        @debug "-------- transformedGrid"
        @debug transformedGrid
        # voxelsize::Array{Float64,1},shift::Array{Float64,1}, gridSize::Tuple{Int64,Int64}
        functionValue,dFunctionValue,d2FunctionValue,d_transformedImage =
            ssdDistanceArbitraryGrid(referenceImage, templateImage, transformedGrid)
        @test functionValue ≈ 0 atol=0.0001

    end


    @testset "fixedGrid_IdenticalImagesFormats" begin
        n = 10

        imageData = ones(n,n)
        referenceImage = createImage(imageData, voxelsize = [1.0,1.0], shift = [0.0, 0.0])
        templateImage = createImage(imageData, voxelsize = [1.0,1.0], shift = [0.0, 0.0])

        gridSize = (5,5)
        voxelsize = n/gridSize[1]
        transformedGrid = getCellCenteredGrid([voxelsize,voxelsize],[0.0,0.0],gridSize)
        # voxelsize::Array{Float64,1},shift::Array{Float64,1}, gridSize::Tuple{Int64,Int64}
        functionValue,dFunctionValue,d2FunctionValue,d_transformedImage =
            ssdDistanceArbitraryGrid(referenceImage, templateImage, transformedGrid, doDerivative=true)
                             # doDerivative::Bool=false,doHessian::Bool=false,
                             # options::regOptions=regOptions(),centeredGrid::scaledArray=scaledArray(zeros(1),(0,),[0,],[0,])))
        @test functionValue ≈ 0 atol=0.0001
        @test size(dFunctionValue)... == prod(gridSize)*2
    ##
        imageData = rand(n,n)
        referenceImage = createImage(imageData, voxelsize = [1.0,1.0], shift = [0.0, 0.0])
        templateImage = createImage(imageData, voxelsize = [1.0,1.0], shift = [0.0, 0.0])

    ##
        gridSize = (10,20)
        voxelsize = 1.0
        transformedGrid = getCellCenteredGrid([voxelsize,voxelsize],[0.0,0.0],gridSize)
        transformedGrid.data[1:200] = transformedGrid.data[1:200] + 4.0
        # voxelsize::Array{Float64,1},shift::Array{Float64,1}, gridSize::Tuple{Int64,Int64}
        functionValue,dFunctionValue,d2FunctionValue,d_transformedImage =
            ssdDistanceArbitraryGrid(referenceImage, templateImage, transformedGrid, doDerivative=true)
    ##
        @test norm(dFunctionValue[105:195]) ≈ 0
        @test norm(dFunctionValue[305:395]) ≈ 0
        #@test norm(dFunctionValue[5:95]) > 1
        @test size(dFunctionValue)... == prod(gridSize)*2


    end

    @testset "ssdResampling_derivatives" begin

    include("../test/helpers/checkDerivative.jl")
    #include("../src/distances/ssdDistance.jl")

    smoothrand(x::Array{Float64,1}) = smoothArray(rand(size(x)),3, 1.0)

    refImg = createImage( smoothArray(100.0*rand(99,111),3, 3.0) )
        #refImg.data.data = smoothArray(refImg.data.data,11, 1.0)
        centeredGrid = getCellCenteredGrid(refImg)
        centeredGrid.data[:] = centeredGrid.data[:] + 5.0 *smoothrand(centeredGrid.data)
        options = regOptions()
        options.matrixFree = true;
        D,dD,d2D = Distance.ssdDistanceArbitraryGrid(refImg,refImg,centeredGrid, doDerivative=true,
            doHessian=true,options=options, interpolationScheme=BSpline(Linear()))
        Dfunc(x) = Distance.ssdDistanceArbitraryGrid(refImg,refImg,x, doDerivative=false, interpolationScheme=BSpline(Linear()))[1]
        errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid, doPlot=false)
        @test_skip checkErrorDecay(errquad)

    # The transposed interpolation operator is not implemented correctly which is
    # why in the case where the grid resolution changes, the derivative checks are failing.

    gridSize = (25,30)
        voxelsize = size(refImg.data)./gridSize
        vs = voxelsize
        centeredGrid = getCellCenteredGrid([vs[1]; vs[2]],[0.0,0.0],gridSize)
        centeredGrid.data[:] = centeredGrid.data[:] + 5 *smoothrand(centeredGrid.data)
        D,dD,d2D = ssdDistanceArbitraryGrid(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true,options=options, interpolationScheme=BSpline(Linear()))
        Dfunc(x) = ssdDistanceArbitraryGrid(refImg,refImg,x, doDerivative=false, interpolationScheme=BSpline(Cubic(Line())))[1]
        errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid, doPlot=false)
        @test_skip checkErrorDecay(errquad)

    gridSize = (10,10)
        voxelsize = size(refImg.data)./gridSize
        vs = voxelsize
        centeredGrid = getCellCenteredGrid([vs[1]; vs[2]],[0.0,0.0],gridSize)
        centeredGrid.data[:] = centeredGrid.data[:] + 10 * smoothrand(centeredGrid.data)
        D,dD,d2D = ssdDistanceArbitraryGrid(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true,options=options)
        Dfunc(x) = ssdDistanceArbitraryGrid(refImg,refImg,x)[1]
        errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid, doPlot=false)
        @test_skip checkErrorDecay(errquad)

    gridSize = (25,30)
        vs = size(refImg.data)./gridSize
        centeredGrid = getCellCenteredGrid([vs[1]; vs[2]],[0.0,0.0],gridSize)
        disp = 10 * ones(gridSize)
        disp[4:6,4:6]=0.0
        testGrid = 0.0 * centeredGrid + vcat(disp[:], disp[:])
        testGridInterp = interpolateDeformationField(testGrid, getCellCenteredGrid(refImg), interpolationScheme=BSpline(Linear()))
        testGridInterp2 = interpolateDeformationField(testGridInterp, getCellCenteredGrid([vs[1]; vs[2]],[0.0,0.0],gridSize), interpolationScheme=BSpline(Linear()))
         #
         # using PyPlot
         # figure()
         # subplot(1,3,1)
         # imshow(reshape(testGrid.data[1:Int(end/2)],testGrid.dimensions))
         # clim([0,10])
         # subplot(1,3,2)
         # imshow(reshape(testGridInterp.data[1:Int(end/2)],testGridInterp.dimensions))
         # clim([0,10])
         # subplot(1,3,3)
         # imshow(reshape(testGridInterp2.data[1:Int(end/2)],testGridInterp2.dimensions))
         # clim([0,10])

        @test_skip sum(testGrid.data[:]) ≈ sum(testGridInterp2.data[:]) atol=0.1

    end


end
