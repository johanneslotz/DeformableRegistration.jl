using Images
 using Logging
 using DeformableRegistration: ImageProcessing, Transformation, Interpolation, Optimization
 using DeformableRegistration.Regularizer
 using DeformableRegistration.regOptions
 using Base.Test
 using DeformableRegistration.Optimization.augmentedLagrangian
 Logging.configure(level=Logging.INFO)
##

include("./fixedGridRegistration.jl")


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
    @test functionValue == 0


    imageData = rand(n,n)
    referenceImage = createImage(imageData, voxelsize = [1.0,1.0], shift = [0.0, 0.0])
    templateImage = createImage(imageData, voxelsize = [1.0,1.0], shift = [0.0, 0.0])
    transformedGrid = getCellCenteredGrid(referenceImage)
    # voxelsize::Array{Float64,1},shift::Array{Float64,1}, gridSize::Tuple{Int64,Int64}
    functionValue,dFunctionValue,d2FunctionValue,d_transformedImage =
        ssdDistanceArbitraryGrid(referenceImage, templateImage, transformedGrid)
    @test functionValue == 0

    imageData = rand(n,n)
    referenceImage = createImage(imageData, voxelsize = [1.0,1.0], shift = [0.0, 0.0])
    templateImage = createImage(imageData, voxelsize = [1.0,1.0], shift = [0.0, 0.0])
    transformedGrid = getCellCenteredGrid(restrictResolutionToLevel(referenceImage,2))
    debug("-------- transformedGrid")
    debug(transformedGrid)
    # voxelsize::Array{Float64,1},shift::Array{Float64,1}, gridSize::Tuple{Int64,Int64}
    functionValue,dFunctionValue,d2FunctionValue,d_transformedImage =
        ssdDistanceArbitraryGrid(referenceImage, templateImage, transformedGrid)
    @test functionValue == 0

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
    @test functionValue == 0
    @test size(dFunctionValue)... == prod(gridSize)*2

end


@testset "smoothing" begin
    for k=[3,5,7,11]
        for x = [ones(10,11)*2, rand(63,26)]
            x = ones(10,11)*2;
            xnew = smoothArray(x,k, Float64(k))
            @test sum(x[:])  == sum(xnew[:])
        end
    end

end

##
@testset "ssdResampling_derivatives" begin

include("../test/helpers/checkDerivative.jl")

smoothrand(x::Array{Float64,1}) = smoothArray(rand(size(x)),5, 5.0)

refImg = createImage(100*rand(50,60))
    refImg.data.data = smoothArray(refImg.data.data,11, 11.0)
    centeredGrid = getCellCenteredGrid(refImg)
    centeredGrid.data[:] = centeredGrid.data[:] + 10 *smoothrand(centeredGrid.data)
    options = regOptions()
    options.matrixFree = true;
    D,dD,d2D = ssdDistanceArbitraryGrid(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true,options=options)
    Dfunc(x) = ssdDistanceArbitraryGrid(refImg,refImg,x)[1]
    errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid)
    @test checkErrorDecay(errquad)


gridSize = (90,60)
    voxelsize = size(refImg.data)./gridSize
    vs = voxelsize
    centeredGrid = getCellCenteredGrid([vs[1]; vs[2]],[0.0,0.0],gridSize)
    centeredGrid.data[:] = centeredGrid.data[:] + 10 *smoothrand(centeredGrid.data)
    D,dD,d2D = ssdDistanceArbitraryGrid(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true,options=options)
    Dfunc(x) = ssdDistanceArbitraryGrid(refImg,refImg,x)[1]
    errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid)
    #@test
    checkErrorDecay(errquad)

gridSize = (25,30)
    voxelsize = size(refImg.data)./gridSize
    vs = voxelsize
    centeredGrid = getCellCenteredGrid([vs[1]; vs[2]],[0.0,0.0],gridSize)
    centeredGrid.data[:] = centeredGrid.data[:] + 10 * smoothrand(centeredGrid.data)
    D,dD,d2D = ssdDistanceArbitraryGrid(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true,options=options)
    Dfunc(x) = ssdDistanceArbitraryGrid(refImg,refImg,x)[1]
    errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid)
    @test checkErrorDecay(errquad)

gridSize = (25,30)
    vs = size(refImg.data)./gridSize
    centeredGrid = getCellCenteredGrid([vs[1]; vs[2]],[0.0,0.0],gridSize)
    disp = 10 * ones(gridSize)
    disp[4:6,4:6]=0.0
    testGrid = 0.0 * centeredGrid + vcat(disp[:], disp[:])
    testGridInterp = interpolateDeformationField(testGrid, getCellCenteredGrid(refImg), interpolationScheme=BSpline(Constant()))
    testGridInterp2 = interpolateDeformationField(testGridInterp, getCellCenteredGrid([vs[1]; vs[2]],[0.0,0.0],gridSize), interpolationScheme=BSpline(Linear()))

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

    @test sum(testGrid.data[:]) ≈ sum(testGridInterp2.data[:]) atol=0.1


gridSize = (45,45)
    voxelsize = size(refImg.data)./gridSize[1]
    vs = voxelsize
    centeredGrid = getCellCenteredGrid([vs[1]; vs[2]],[0.0,0.0],gridSize)
    x= 0.3 * smoothrand(centeredGrid.data)
    centeredGrid.data[:] = centeredGrid.data[:] + x
    @profile D,dD,d2D = ssdDistanceArbitraryGrid(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true,options=options)
    D,dD,d2D = ssdDistanceArbitraryGrid(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true,options=options)
    Dfunc(x) = ssdDistanceArbitraryGrid(refImg,refImg,x)[1]
    @time errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid, smoothOffset=true )
    @test checkErrorDecay(errquad)

gridSize = (10,10)
    voxelsize = size(refImg.data)./gridSize[1]
    vs = voxelsize
    centeredGrid = getCellCenteredGrid([vs[1]; vs[2]],[0.0,0.0],gridSize)
    centeredGrid.data[:] = centeredGrid.data[:] + 1 *smoothrand(centeredGrid.data)
    D,dD,d2D = ssdDistanceArbitraryGrid(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true,options=options)
    Dfunc(x) = ssdDistanceArbitraryGrid(refImg,refImg,x)[1]
    errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid)
    @test checkErrorDecay(errquad)

end
