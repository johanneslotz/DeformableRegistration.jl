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
    for x = [ones(10,11)*2, rand(63,26)]
        x = ones(10,11)*2;
        xnew = smoothArray(x)
        @test sum(x[:])  == sum(xnew[:])
    end

end

@testset "ssdResampling_derivatives" begin

include("../test/helpers/checkDerivative.jl")


refImg = createImage(100*rand(50,60))
 refImg.data.data = smoothArray(refImg.data.data,11)
# measure distance and check derivative (nonparametric)
centeredGrid = getCellCenteredGrid(refImg)
 centeredGrid.data[:] = centeredGrid.data[:] + 0.3 *rand(size(centeredGrid.data))
 options = regOptions()
 options.matrixFree = true;
 D,dD,d2D = ssdDistanceArbitraryGrid(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true,options=options)
 Dfunc(x) = ssdDistanceArbitraryGrid(refImg,refImg,x)[1]
 errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid)
 @test checkErrorDecay(errquad)

gridSize = (30,30)
 voxelsize = size(refImg.data)./gridSize[1]
 vs = voxelsize
 centeredGrid = getCellCenteredGrid([vs[1]; vs[2]],[0.0,0.0],gridSize)
 centeredGrid.data[:] = centeredGrid.data[:] + 0.3 *rand(size(centeredGrid.data))
  D,dD,d2D = ssdDistanceArbitraryGrid(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true,options=options)
  Dfunc(x) = ssdDistanceArbitraryGrid(refImg,refImg,x)[1]
  errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid)
  @test checkErrorDecay(errquad)

gridSize = (10,10)
  voxelsize = size(refImg.data)./gridSize[1]
  vs = voxelsize
  centeredGrid = getCellCenteredGrid([vs[1]; vs[2]],[0.0,0.0],gridSize)
  centeredGrid.data[:] = centeredGrid.data[:] + 0.1/vs[1] *rand(size(centeredGrid.data))
   D,dD,d2D = ssdDistanceArbitraryGrid(refImg,refImg,centeredGrid,doDerivative=true,doHessian=true,options=options)
   Dfunc(x) = ssdDistanceArbitraryGrid(refImg,refImg,x)[1]
   errlin,errquad = checkDerivative(Dfunc,dD',centeredGrid)
   @test checkErrorDecay(errquad)


end
