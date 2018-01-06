using Images
 using Logging
 using DeformableRegistration: ImageProcessing, Transformation, Interpolation, Optimization
 using DeformableRegistration.Regularizer
 using DeformableRegistration.regOptions
 using Base.Test
 using DeformableRegistration.Optimization.augmentedLagrangian
 Logging.configure(level=Logging.DEBUG)
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
    @test size(dFunctionValue) == prod(gridSize)*2

end


@testset "constrainmts are met" begin
    x = ones(10)*2;
    index = zeros(size(x))
    index[1:3] = 1
    #x[index.==1]=0
    function zerosConstraintAtIndex(x::Array{Float64,1}, index)
        f = 0.5 * x .* x
        f[.!(index.==1)] = 0
        dF = spdiagm(x .* index,0)
        return [f, dF]
    end
    c(x) = zerosConstraintAtIndex(x, index)

    function J(x; doDerivative=true, doHessian=true)
        r = (x - Array(1:10))
        f = 0.5 * (r'*r)
        df = r
        d2F = speye(10)
        return [f, df, d2F]
    end

    options = regOptions()
    options.maxIterCG =100
    options.maxIterGaussNewton =100
    options.stopping["tolQ"]= 1e-20
    yopt = optimizeGaussNewtonAugmentedLagrangian(J, x, options, constraint=c, constraintObjectiveFunction=augmentedLagrangian)
    debug(@sprintf("... constraint norm:         %1.3e", norm(c(yopt)[1])))
    debug(@sprintf("... objective function norm: %1.3e", J(yopt)[1]))
    y_true = Array(1:10)
    y_true[index.==1] = 0
    debug(         "... error:                   ",yopt - y_true)
    debug(@sprintf("... mean error:              %1.3e",mean(abs.(yopt - y_true))))
    debug(@sprintf("... error norm:              %1.3e",norm(yopt - y_true)))
    @test norm(yopt - y_true)^2 < 0.15
end
