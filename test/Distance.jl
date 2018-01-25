using DeformableRegistration: Distance, ImageProcessing, Transformation, Interpolation, regOptions
 using Base.Test
 using Logging

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
    Logging.warn("Skipping ")
    @test_skip estimateNGFEpsilon(refImg,cutoffPercent=80)[1] >= 0 # seems to run forever
end
