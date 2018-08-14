using DeformableRegistration: ImageProcessing, Distance, Transformation, Interpolation, Examples, Distance, regOptions
 using Base.Test
#  using Logging

@testset "Parametric Registration" begin
@testset "SSD/NGF shifted and scaled"  for D=[ssdDistance, ngfDistance]
##
    testimage = dirname(Base.source_path()) * "/testdata/luebeck.jpg"
    refImg = loadImage(testimage)
    # define a cell centered grid, transform it and create a template image
    centeredGrid = getCellCenteredGrid(refImg)
    affineParametersInitial = [1.2,0.0,-50,0.0,1.2,-0]
    transformedGrid = transformGridAffine(centeredGrid,affineParametersInitial)
    temImg = interpolateImage(refImg,transformedGrid,interpolationScheme=InterpLinearFast)[1]
    temImg = createImage(temImg)
    options = regOptions()
    options.levels = [4,3]
    options.parametricOnly = true;
    options.matrixFree = true;
    options.tolCG = 1e-3
    options.useEdgeParameterInNumerator = true
    edgeParameterR = 0.1
    edgeParameterT = 0.1

    @info D

    affineParameters = registerImagesParametric(refImg,temImg, options, measureDistance=D)
    if D==DeformableRegistration.Distance.ssdDistance
        @test affineParameters[1] ≈ 1/affineParametersInitial[1] atol=0.1
        @test affineParameters[2] ≈ 0  atol=1e-2
        @test affineParameters[3] ≈ -affineParametersInitial[3]/affineParametersInitial[1]  atol=5
        @test affineParameters[4] ≈ 0  atol=1e-2
        @test affineParameters[5] ≈ 1/affineParametersInitial[5]  atol=0.1
    else
        @test affineParameters[1] ≈ 1/affineParametersInitial[1] atol=0.1
        @test affineParameters[2] ≈ 0  atol=5e-2
        @test affineParameters[3] ≈ -affineParametersInitial[3]/affineParametersInitial[1]  atol=5
        @test affineParameters[4] ≈ 0  atol=1e-1
        @test affineParameters[5] ≈ 1/affineParametersInitial[5]  atol=0.5
    end
end
end
