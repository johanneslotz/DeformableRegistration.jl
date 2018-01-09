using Images
 using Logging
 using DeformableRegistration: Visualization, Examples, ImageProcessing, Transformation, Interpolation, Distance
 using DeformableRegistration.Regularizer
 using DeformableRegistration.regOptions
 using Base.Test
 Logging.configure(level=Logging.INFO)
 include("./constraint.jl")
 include("./test-helpers.jl")
 include("./fixedGridRegistration.jl")

function testFixedGridRegistration(;plot=true,showDeformationImage=false)
    refImg, temImg, options = constructTestImages()
    options.stopping["tolQ"] = 1e-10

    reg=createDiffusiveOperatorCentered

    interpolationScheme=BSpline(Linear())
    #interpolationScheme=InterpLinearFast

    options.maxIterCG = 10000
    plotRefImg, plotTemImg, plotOptions = constructTestImages()
    plotTemImg.data[13:12:end,13:12:end] = 1
    plotTemImg.data[14:12:end,13:12:end] = 1
    plotTemImg.data[13:12:end,14:12:end] = 1
    plotTemImg.data[14:12:end,14:12:end] = 1

    #"FIXED GRID RESOLUTION"
    close("all")
    options.levels = [4, 3, 2, 1, 0]
    options.regularizerWeight = 2
    @time initialDisplacement = registerNonParametricFixedGridResolution(refImg, temImg, options,regularizerOperator=reg, interpolationScheme=interpolationScheme) # regularizerOperator=createDiffusiveOperatorCentered
    plot ? (figure();  visualizeResults(refImg, plotTemImg, displacement= initialDisplacement, showDeformationImage=showDeformationImage, suptitle="fixed grid resolution at level 4", filename = "");) : 0
    ssdCoarse = ssdDistance(refImg, temImg, (initialDisplacement+getCellCenteredGrid(refImg)).data, doDerivative=true)[1]


    #"ALL MULTILEVEL"
    options.levels = [4, 3, 2, 1, 0]
    options.regularizerWeight = 1
    @time fineDisplacement = registerNonParametricConstraint(refImg, temImg, options,regularizerOperator=reg, interpolationScheme=interpolationScheme) # regularizerOperator=createDiffusiveOperatorCentered
    plot ? (figure();  visualizeResults(refImg, plotTemImg, displacement= fineDisplacement, showDeformationImage=showDeformationImage, suptitle="Fine full registration", filename = "");) : 0
    ssdFine = ssdDistance(refImg, temImg, (fineDisplacement+getCellCenteredGrid(refImg)).data, doDerivative=true)[1]

    @printf("\n")

    @printf("\n
    SSD values
    ---------------------------------
    coarse registration:   %3.3e
    fine full registration %3.3e\n", ssdCoarse, ssdFine)


    @printf("\n
    error norm (y-y_fine)
    ---------------------------------
    coarse registration:   %3.3e
    fine full registration %3.3e\n",
    norm(initialDisplacement.data - fineDisplacement.data),
    norm(fineDisplacement.data - fineDisplacement.data))

    return initialDisplacement, fineDisplacement

end

initialDisplacement, fineDisplacement = testFixedGridRegistration(plot=true, showDeformationImage=true)
