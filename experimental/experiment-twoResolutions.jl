using Images
 using Logging
 using DeformableRegistration: Visualization, Examples, ImageProcessing, Transformation, Interpolation, Distance
 using DeformableRegistration.Regularizer
 using DeformableRegistration.regOptions
 using Base.Test
 Logging.configure(level=Logging.WARNING)
 include("./constraint.jl")
 include("./test-helpers.jl")
 include("./fixedGridRegistration.jl")

plot=true
showDeformationImage = false
function testFixedGridRegistration(;plot=true,showDeformationImage=false)

refImg, temImg, options = constructTestImages()
    options.stopping["tolQ"] = 1e-10
    options.stopping["tolJ"] = 1e-10
    options.stopping["tolY"] = 1e-10
    options.stopping["tolG"] = 1e-10

    reg=createDiffusiveOperatorCentered

    #interpolationScheme=BSpline(Linear())
    interpolationScheme=InterpLinearFast

    m = getCellCenteredGrid(refImg).dimensions
    B = createDiffusiveOperatorCentered([1.0,1.0],[m[1],m[2]])
    B=B'*B

    options.interpolateToReferenceImage = true

    options.maxIterCG = 10000
    plotRefImg, plotTemImg, plotOptions = constructTestImages()
    plotTemImg.data[13:12:end,13:12:end] = 1
    plotTemImg.data[14:12:end,13:12:end] = 1
    plotTemImg.data[13:12:end,14:12:end] = 1
    plotTemImg.data[14:12:end,14:12:end] = 1

    #"FIXED GRID RESOLUTION"
    close("all")
options.levels = [4, 3, 2, 1, 0]
    options.regularizerWeight = 100
    @time singleGridDisplacement = registerNonParametricFixedGridResolution(refImg, temImg, options,
                                    regularizerOperator=reg, interpolationScheme=interpolationScheme, gridLevel=4) # regularizerOperator=createDiffusiveOperatorCentered
    plot ? (figure();  visualizeResults(refImg, plotTemImg, displacement= singleGridDisplacement, showDeformationImage=showDeformationImage, suptitle="fixed grid resolution at level", filename = "");) : 0
    ssdCoarse = ssdDistance(refImg, temImg, (singleGridDisplacement+getCellCenteredGrid(refImg)).data, doDerivative=true)[1]
    regularizerCoarse = singleGridDisplacement.data[:]'*B*singleGridDisplacement.data[:]


    #"ALL MULTILEVEL"
options.levels = [4, 3, 2, 1]
    options.regularizerWeight = 1
    @time fineDisplacement = registerNonParametricConstraint(refImg, temImg, options,regularizerOperator=reg, interpolationScheme=interpolationScheme) # regularizerOperator=createDiffusiveOperatorCentered
    plot ? (figure();  visualizeResults(refImg, plotTemImg, displacement= fineDisplacement, showDeformationImage=showDeformationImage, suptitle="Fine full registration", filename = "");) : 0
    ssdFine = ssdDistance(refImg, temImg, (fineDisplacement+getCellCenteredGrid(refImg)).data, doDerivative=true)[1]
    regularizerFine = fineDisplacement.data[:]'*B*fineDisplacement.data[:]


@printf("\n")

    @printf("\n
    SSD values
    ---------------------------------
    single grid:           %3.3e
    fine full registration %3.3e\n", ssdCoarse, ssdFine)


    @printf("\n
    error norm (y-y_fine)
    ---------------------------------
    single Grid:           %3.3e
    fine full registration %3.3e\n",
    norm(singleGridDisplacement.data - fineDisplacement.data),
    norm(fineDisplacement.data - fineDisplacement.data))


    @printf("\n
    REGULARIZER values
    ---------------------------------
    single Grid           %3.3e
    full registratio      %3.3e
    gradient descent      %3.3e\n",  options.regularizerWeight * regularizerCoarse,  options.regularizerWeight * regularizerFine,0)

    @printf("\n
    J values
    ---------------------------------
    single Grid           %3.3e
    full registratio      %3.3e
    gradient descent      %3.3e\n", ssdCoarse + options.regularizerWeight * regularizerCoarse, ssdFine + options.regularizerWeight * regularizerFine,0)


    return singleGridDisplacement, fineDisplacement

end

initialDisplacement, fineDisplacement = testFixedGridRegistration(plot=true, showDeformationImage=true)
