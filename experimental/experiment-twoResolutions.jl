using Images
 using Logging
 using DeformableRegistration: Visualization, Examples, ImageProcessing, Transformation, Interpolation, Distance
 using DeformableRegistration.Regularizer
 using DeformableRegistration.regOptions
 using Base.Test
 Logging.configure(level=Logging.INFO)
 include("./constraint.jl")
 include("./test-helpers.jl")
 include("./twoResolutionRegistration.jl")

plot=true
showDeformationImage = false
function testTwoResolutionsRegistration(;plot=true,showDeformationImage=false)

##
refPatches, temPatches, refImgCoarse, temImgCoarse, options, patchLevel, imageLevel =
        constructTestImagesAndPatches(patchLevel=1, imageLevel=3)
    targetGridSize = (16,32)
    targetVoxelSize = size(refImgCoarse.data) .* refImgCoarse.voxelsize ./ targetGridSize
    targetGrid = getCellCenteredGrid(targetVoxelSize, refImgCoarse.shift, targetGridSize)
    refImgFine, temImgFine, _o = constructTestImages()
    options.stopping["tolQ"] = 1e-12
    options.stopping["tolJ"] = 1e-12
    options.stopping["tolY"] = 1e-12
    options.stopping["tolG"] = 1e-12

    reg=createDiffusiveOperatorCentered

    interpolationScheme=BSpline(Linear())
    #interpolationScheme=InterpLinearFast

    options.interpolateToReferenceImage = true

    options.maxIterCG = 10000
    plotRefImg, plotTemImg, plotOptions = constructTestImages()
    plotTemImg.data[13:12:end,13:12:end] = 1
    plotTemImg.data[14:12:end,13:12:end] = 1
    plotTemImg.data[13:12:end,14:12:end] = 1
    plotTemImg.data[14:12:end,14:12:end] = 1


    m = getCellCenteredGrid(refImgFine).dimensions
    B = createDiffusiveOperatorCentered([1.0,1.0],[m[1],m[2]])
    B=B'*B


## "PRE-REG"
plot ? close("all") : 0
options.regularizerWeight = 1
        options.levels = [1, 0]
        preRegResult = registerNonParametricOnTargetGrid(refImgCoarse, temImgCoarse, targetGrid,
            options,regularizerOperator=reg,
            interpolationScheme=interpolationScheme) # regularizerOperator=createDiffusiveOperatorCentered
        tempDisplacement = interpolateDeformationField(
            preRegResult, getCellCenteredGrid(refImgFine))
        plot ? (
            figure();
            visualizeResults(refImgFine, plotTemImg,
                displacement= tempDisplacement,
                showDeformationImage=showDeformationImage,
                suptitle="pre registration", filename = "");) : 0
        ssdCoarse = ssdDistance(refImgFine, temImgFine,
            (tempDisplacement+getCellCenteredGrid(refImgFine)).data, doDerivative=true)[1]
        regularizerCoarse = tempDisplacement.data[:]'*B*tempDisplacement.data[:]

## "COARSE + PATCH COMBINED"
options.regularizerWeight = 0.1
displacements = Array{scaledArray}(2)
    @time for i = 1:2
        displacements[i] = registerInTwoResolutions(refPatches[i], temPatches[i],
                    refImgCoarse, temImgCoarse, targetGrid, options,
                    initialDisplacement = preRegResult,
                    regularizerOperator=reg, interpolationScheme=interpolationScheme,
                    patchWeight = 1,
                    imageWeight = 1,
                    imageLevel=imageLevel,
                    patchLevel=patchLevel,
                    gradientDescentOnly=false) # regularizerOperator=createDiffusiveOperatorCentered
        tempDisplacement = interpolateDeformationField(displacements[i], getCellCenteredGrid(refImgFine))
        plot ? (figure();  visualizeResults(refImgFine, plotTemImg, displacement= tempDisplacement, showDeformationImage=showDeformationImage, suptitle="fixed grid resolution at level", filename = "");) : 0
    end
    #singleGridDisplacement = interpolateDeformationField(singleGridDisplacement, getCellCenteredGrid(refImgFine))
    #plot ? (close("all"); figure();  visualizeResults(refImgFine, plotTemImg, displacement= singleGridDisplacement, showDeformationImage=showDeformationImage, suptitle="fixed grid resolution at level", filename = "");) : 0
    #ssdCoarse = ssdDistance(refImgFine, temImgFine, (singleGridDisplacement+getCellCenteredGrid(refImgFine)).data, doDerivative=true)[1]
    #regularizerCoarse = singleGridDisplacement.data[:]'*B*singleGridDisplacement.data[:]

    displacements[1] = stripDisplacement(displacements[1],(1:16,1:16))
    displacements[2] = stripDisplacement(displacements[2],(1:16,17:32))
    combinedDisplacement = combineDisplacementsNaive(displacements,(1,2))



    combinedDisplacement = interpolateDeformationField(combinedDisplacement, getCellCenteredGrid(refImgFine))
    plot ? (figure();  visualizeResults(refImgFine, plotTemImg, displacement= combinedDisplacement, showDeformationImage=showDeformationImage, suptitle="combined patches", filename = "");) : 0
    ssdCombined = ssdDistance(refImgFine, temImgFine, (combinedDisplacement+getCellCenteredGrid(refImgFine)).data, doDerivative=true)[1]
    regularizerCombined = combinedDisplacement.data[:]'*B*combinedDisplacement.data[:]


## "ALL MULTILEVEL"
options.regularizerWeight = 1
    options.levels = [4, 3, 2, 1]
    @time fineDisplacement = registerNonParametricConstraint(
        refImgFine, temImgFine, options,regularizerOperator=reg,
        interpolationScheme=interpolationScheme)
    plot ? (figure();  visualizeResults(refImgFine, plotTemImg, displacement= fineDisplacement, showDeformationImage=showDeformationImage, suptitle="Fine full registration", filename = "");) : 0
    ssdFine = ssdDistance(refImgFine, temImgFine, (fineDisplacement+getCellCenteredGrid(refImgFine)).data, doDerivative=true)[1]
    regularizerFine = fineDisplacement.data[:]'*B*fineDisplacement.data[:]


@printf("\n")

    @printf("\n
    SSD values
    ---------------------------------
    pre-reg:               %3.3e
    combined patches       %3.3e
    fine full registration %3.3e\n", ssdCoarse, ssdCombined,  ssdFine)


    @printf("\n
    error norm (y-y_fine)
    ---------------------------------
    pre-reg:               %3.3e
    combined patches       %3.3e
    fine full registration %3.3e\n",
    norm(tempDisplacement.data - fineDisplacement.data),
    norm(combinedDisplacement.data - fineDisplacement.data),
    norm(fineDisplacement.data - fineDisplacement.data))


    @printf("\n
    REGULARIZER values
    ---------------------------------
    pre-reg:               %3.3e
    combined patches       %3.3e
    fine full registration %3.3e\n",
        options.regularizerWeight * regularizerCoarse,
        options.regularizerWeight * regularizerCombined,
        options.regularizerWeight * regularizerFine)

    @printf("\n
    J values
    ---------------------------------
    pre-reg:               %3.3e
    combined patches       %3.3e
    fine full registration %3.3e\n",
        ssdCoarse + options.regularizerWeight * regularizerCoarse,
        ssdCombined + options.regularizerWeight * regularizerCombined,
        ssdFine + options.regularizerWeight * regularizerFine)
end

testTwoResolutionsRegistration(plot=true, showDeformationImage=false)
