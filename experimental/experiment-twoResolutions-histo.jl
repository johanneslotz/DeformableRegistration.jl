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

 include("../histo/HistokatControllerImageLoader.jl")
 using HistokatControllerImageLoader: getTile, getImageObject
 using PyPlot
##
 filenameR = "/Users/jo/data/example-data-LL1_1_CD146-2014.tif"
 filenameT = "/Users/jo/data/example-data-LL1_4_KL1-2014.tif"

 imageR = getImageObject(filenameR)
 imageT = getImageObject(filenameT)

function removeBlack!(arr::Array)
     for i = 1:length(arr)
         if arr[i]<5
             arr[i]=255
         end
     end
 end

function normalizeImage!(arr::Array)
    removeBlack!(arr)
    arr[:] = arr[:] .* (arr[:] .< 200)
end

refImgCoarse = getTile(imageR, 9, 0,0,0)
temImgCoarse = getTile(imageT, 9, 0,0,0)

normalizeImage!(refImgCoarse.data.data)
normalizeImage!(temImgCoarse.data.data)

targetGridSize = (16,32)
targetVoxelSize = size(refImgCoarse.data) .* refImgCoarse.voxelsize ./ targetGridSize
targetGrid = getCellCenteredGrid(targetVoxelSize, refImgCoarse.shift, targetGridSize)

options = getSomeStandardOptions()
reg=createDiffusiveOperatorCentered
interpolationScheme=BSpline(Linear())

plot = true
# "PRE-REG"
plot ? close("all") : 0
options.regularizerWeight = 1000
        options.levels = [5,4,3,2, 1]
        preRegResult = registerNonParametricOnTargetGrid(refImgCoarse, temImgCoarse, targetGrid,
            options,regularizerOperator=reg,
            interpolationScheme=interpolationScheme)
            # regularizerOperator=createDiffusiveOperatorCentered
            tempDisplacement = interpolateDeformationField(
            preRegResult, getCellCenteredGrid(refImgCoarse))
        plot ? (
            figure();
            visualizeResults(refImgCoarse, temImgCoarse,
                displacement= tempDisplacement,
                showDeformationImage=false,
                suptitle="pre registration", filename = "");) : 0

# "COARSE + PATCH COMBINED"
options.regularizerWeight = 3
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
        #plot ? (figure();  visualizeResults(refImgFine, plotTemImg, displacement= tempDisplacement, showDeformationImage=showDeformationImage, suptitle="fixed grid resolution at level", filename = "");) : 0
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


# "ALL MULTILEVEL"
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
