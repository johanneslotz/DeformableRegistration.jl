using Images
 using Logging
 using DeformableRegistration: Visualization, Examples, ImageProcessing, Transformation, Interpolation, Distance
 using DeformableRegistration.Regularizer
 using DeformableRegistration.regOptions
 using Base.Test
 Logging.configure(level=Logging.INFO)
 include("./constraint.jl")
 #include("./aspin.jl")
 include("./test-helpers.jl")
 include("./twoResolutionRegistration.jl")
 plot=true
 showDeformationImage = false

#
#function testTwoResolutionsRegistration(;plot=true,showDeformationImage=false)

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



# "ALL MULTILEVEL"
options.regularizerWeight = 1
    options.levels = [4, 3, 2, 1]
    @time fineDisplacement = registerNonParametricConstraint(
        refImgFine, temImgFine, options,regularizerOperator=reg,
        interpolationScheme=interpolationScheme)
    plot ? (figure();  visualizeResults(refImgFine, plotTemImg, displacement= fineDisplacement, showDeformationImage=showDeformationImage, suptitle="Fine full registration", filename = "artificial-example-DDS-fine-level-$(options.levels[end]).png");) : 0
    ssdFine = ssdDistance(refImgFine, temImgFine, (fineDisplacement+getCellCenteredGrid(refImgFine)).data, doDerivative=true)[1]
    regularizerFine = fineDisplacement.data[:]'*B*fineDisplacement.data[:]


# "ASPIN"
options.regularizerWeight = 1
        options.levels = [0]
        options.maxIterGaussNewton = 5
        imgMask = 0*refImgFine.data
        imgMask[1:120,1:120] = 1
        S = zeros(2,length(imgMask[:]))
        S[1,:] = imgMask[:]

        imgMask = 0*refImgFine.data
        imgMask[1:120,121:240] = 1
        S[2,:] = imgMask[:]

        @time aspinDisplacement = registerNonParametricConstraint(
            refImgFine, temImgFine, options,regularizerOperator=reg,initialDisplacement = fineDisplacement,
            interpolationScheme=interpolationScheme, subspaces = S)
        plot ? (figure();  visualizeResults(refImgFine, plotTemImg, displacement= aspinDisplacement, showDeformationImage=showDeformationImage, suptitle="Fine full registration", filename = "artificial-example-DDS-fine-level-$(options.levels[end]).png");) : 0
        ssdFine = ssdDistance(refImgFine, temImgFine, (aspinDisplacement+getCellCenteredGrid(refImgFine)).data, doDerivative=true)[1]
        regularizerFine = aspinDisplacement.data[:]'*B*aspinDisplacement.data[:]


#
# @printf("\n")
#
#     @printf("\n
#     SSD values
#     ---------------------------------
#     pre-reg:               %3.3e
#     combined patches       %3.3e
#     fine full registration %3.3e\n", ssdCoarse, ssdCombined,  ssdFine)
#
#
#     @printf("\n
#     error norm (y-y_fine)
#     ---------------------------------
#     pre-reg:               %3.3e
#     combined patches       %3.3e
#     fine full registration %3.3e\n",
#     norm(tempDisplacement.data - fineDisplacement.data),
#     norm(combinedDisplacement.data - fineDisplacement.data),
#     norm(fineDisplacement.data - fineDisplacement.data))
#
#
#     @printf("\n
#     REGULARIZER values
#     ---------------------------------
#     pre-reg:               %3.3e
#     combined patches       %3.3e
#     fine full registration %3.3e\n",
#         options.regularizerWeight * regularizerCoarse,
#         options.regularizerWeight * regularizerCombined,
#         options.regularizerWeight * regularizerFine)
#
#     @printf("\n
#     J values
#     ---------------------------------
#     pre-reg:               %3.3e
#     combined patches       %3.3e
#     fine full registration %3.3e\n",
#         ssdCoarse + options.regularizerWeight * regularizerCoarse,
#         ssdCombined + options.regularizerWeight * regularizerCombined,
#         ssdFine + options.regularizerWeight * regularizerFine)
# end
#
# testTwoResolutionsRegistration(plot=true, showDeformationImage=false)
