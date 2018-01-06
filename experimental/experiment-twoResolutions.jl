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
    reg = createDiffusiveOperatorCentered
    refImg, temImg, options = constructTestImages()
    options.stopping["tolQ"] = 1e-10

    options.maxIterCG = 10000
    plotRefImg, plotTemImg, plotOptions = constructTestImages()
    plotTemImg.data[13:12:end,13:12:end] = 1
    plotTemImg.data[14:12:end,13:12:end] = 1
    plotTemImg.data[13:12:end,14:12:end] = 1
    plotTemImg.data[14:12:end,14:12:end] = 1

    #"COARSE PRE-REG"
    close("all")
    options.levels = [5, 4]
    options.regularizerWeight = 10
    @time initialDisplacement = registerNonParametricFixedGridResolution(refImg, temImg, options,regularizerOperator=reg) # regularizerOperator=createDiffusiveOperatorCentered
    plot ? (figure();  visualizeResults(refImg, plotTemImg, displacement= initialDisplacement, showDeformationImage=showDeformationImage, suptitle="Coarse pre-registration", filename = "resultImg-coarse.png");) : 0
    ssdCoarse = ssdDistance(refImg, temImg, (initialDisplacement+getCellCenteredGrid(refImg)).data, doDerivative=true)[1]

    nPatches = (1,2)
    refImgs, temImgs, iDs = constructPatches(nPatches, refImg, temImg, initialDisplacement)





    #"FINE FULL REG"
    options.levels = [4, 3, 2]
    options.regularizerWeight = 10
    @time fineDisplacement = registerNonParametricConstraint(refImg, temImg, options,regularizerOperator=reg) # regularizerOperator=createDiffusiveOperatorCentered
    plot ? (figure();  visualizeResults(refImg, plotTemImg, displacement= fineDisplacement, showDeformationImage=showDeformationImage, suptitle="Fine full registration", filename = "resultImg-fine.png");) : 0
    ssdFine = ssdDistance(refImg, temImg, (fineDisplacement+getCellCenteredGrid(refImg)).data, doDerivative=true)[1]

    @printf("\n")

    @printf("\n
    SSD values
    ---------------------------------
    coarse registration:   %3.3e
    patch unconstraint:    %3.3e
    patch w/ constraints:  %3.3e
    fine full registration %3.3e\n", ssdCoarse, ssdFree, ssdConstraint, ssdFine)


    @printf("\n
    error norm (y-y_fine)
    ---------------------------------
    coarse registration:   %3.3e
    patch unconstraint:    %3.3e
    patch w/ constraints:  %3.3e
    fine full registration %3.3e\n",
    norm(initialDisplacement.data - fineDisplacement.data),
    norm(displacement_free.data - fineDisplacement.data),
    norm(displacement_constraint.data - fineDisplacement.data),
    norm(fineDisplacement.data - fineDisplacement.data))

    return initialDisplacement, displacement_free, displacement_constraint, fineDisplacement

end

initialDisplacement, displacement_constraint, displacement_free, fineDisplacement = testFixedGridRegistration(plot=true, showDeformationImage=true)

STOP
##
@printf("\n
error norm (y-y_fine)
---------------------------------
coarse registration:   %3.3e
patch w/ constraints:  %3.3e
patch free:            %3.3e
fine full registration %3.3e\n",
norm(initialDisplacement.data - fineDisplacement.data),
norm(displacement_constraint.data - fineDisplacement.data),
norm(displacement_free.data - fineDisplacement.data),
norm(fineDisplacement.data - fineDisplacement.data))

##
refImg, temImg, options = constructTestImages()


close("all"); f = figure();  visualizeResults(refImg, temImg, displacement= initialDisplacement, showDeformationImage=false, suptitle="test");

function reshapeGrid(g::scaledArray)
    fX = reshape(g.data[1:Int(length(g.data)/2)],g.dimensions)
    fY = reshape(g.data[Int(length(g.data)/2+1):end],g.dimensions)
    return fX, fY
end

cX, cY = reshapeGrid(initialDisplacement)
fX, fY = reshapeGrid(fineDisplacement)

figure(); PyPlot.imshow(cY); PyPlot.colorbar()

Ty = interpolateImage(temImg,getCellCenteredGrid(refImg),doDerivative=false)
 figure(); PyPlot.imshow(Ty)
 norm(Ty[:]-refImg.data[:])^2


ssdDistance(refImg, temImg, (fineDisplacement+getCellCenteredGrid(refImg)).data, doDerivative=true)


refImg.voxelsize

##
    #figure(); dx,dy = plotGrid(dispalcements_combined-initialDisplacement+getCellCenteredGrid(refImg))
    #m = initialDisplacement.dimensions
    #figure(); PyPlot.imshow(reshape(dispalcements_combined-initialDisplacement,(2*m[1],m[2])))
##

## ARTIFICIAL PRE-REG
initialDisplacement = getCellCenteredGrid(refImg)
 initialDisplacement.data[:] = 0.0
 initialDisplacement.data[1+Int(end/2):end] = -15.0

## CONSTRAINT TEST / PATCH-REG CONSTRAINT






## PATCH-REG
nPatches = (1,2)
refImgs, temImgs, iDs = constructPatches(nPatches, refImg, temImg, initialDisplacement)

displacements = Array{scaledArray}(2)
options.levels = [2]
 for i=1:2
 @time displacements[i] = registerImagesNonParametric(refImgs[i], temImgs[i], options,regularizerOperator=createCurvatureOperatorCentered, initialDisplacement=iDs[i]) # regularizerOperator=createDiffusiveOperatorCentered
 figure(); visualizeResults(refImgs[i], temImgs[i], displacement=displacements[i],showDeformationImage=true, suptitle="Patch registration , patch $(i)");
 end

 dispalcements_combined = combineDisplacementsNaive(displacements, nPatches)
 temImg.data[13:12:end,13:12:end] = 1

 figure(); visualizeResults(refImg, temImg, displacement=dispalcements_combined,showDeformationImage=true, suptitle="Patch registration , combined naively");
