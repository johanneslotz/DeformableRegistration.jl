using Images
 using Logging
 using DeformableRegistration: Visualization, Examples, ImageProcessing, Transformation, Interpolation, Distance
 using DeformableRegistration.Regularizer
 using DeformableRegistration.regOptions
 using Base.Test
 Logging.configure(level=Logging.INFO)
 include("./constraint.jl")
 include("./test-helpers.jl")


# function matchInitialDisplacementAtIndex(x::Array{Float64,1}, initialDisplacement::scaledArray, side; β = 0.01)
#     index = spzeros(initialDisplacement.dimensions...)
#     if (side == "left") || (side == 1)
#         index[:,end-1:end] = 1
#     elseif (side == "right") || (side == 2)
#         index[:,1:2] = 1
#     elseif (side == "boundary")
#         index[[1, end],:] = 1
#         index[:,[1, end]] = 1
#     else
#         throw("side needs to be left or right")
#     end
#     index = vcat(index[:], index[:])
#     f = β * 0.5 * index .* (x - initialDisplacement.data).^2
#     dF = β * spdiagm(index .* (x - initialDisplacement.data),0)
#     return [f, dF]
# end

function preRegistrationAsBoundaryConstraint(;plot=true,showDeformationImage=false)
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
    @time initialDisplacement = registerNonParametricConstraint(refImg, temImg, options,regularizerOperator=reg) # regularizerOperator=createDiffusiveOperatorCentered
    plot ? (figure();  visualizeResults(refImg, plotTemImg, displacement= initialDisplacement, showDeformationImage=showDeformationImage, suptitle="Coarse pre-registration", filename = "resultImg-coarse.png");) : 0
    ssdCoarse = ssdDistance(refImg, temImg, (initialDisplacement+getCellCenteredGrid(refImg)).data, doDerivative=true)[1]

    nPatches = (1,2)
    refImgs, temImgs, iDs = constructPatches(nPatches, refImg, temImg, initialDisplacement)

    #"CONSTRAINT PATCH-REG"
    displacements = Array{scaledArray}(2)
    options.levels = [2]
    options.regularizerWeight = 0.01
    for i=1:2
       c(x, initial) = matchInitialDisplacementAtIndex(x, initial, i)
       @time displacements[i] = registerNonParametricConstraint(refImgs[i], temImgs[i], options, regularizerOperator=reg, initialDisplacement=iDs[i], constraint = c) # regularizerOperator=createDiffusiveOperatorCentered
       plot ? (figure(); visualizeResults(refImgs[i], temImgs[i], displacement=displacements[i],showDeformationImage=showDeformationImage, suptitle="Constraint Patch registration , patch $(i)", filename = "resultImg-constraint-patch-$(i).png");) : 0
    end

    displacement_constraint = combineDisplacementsNaive(displacements, nPatches)
    plot ? (figure(); visualizeResults(refImg, plotTemImg, displacement=displacement_constraint,showDeformationImage=showDeformationImage, suptitle="Constraint Patch registration , combined by padding", filename = "resultImg-constraint-padding.png");) : 0
    ssdConstraint = ssdDistance(refImg, plotTemImg, (displacement_constraint+getCellCenteredGrid(refImg)).data, doDerivative=true)[1]
    #"FREE PATCH-REG"
    displacements = Array{scaledArray}(2)
    options.levels = [2]
    options.regularizerWeight = 10
    for i=1:2
       c(x, initial) = matchInitialDisplacementAtIndex(x, initial, i)
       @time displacements[i] = registerNonParametricConstraint(refImgs[i], temImgs[i], options, regularizerOperator=reg, initialDisplacement=iDs[i]) # regularizerOperator=createDiffusiveOperatorCentered
       plot ? (figure(); visualizeResults(refImgs[i], temImgs[i], displacement=displacements[i],showDeformationImage=showDeformationImage, suptitle="Unconstraint Patch registration , patch $(i)", filename = "resultImg-free-patch-$(i).png");) : 0
    end

    displacement_free = combineDisplacementsNaive(displacements, nPatches)
    plot ? (figure(); visualizeResults(refImg, plotTemImg, displacement=displacement_free,showDeformationImage=showDeformationImage, suptitle="Unconstraint Patch registration , combined by padding", filename = "resultImg-free-padding.png");) : 0
    ssdFree = ssdDistance(refImg, temImg, (displacement_free+getCellCenteredGrid(refImg)).data, doDerivative=true)[1]

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

initialDisplacement, displacement_constraint, displacement_free, fineDisplacement = preRegistrationAsBoundaryConstraint(plot=true, showDeformationImage=true)

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


function compareRegistrationWithOneLevelFiner()
## COMPARE REGISTRATION WITH ONE FINER LEVEL
options.levels = [3]
 @time displacement1 = registerImagesNonParametric(refImg, temImg, options,regularizerOperator=createCurvatureOperatorCentered) # regularizerOperator=createDiffusiveOperatorCentered
 ssd = ssdDistance(refImg, temImg, (displacement1+getCellCenteredGrid(displacement)).data, doDerivative=false)[1]
 println(ssd)

options.levels = [3,2]
 @time displacement2 = registerImagesNonParametric(refImg, temImg, options,regularizerOperator=createCurvatureOperatorCentered) # regularizerOperator=createDiffusiveOperatorCentered
 ssd = ssdDistance(refImg, temImg, (displacement2+getCellCenteredGrid(displacement)).data, doDerivative=false)[1]
 println(ssd)

figure()
clf(); visualizeResults(refImg,temImg;
                          displacement=displacement2)

diffImg = deepcopy(refImg)
 diffImg.data[:] = refImg.data - Array(temImg.data)
 figure()
 subplot(121)
 imshow(restrictResolutionToLevel(diffImg,3).data)
 subplot(122)
 imshow(restrictResolutionToLevel(diffImg,2).data)
end
##
# Profile.print()
# open("/tmp/prof.txt", "w") do s
#     Profile.print(IOContext(s, :displaysize => (24, 500)))
# end

refImg, temImg, options = constructTestImages()
 refImg.data[:] = 0.0
 refImg.data[30:90,:30:90] = 1.0
 temImg.data[:] = 0.0
 temImg.data[30:90,:40:100] = 1.0
 referenceGrid = getCellCenteredGrid(refImg)
 ## ARTIFICIAL PRE-REG
 testDisplacement = getCellCenteredGrid(restrictResolutionToLevel(refImg,2))
 testDisplacement.data[:] = 0.0
 testDisplacementHD = interpolateDeformationField(testDisplacement, referenceGrid, interpolationScheme=BSpline(Linear()))

 ssdTest1 = ssdDistance(refImg, temImg, (testDisplacementHD+getCellCenteredGrid(refImg)).data, doDerivative=false)[1]
 testDisplacement.data[1+Int(end/2):end] = 10.0
 testDisplacementHD = interpolateDeformationField(testDisplacement, referenceGrid, interpolationScheme=BSpline(Linear()))
 ssdTest2 = ssdDistance(refImg, temImg, (testDisplacementHD+getCellCenteredGrid(refImg)).data, doDerivative=false)[1]

 println(ssdTest1)
 println(ssdTest2)

 close("all");figure();visualizeResults(refImg, temImg, displacement=testDisplacementHD,showDeformationImage=false, suptitle="TEST");
