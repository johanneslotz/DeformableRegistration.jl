using Images
 using Logging
 using DeformableRegistration: ImageProcessing, Transformation, Interpolation, Optimization
 using DeformableRegistration.Regularizer
 using DeformableRegistration.Distance
 using DeformableRegistration.regOptions
 using Base.Test
 using DeformableRegistration.Optimization.augmentedLagrangian
 using Interpolations
 Logging.configure(level=Logging.INFO)
 include("./twoResolutionRegistration.jl")


##
@testset "fixedGrid" begin







@testset "speed ssd arbitrary grid" begin
##
    n=128
    offset = Int(round(n/10))
    Rdata = zeros(n,n)
    left = Int(n/4)
    right = Int(n/2)
    Rdata[left:right,left:right] = 1.0

    Tdata = zeros(n,n)
    Tdata[1:end-offset,:] = Rdata[offset+1:end,:]
    Rdata[n-2*offset:n-offset, n-2*offset:n-offset] = 1
    Tdata[n-2*offset:n-offset, n-2*offset:n-offset] = 1

    R = createImage(Rdata)
    T = createImage(Tdata)

    targetGridSize = (128,128)
    targetVoxelSize = size(R.data) .* R.voxelsize ./ targetGridSize
    targetGrid = getCellCenteredGrid(targetVoxelSize, R.shift, targetGridSize)

    options = regOptions()
    options.matrixFree = true;
    options.interpolateToReferenceImage = true
    options.regularizerWeight = 1
    options.stopping["tolQ"] = 1e-4
    options.maxIterCG = 4000
    options.levels=[1,0]

    @time displacement = registerNonParametricOnTargetGrid(R, T, targetGrid,
    options,regularizerOperator=createDiffusiveOperatorCentered,
    interpolationScheme=BSpline(Linear()), gradientDescentOnly=false)
##
    #tempDisplacement = interpolateDeformationField(
#        displacement, getCellCenteredGrid(R))
    # using PyPlot
    # using DeformableRegistration: Visualization
    # figure();
    # visualizeResults(R, T,
    #         displacement= displacement,
    #         showDeformationImage=true,
    #         suptitle="registration", filename = "")

end

##

end
