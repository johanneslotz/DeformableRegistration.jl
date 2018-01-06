using Images
 using Logging
 using DeformableRegistration: Visualization, Examples, ImageProcessing, Transformation, Interpolation, Distance
 using DeformableRegistration.Regularizer
 using DeformableRegistration.regOptions
 using Base.Test
 Logging.configure(level=Logging.WARNING)
 #Logging.configure(level=Logging.INFO)
 include("./constraint.jl")
 include("./test-helpers.jl")
 include("./ssdDistance-noHessian.jl")

 function createDiffusiveOperatorDirichlet(h::Array{Float64,1},m::Array{Int64,1})

   d(k) = spdiagm((-ones(m[k],1),ones(m[k],1)),[-1 0],m[k]+1,m[k])/(h[k])
   dx = d(1)
   dy = d(2)
   D1 = kron(speye(m[2]),dx)
   D2 = kron(dy,speye(m[1]))
   p1,p2 = size(D1)
   p3,p4 = size(D2)
   B = [ D1 spzeros(p1,p2);
         D2 spzeros(p3,p4);
         spzeros(p1,p2) D1;
         spzeros(p3,p4) D2
       ]
 end

#B = createDiffusiveOperatorDirichlet([1.0,1.0],[6,6]); imshow(B'*B)
function registerFineGlobal(;plot=true,showDeformationImage=false, reg = createDiffusiveOperatorDirichlet)
    refImg, temImg, options = constructTestImages()
    options.stopping["tolQ"] = 1e-6
    options.stopping["tolY"] = 1e-6
    options.stopping["tolG"] = 1e-6
    options.stopping["tolJ"] = 1e-6

    m = getCellCenteredGrid(refImg).dimensions
    B = createDiffusiveOperatorDirichlet([1.0,1.0],[m[1],m[2]])
    B=B'*B

    options.maxIterCG = 10000
    options.maxIterGaussNewton = 100000;

    plotRefImg, plotTemImg, plotOptions = constructTestImages()
    plotTemImg.data[13:12:end,13:12:end] = 1
    plotTemImg.data[14:12:end,13:12:end] = 1
    plotTemImg.data[13:12:end,14:12:end] = 1
    plotTemImg.data[14:12:end,14:12:end] = 1

    options.levels = [4, 3, 2, 1]
    options.regularizerWeight = 1

    print("""

    ===============   FULL SSD   =================

    """)
    #"FINE FULL REG"
    @time fineDisplacementFull = registerNonParametricConstraint(refImg, temImg, options,regularizerOperator=reg)
    plot ? (figure();  visualizeResults(refImg, plotTemImg,  numberOfGridLines = 30, displacement= fineDisplacementFull, showDeformationImage=showDeformationImage, suptitle="full SSD", filename = "resultImg-fullssd.png");) : 0
    ssdFineFull = ssdDistance(refImg, temImg, (fineDisplacementFull+getCellCenteredGrid(refImg)).data, doDerivative=true)[1]
    regularizerFull = fineDisplacementFull.data[:]'*B*fineDisplacementFull.data[:]

    print("""

    =============== REDUCED SSD =================

    """)
    @time fineDisplacementReduced = registerNonParametricConstraint(refImg, temImg, options,regularizerOperator=reg, measureDistance=ssdDistanceNoHessian)
    plot ? (figure();  visualizeResults(refImg, plotTemImg, numberOfGridLines = 30, displacement= fineDisplacementReduced, showDeformationImage=showDeformationImage, suptitle="SSD wo/ Hessian", filename = "resultImg-ssdWOHessian.png");) : 0
    ssdFineReduced = ssdDistance(refImg, temImg, (fineDisplacementReduced+getCellCenteredGrid(refImg)).data, doDerivative=true)[1]
    regularizerReduced = fineDisplacementReduced.data[:]'*B*fineDisplacementReduced.data[:]

    plot ? (figure();  visualizeResults(refImg, plotTemImg, displacement= fineDisplacementReduced-fineDisplacementFull, showDeformationImage=true, suptitle="difference reduced operator", filename = "resultImg-diff-reduced-operator.png");) : 0


    print("""

    =============== GRADIENT DESCENT SSD =================

    """)
    ssdGradientDescent = 0
    # @time fineDisplacementGD = registerNonParametricConstraint(refImg, temImg, options,regularizerOperator=reg, measureDistance=ssdDistanceNoHessian, gradientDescentOnly=true)
    # plot ? (figure();  visualizeResults(refImg, plotTemImg, displacement= fineDisplacementGD, showDeformationImage=showDeformationImage, suptitle="SSD gradient descent", filename = "resultImg-ssdWOHessian.png");) : 0
    # ssdGradientDescent = ssdDistance(refImg, temImg, (fineDisplacementGD+getCellCenteredGrid(refImg)).data, doDerivative=true)[1]
    #
    # plot ? (figure();  visualizeResults(refImg, plotTemImg, displacement= fineDisplacementGD-fineDisplacementFull, showDeformationImage=true, suptitle="difference gradient descent", filename = "resultImg-diff-reduced-operator.png");) : 0




    # print("""
    #
    # =============== REDUCED SSD, imhom. Dirichlet =================
    #
    # """)
    #
    # #"FINE FULL REG reduced OP"
    # options.levels = [4, 3, 2]
    # options.regularizerWeight = 10
    #
    # initialDisplacement = getCellCenteredGrid(refImg)
    # initialDisplacement.data[:] = 0.0
    # initialDisplacement.data[Int(end/2)+1:end] = 10.0
    #
    # c(x, initial) = matchInitialDisplacementAtIndex(x, initial, "boundary",β = 0.1)
    # @time fineDisplacementReducedInHom = registerNonParametricConstraint(refImg, temImg, options, initialDisplacement=initialDisplacement, regularizerOperator=reg, measureDistance=ssdDistanceNoHessian, constraint = c)
    # plot ? (figure();  visualizeResults(refImg, plotTemImg, displacement= fineDisplacementReducedInHom, showDeformationImage=true, suptitle="SSD wo/ Hessian inhom dirichlet", filename = "resultImg-ssdWOHessianInHom.png");) : 0
    # ssdFineReducedInHom = ssdDistance(refImg, temImg, (fineDisplacementReducedInHom+getCellCenteredGrid(refImg)).data, doDerivative=true)[1]
    #
    #



    @printf("\n
    SSD values
    ---------------------------------
    full registration     %3.3e
    reduced operator      %3.3e
    gradient descent      %3.3e\n", ssdFineFull, ssdFineReduced,ssdGradientDescent)

    @printf("\n
    REGULARIZER values
    ---------------------------------
    full registration     %3.3e
    reduced operator      %3.3e
    gradient descent      %3.3e\n",  options.regularizerWeight * regularizerFull,  options.regularizerWeight * regularizerReduced,0)

    @printf("\n
    J values
    ---------------------------------
    full registration     %3.3e
    reduced operator      %3.3e
    gradient descent      %3.3e\n", ssdFineFull + options.regularizerWeight * regularizerFull, ssdFineReduced + options.regularizerWeight * regularizerReduced,0)

    return fineDisplacementFull, fineDisplacementReduced

end

##
close("all")
fineDisplacementFull, fineDisplacementReduced = registerFineGlobal(plot=true, showDeformationImage=false, reg=createDiffusiveOperatorDirichlet)



##
