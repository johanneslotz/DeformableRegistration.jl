using Images
 using Logging
 using DeformableRegistration: ImageProcessing, Transformation, Interpolation, Optimization
 using DeformableRegistration.Regularizer
 using DeformableRegistration.regOptions
 using Base.Test
 using DeformableRegistration.Optimization.augmentedLagrangian
 #Logging.configure(level=Logging.DEBUG)
##

@testset "constraint derivatives" begin
    include("../test/helpers/checkDerivative.jl")


    function zerosConstraintAtIndex(x::Array{Float64,1}, index)
        f = 0.5 * x .* x
        f[.!(index.==1)] = 0
        dF = spdiagm(x .* index,0)
        return [f, dF]
    end
    c(x) = zerosConstraintAtIndex(x, index)

    function J(x; doDerivative=true, doHessian=true)
        r = (x - Array(1:10))
        f = 0.5 * (r'*r)
        df = r
        d2F = speye(10)
        return [f, df, d2F]
    end

    xtest = 10*rand(10)
    index = zeros(size(xtest)); index[1:3] = 1
    Ifunc(x) = zerosConstraintAtIndex(x,index)
    errlin,errquad = checkDerivative(Ifunc,xtest)
    @test checkErrorDecay(errquad)

    xtest = rand(10)
    errlin,errquad = checkDerivative(J,xtest)
    @test checkErrorDecay(errquad)


    λ = 10*[-20, -30, -50, 0,0,0,0,0,0,0]
    μ = 32
    c(x) = zerosConstraintAtIndex(x, index)
    testf(x) = augmentedLagrangian(x, λ, μ, c)
    errlin,errquad = checkDerivative(testf,xtest)
    @test checkErrorDecay(errquad)

end

@testset "constrainmts are met" begin
    x = ones(10)*2;
    index = zeros(size(x))
    index[1:3] = 1
    #x[index.==1]=0
    function zerosConstraintAtIndex(x::Array{Float64,1}, index)
        f = 0.5 * x .* x
        f[.!(index.==1)] = 0
        dF = spdiagm(x .* index,0)
        return [f, dF]
    end
    c(x) = zerosConstraintAtIndex(x, index)

    function J(x; doDerivative=true, doHessian=true)
        r = (x - Array(1:10))
        f = 0.5 * (r'*r)
        df = r
        d2F = speye(10)
        return [f, df, d2F]
    end

    options = regOptions()
    options.maxIterCG =100
    options.maxIterGaussNewton =100
    options.stopping["tolQ"]= 1e-20
    yopt = optimizeGaussNewtonAugmentedLagrangian(J, x, options, constraint=c, constraintObjectiveFunction=augmentedLagrangian)
    debug(@sprintf("... constraint norm:         %1.3e", norm(c(yopt)[1])))
    debug(@sprintf("... objective function norm: %1.3e", J(yopt)[1]))
    y_true = Array(1:10)
    y_true[index.==1] = 0
    debug(         "... error:                   ",yopt - y_true)
    debug(@sprintf("... mean error:              %1.3e",mean(abs.(yopt - y_true))))
    debug(@sprintf("... error norm:              %1.3e",norm(yopt - y_true)))
    @test norm(yopt - y_true)^2 < 0.15
end
