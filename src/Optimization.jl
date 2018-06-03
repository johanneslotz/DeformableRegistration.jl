module Optimization

using MicroLogging
using KrylovMethods

using DeformableRegistration
include("helpers/objectiveFunctionCreation.jl")
include("../experimental/aspin.jl")

export checkStoppingCriteria, ArmijoLineSearch, optimizeGaussNewton, optimizeGaussNewtonAugmentedLagrangian


# function augmentedLagrangian(x, λ, μ, c::Function)
#     m = length(x)
#     F = - λ'*c(x)[1] + μ/2 * c(x)[1]'*c(x)[1]
#     dF = - c(x)[2]*λ   + μ * c(x)[2]'*c(x)[1]
#     d2F = spzeros(m,m)#pdiagm(-λ, 0)  +  μ * c(x)[2]'*c(x)[2] + μ * spdiagm(c(x)[1],0)
#     return [F, dF, d2F]
# end

function augmentedLagrangian(x, λ, μ, c::Function)
    m = length(x)
    cx, ∇c = c(x)
    F = - λ' * cx + μ/2 * cx' * cx
    ∇F = - ∇c * λ   + μ * ∇c' * cx
    ∇2F(y) = μ * ∇c' * ∇c * y
    return [F, ∇F, ∇2F]
end

"""
optimizeGaussNewtonAugmentedLagrangian(Jfunc, y, yReference, options; constraint)

Augmented lagrangian implementation following Nocelda & Wright, 2006, p. 515 (Framework 17.3)
"""
function optimizeGaussNewtonAugmentedLagrangian(Jfunc::Function,  # objective Function
                             y::Array{Float64,1}, yReference::Array{Float64,1}, options::regOptions;
                             constraint::Function = x -> [0, 0, 0], printFunction=x->[], gradientDescentOnly=false, displayCallback=x->x) #no constraint by default,

    c = constraint
    computeConstraint = ! (c == x -> [0, 0, 0])

    L(x, λ, μ; doDerivative=true, doHessian=true) = (c(y) != [0, 0, 0]) ? Jfunc(x, doDerivative=doDerivative, doHessian=doHessian) +
     augmentedLagrangian(x, λ, μ, c) : Jfunc(x, doDerivative=doDerivative, doHessian=doHessian)

    λ = spzeros(size(y)...)
    μ = 1

    JRef = L(yReference, λ, μ)[1]
    JOld = NaN
    oldNormCy = 99 * norm(c(y)[1], Inf)


    for iter = 1:options.maxIterGaussNewton
     if computeConstraint
         λ = λ - μ * c(y)[1]

         "Enforce reduction of infeasability"
         τReduction = 0.8
         if norm(c(y)[1], Inf) >  τReduction * oldNormCy;
             @debug @sprintf("... norm(c(y)[1], Inf) >  %1.1f * oldNormCy: %3.3e > %3.3e", τReduction, norm(c(y)[1], Inf), oldNormCy)
             @debug "... updating μ = 2 * μ = ", 2*μ
             μ = 2 * μ;
         end
         oldNormCy = norm(c(y)[1], Inf)

     end

     J, dJ, d2J = L(y, λ, μ, doDerivative=true, doHessian=true)

     @debug d2J(-dJ)
     if gradientDescentOnly
         dy = -dJ
         cgIterations = -1
     else
         dy,flag,resvec,cgIterations = KrylovMethods.cg(d2J,-dJ,maxIter=options.maxIterCG, tol=1e-5)[1:4] #2/3 of the time of this function is here
     end

        # check descent direction
        if( (dJ'*dy)[1] > 0)
         dy = -dy
         warn("Changing sign of computed descent direction. This is a suspicious move.")
        end

     # armijo line search (LS) method
     stepLength,LSiter,LSfailed = ArmijoLineSearch(x -> L(x, λ, μ), J, dJ, y, dy, tolLS=1e-4)

     if(LSfailed)
         info("STOPPING after line search failed.")
         break
     end

     # output
     # is this a IJulia notebook? -> delete previous iteration's output
     if isdefined(Main, :IJulia) && Main.IJulia.inited
        displayCallback(true)
        @info size(y)
     end
     D, α, S = printFunction(y)
     s = @sprintf("%03d | D %2.2e | S %2.2e | λ!=0: %d | |c(y)|= %2.2e | μ=%2.2e | LSiter: %2d | CGiter: %3d | J/Jref: %2.2e",
         iter, D[1], S[1], sum(λ.!=0), norm(c(y)[1], Inf), μ,  LSiter, cgIterations, J[1]/JRef[1],
         )
     @info(s)
     @debug "  λ=", λ
     @debug "  c=", c(y)[1]


     # update parameter y
     y = y + stepLength .* dy

     # stopping criteria
     if (iter>1)
         if (checkStoppingCriteria(J, JOld, JRef, dJ, y, stepLength*dy,
                 tolJ = options.stopping["tolJ"],    # tolerance: change of the objective function
                 tolY = options.stopping["tolY"],    # tolerance: change of the variables
                 tolG = options.stopping["tolG"],    # tolerance: change of the gradient
                 tolQ = options.stopping["tolQ"]))   # tolerance: change of quotient)
             break
         end
     end

     # updating JOld
     JOld = J[1]
     info("\n")
    end

    return y

end



function optimizeGaussNewton(Jfunc::Function,  # objective Function
                             y::Array{Float64,1}, yInitial::Array{Float64,1}, options::regOptions)

    # Initilization of JRef: reference value of objective Function
    JRef = Jfunc(yInitial)[1]
    JOld = NaN

    for iter = 1:options.maxIterGaussNewton

        # get derivatives of objective function
        J, dJ, d2J = Jfunc(y,doDerivative=true,doHessian=true)

        # solve system of linear equations d2J*dy=-dJ
          # (dy: change of the variables / search direction)
        # number of parameters < 10 => use backslash operator
        # else (large systems) => use conjugate gradient (cg) method

        cgIterations=0

        if(length(y)<10 && !(typeof(d2J) <: Function))
            dy=d2J\-dJ
        else
            dy,flag,resvec,cgIterations = KrylovMethods.cg(d2J,-dJ,maxIter=options.maxIterCG, tol=1e-5)[1:4]
        end

	    # check descent direction
	    if( (dJ'*dy)[1] > 0)
            dy = -dy
            warn("Changing sign of computed descent direction. This is a suspicious move.")
	    end

        # armijo line search (LS) method
        stepLength,LSiter,LSfailed = ArmijoLineSearch(Jfunc,J,dJ,y,dy)
        if(LSfailed)
            break
        end
        displayCallback()
        # output
        if(cgIterations==0)
          s = @sprintf("%3d: J %8.4e     LSiter: %2d    J/Jref: %1.2f \n",iter, J[1], LSiter, J[1]/JRef[1])
          info(s)
        else
          s = @sprintf("%3d: J %8.4e     LSiter: %2d     CGiter: %3d     J/Jref: %1.2f \n",iter, J[1], LSiter,cgIterations, J[1]/JRef[1])
          info(s)
        end

        # update parameter y
        y = y + stepLength.*dy

        # stopping criteria
        if (iter>1)
            if checkStoppingCriteria(J, JOld, JRef, dJ, y, stepLength*dy,
                    tolJ = options.stopping["tolJ"],    # tolerance: change of the objective function
                    tolY = options.stopping["tolY"],    # tolerance: change of the variables
                    tolG = options.stopping["tolG"],    # tolerance: change of the gradient
                    tolQ = options.stopping["tolQ"])    # tolerance: change of quotient)
                break
            end
        end

        # updating JOld
        JOld = J[1]

    end

    return y

end

function checkArmijoStepValid(Jfunc::Function, Jc::Float64, y::Array{Float64,1},
        stepLength::Float64, dy::Array{Float64,1}, tolLS::Float64,
        dJ::Array{Float64,1})
    @debug "Jfunc(y+stepLength*dy)[1] < (Jc + (stepLength.*tolLS.*(dJ'*dy))[1])" Jfunc(y+stepLength*dy)[1]  (Jc + (stepLength.*tolLS.*(dJ'*dy))[1])
    validStep = Jfunc(y+stepLength*dy)[1] < (Jc + (stepLength.*tolLS.*(dJ'*dy))[1])
    return validStep
end

function ArmijoLineSearch(Jfunc::Function,         # objective function
                          Jc::Float64,             # Jc = Jfunc(y)
                          dJ::Array{Float64,1},    # gradient of the objective funtion
                          y::Array{Float64,1},     # old variables
                          dy::Array{Float64,1};    # change of the variables / search direction
                          tolLS::Float64 = 1e-4)   # tolerance

    stepLength = 1.0; LSiter = 1; LSfailed = false;

    while( true )
        validStep = checkArmijoStepValid(Jfunc, Jc, y, stepLength, dy, tolLS, dJ)
        if validStep
            break
        end
        if(stepLength < 1e-4)
            LSfailed = true
            break
        end

        stepLength = 0.5 * stepLength
        LSiter = LSiter + 1

    end

    if(LSfailed)
        s = @sprintf("   Line search failed after %2d iterations.",LSiter)
        @info s
    end

    return stepLength,LSiter,LSfailed

end

# somehow slower with this signature
# function checkStoppingCriteria(J::Float64,JOld::Float64,JRef::Float64,    # value of the current objective function J(y+dy), the old J(y) and the reference J(yInital)
#                                dJ::Array{Float64,1},             # gradient of the objective funtion
#                                y::Array{Float64,1},              # old variables
#                                dy::Array{Float64,1};             # change of the variables / search direction
#                                tolJ::Float64 = 1e-3,    # tolerance: change of the objective function
#                                tolY::Float64 = 1e-2,    # tolerance: change of the variables
#                                tolG::Float64 = 1e-2,    # tolerance: change of the gradient
#                                tolQ::Float64 = 1e-4     # tolerance: change of quotient
#                                )

function checkStoppingCriteria(J,JOld,JRef,    # value of the current objective function J(y+dy), the old J(y) and the reference J(yInital)
                               dJ,             # gradient of the objective funtion
                               y,              # old variables
                               dy;             # change of the variables / search direction
                               tolJ = 1e-3,    # tolerance: change of the objective function
                               tolY = 1e-2,    # tolerance: change of the variables
                               tolG = 1e-2,    # tolerance: change of the gradient
                               tolQ = 1e-4     # tolerance: change of quotient
                               )


    # Initialize STOP-Array with false
    STOP = Array{Bool}(5); STOP[:] = false

    # Either STOP[1:3] or STOP[4]
    STOP[1] = abs(JOld-J) <= tolJ  * (1+abs(JOld))
    STOP[2] = norm(dy)    <= tolY  * (1+norm(y))
    STOP[3] = norm(dJ)    <= tolG  * (1+abs(JOld))
    STOP[4] = norm(dJ)    <= 1e6 * eps()
    STOP[5] = abs((JOld-J)/(JRef-J)) <= tolQ

    if(all(STOP[1:3]))
        s = @sprintf("STOPPING CRITERIA:\n")
        @info s
        s = @sprintf("   1. |JOld-J| = %3e <= %3e (tolJ * (1 + |JOld|)) \n",abs(JOld-J), tolJ*(1+abs(JOld)) )
        @info s
        s = @sprintf("   2. ||dy||   = %3e <= %3e (tolY * (1 + ||y||))\n",norm(dy),    tolY*(1+norm(y))  )
        @info s
        s = @sprintf("   3. ||dJ||   = %3e <= %3e (tolG * (1 + |JOld|)) \n",norm(dJ),    tolG*(1+abs(JOld)) )
        @info s
    end

    if(STOP[4])
        s = @sprintf("STOPPING CRITERION:\n")
        @info s
        s = @sprintf("   4. ||dJ|| = %3e <= %3e (1e6*eps)\n",norm(dJ), 1e6*eps())
        @info s
    end

    if(STOP[5])
        s = @sprintf("STOPPING CRITERION:\n")
        @info s
        s = @sprintf("   5. |(JOld-J)/(JRef-J)|= %3e <= %3e (tolQ)\n",abs((JOld-J)/(JRef-J)), tolQ)
        @info s
    end

    return all(STOP[1:3]) | STOP[4] | STOP[5]

end

end
