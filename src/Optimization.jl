module Optimization

using Logging
using KrylovMethods

using ImageRegistration

export checkStoppingCriteria, ArmijoLineSearch, optimizeGaussNewton

function optimizeGaussNewton(Jfunc::Function,
                             y::Array{Float64,1}, options::regOptions)

    JRef = Jfunc(y)[1]; y0=y; JOld = Inf

    output = (Logging.LogLevel == Logging.DEBUG) |  (Logging.LogLevel == Logging.INFO)

    for iter = 1:options.maxIterGaussNewton

        # get derivatives of objective function
        J, dJ, d2J = Jfunc(y,doDerivative=true,doHessian=true)

        # solve system of linear equations d2J*dy=-dJ
        # number of parameters < 10 => use backslash operator
        # else use cg method
        cgIterations=0
        Logging.debug("dJ: ",dJ)
        Logging.debug("d2J: ",d2J)

        if(length(y)<10 && (typeof(d2J)!=Function))
            dy=d2J\-dJ
        else
            dy,flag,resvec,cgIterations = KrylovMethods.cg(d2J,-dJ,maxIter=options.maxIterCG)[1:4]
        end

		    # check descent direction
		    if( (dJ'*dy)[1] > 0)
            dy = -dy
            Logging.warn("Changing sign of computed descent direction. This is a suspicious move.")
		    end

        # armijo line search method
        stepLength,LSiter,LSfailed = ArmijoLineSearch(Jfunc,J,dJ,y,dy,printFailedLineSearch=output)
        if(LSfailed)
            break
        end

        # output

        if(cgIterations==0)
          s = @sprintf("%3d: J %8.4e     LSiter: %2d    J/Jref: %1.2f \n",iter, J[1], LSiter, J[1]/JRef[1])
          Logging.info(s)
        else
          s = @sprintf("%3d: J %8.4e     LSiter: %2d     CGiter: %3d     J/Jref: %1.2f \n",iter, J[1], LSiter,cgIterations, J[1]/JRef[1])
          Logging.info(s)
        end

        # update parameter y
        y = y + stepLength.*dy
        Logging.debug("y: ",y)

        # stopping criteria
        if checkStoppingCriteria(J[1],JOld[1],JRef[1],dJ,y0,y,stepLength*dy,printActiveStoppingCirteria=output)
            break
        end

        # updating JOld
        JOld = J[1]

    end

    return y

end

function ArmijoLineSearch(Jfunc::Function,         # objective function
                          Jc::Float64,             # Jc = Jfunc(y)
                          dJ::Array{Float64,1},    # gradient of the objective funtion
                          y::Array{Float64,1},     # old variables
                          dy::Array{Float64,1};    # search direction
                          tolLS::Float64 = 1e-4,
                          printFailedLineSearch = false)

    stepLength = 1.0; LSiter = 1; LSfailed = false;

    while( (Jfunc(y+stepLength*dy)[1]) >= (Jc + (stepLength.*tolLS.*dJ'*dy)[1]) )

        if(stepLength < 1e-5)
            LSfailed = true
            break
        end

        stepLength = 0.5 * stepLength
        LSiter = LSiter + 1

    end

    if(printFailedLineSearch & LSfailed)
        @printf("STOPPING:\n")
        @printf("   0. Line search failed after %2d iterations.\n",LSiter)
    end

    return stepLength,LSiter,LSfailed

end

function checkStoppingCriteria(J,JOld,JRef,    # value of the current objective function J(y+dy), the old J(y) and the reference J(y0)
                               dJ,             # gradient of the objective funtion
                               y0,             # intial guess
                               y,              # old variables
                               dy;             # change of the variables / search direction
                               tolJ = 1e-3,    # tolerance: change of the objective function
                               tolY = 1e-2,    # tolerance: change of the variables
                               tolG = 1e-2,    # tolerance: change of the gradient
                               printActiveStoppingCirteria = false)

    STOP = Array(Bool,4); STOP[:] = false

    STOP[1] = abs(JOld-J) <= tolJ  * (1+abs(JRef))
    STOP[2] = norm(dy)    <= tolY  * (1+norm(y0))
    STOP[3] = norm(dJ)    <= tolG  * (1+abs(JRef))
    STOP[4] = norm(dJ)    <= 1e6 * eps()

    if( printActiveStoppingCirteria & all(STOP[1:3]) )
	      @printf("STOPPING CRITERIA:\n")
        @printf("   1. |JOld-J| = %5e <= %5e \n",abs(JOld-J), tolJ*(1+abs(JRef)) )
	      @printf("   2. norm(dy) = %5e <= %5e \n",norm(dy),    tolY*(1+norm(y0))  )
	      @printf("   3. norm(dJ) = %5e <= %5e \n",norm(dJ),    tolG*(1+abs(JRef)) )
    end
    if( printActiveStoppingCirteria & STOP[4] )
	      @printf("STOPPING CRITERIA:\n")
        @printf("   4. norm(dJ) = %5e <= %5e \n",norm(dJ), 1e6*eps())
    end

    return all(STOP[1:3]) | STOP[4]

end

end
