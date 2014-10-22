module Optimization

using KrylovMethods

export checkStoppingCriteria, ArmijoLineSearch, optimizeGaussNewton

function optimizeGaussNewton(Jfunc::Function,
                             JfuncWithDerivative::Function,
                             y::Array{Float64,1};
                             maxIterGaussNewton=10,
                             maxIterCG=2000,
                             output=false)

    JRef = Jfunc(y)[1]; y0=y; JOld = Inf

    for iter = 1:maxIterGaussNewton

        # get derivatives of objective function
        J, dJ, d2J = JfuncWithDerivative(y)

        # solve system of linear equations d2J*dy=-dJ
        # number of parameters < 10 => use backslash operator
        # else use cg method
        cgIterations=0
        if(length(y)<10)
            dy=d2J\-dJ
        else
            dy,flag,resvec,cgIterations = KrylovMethods.cg(d2J,-dJ,maxIter=maxIterCG)[1:4]
        end

		    # check descent direction
		    if( (dJ'*dy)[1] > 0)
            dy = -dy
		    end

        # armijo line search method
        stepLength,LSiter,LSfailed = ArmijoLineSearch(Jfunc,dJ,y,dy,printFailedLineSearch=output)
        if(LSfailed)
            break
        end

        # output
        if(output)
            if(cgIterations==0)
                @printf("%3d: J %8.4e     LSiter: %2d    J/Jref: %1.2f \n",iter, J[1], LSiter, J[1]/JRef[1])
            else
                @printf("%3d: J %8.4e     LSiter: %2d     CGiter: %3d     J/Jref: %1.2f \n",iter, J[1], LSiter,cgIterations, J[1]/JRef[1])
            end
        end

        # update parameter y
        y = y + stepLength.*dy

        # stopping criteria
        if checkStoppingCriteria(J[1],JOld[1],JRef[1],dJ,y0,y,stepLength*dy,printActiveStoppingCirteria=output)
            break
        end

        # updating JOld
        JOld = J[1]

    end

    return y

end

function ArmijoLineSearch(J,   # objective function
                          dJ,  # gradient of the objective funtion
                          y,   # old variables
                          dy;  # search direction
                          tolLS = 1e-4,
                          printFailedLineSearch = false)

    stepLength = 1.0; LSiter = 1; LSfailed = false; Jc = J(y)[1]

    while( (J(y+stepLength*dy)[1]) >= (Jc + (stepLength.*tolLS.*dJ'*dy)[1]) )

        if(stepLength < 1e-5)
            LSfailed = true
            break
        end

        stepLength = stepLength / 2
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
