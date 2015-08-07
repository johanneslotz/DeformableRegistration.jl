module Optimization

using Logging
using KrylovMethods

using ImageRegistration

export checkStoppingCriteria, ArmijoLineSearch, optimizeGaussNewton

function optimizeGaussNewton(Jfunc::Function,  # objective Function
                             y::Array{Float64,1}, yInitial::Array{Float64,1}, options::regOptions)

    # Initilization of JRef: reference value of objective Function
    JRef = Jfunc(yInitial)[1]


    output = (Logging.LogLevel == Logging.DEBUG) |  (Logging.LogLevel == Logging.INFO)

    for iter = 1:options.maxIterGaussNewton

        # get derivatives of objective function
        J, dJ, d2J = Jfunc(y,doDerivative=true,doHessian=true)

        # solve system of linear equations d2J*dy=-dJ
          # (dy: change of the variables / search direction)
        # number of parameters < 10 => use backslash operator
        # else (large systems) => use conjugate gradient (cg) method

        cgIterations=0

        if(length(y)<10 && (typeof(d2J)!=Function))
            dy=d2J\-dJ
        else
            dy,flag,resvec,cgIterations = KrylovMethods.cg(d2J,-dJ,maxIter=options.maxIterCG, tol=1e-5)[1:4]
        end

		    # check descent direction
		    if( (dJ'*dy)[1] > 0)
            dy = -dy
            Logging.warn("Changing sign of computed descent direction. This is a suspicious move.")
		    end

        # armijo line search (LS) method
        stepLength,LSiter,LSfailed = ArmijoLineSearch(Jfunc,J,dJ,y,dy)
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
        #Logging.debug("y: ",y)

        # stopping criteria
        if ((iter>1) && checkStoppingCriteria(J[1],JOld[1],JRef[1],dJ,y,stepLength*dy)) #printActiveStoppingCirteria=output
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
                          dy::Array{Float64,1};    # change of the variables / search direction
                          tolLS::Float64 = 1e-4)   # tolerance

    stepLength = 1.0; LSiter = 1; LSfailed = false;

    while( (Jfunc(y+stepLength*dy)[1]) >= (Jc + (stepLength.*tolLS.*dJ'*dy)[1]) )

        if(stepLength < 1e-5)
            LSfailed = true
            break
        end

        stepLength = 0.5 * stepLength
        LSiter = LSiter + 1

    end

    if(LSfailed)
        s = @sprintf("   0. Line search failed after %2d iterations.",LSiter)
        Logging.debug(s)
    end

    return stepLength,LSiter,LSfailed

end

function checkStoppingCriteria(J,JOld,JRef,    # value of the current objective function J(y+dy), the old J(y) and the reference J(yInital)
                               dJ,             # gradient of the objective funtion
                               y,              # old variables
                               dy;             # change of the variables / search direction
                               tolJ = 1e-3,    # tolerance: change of the objective function
                               tolY = 1e-2,    # tolerance: change of the variables
                               tolG = 1e-2,    # tolerance: change of the gradient
                               tolQ = 1e-4,    # tolerance: change of quotient
                               printActiveStoppingCirteria = true)


    # Initialize STOP-Array with false
    STOP = Array(Bool,5); STOP[:] = false

    # Either STOP[1:3] or STOP[4]
    STOP[1] = abs(JOld-J) <= tolJ  * (1+abs(JOld))
    STOP[2] = norm(dy)    <= tolY  * (1+norm(y))
    STOP[3] = norm(dJ)    <= tolG  * (1+abs(JOld))
    STOP[4] = norm(dJ)    <= 1e6 * eps()
    STOP[5] = abs((JOld-J)/(JRef-J)) <= tolQ


    if( printActiveStoppingCirteria & all(STOP[1:3]) )
	    s = @sprintf("STOPPING CRITERIA:\n")
      Logging.info(s)
      s = @sprintf("   1. |JOld-J| = %5e <= %5e \n",abs(JOld-J), tolJ*(1+abs(JOld)) )
      Logging.info(s)
	    s = @sprintf("   2. ||dy|| = %5e <= %5e \n",norm(dy),    tolY*(1+norm(y))  )
      Logging.info(s)
	    s = @sprintf("   3. ||dJ|| = %5e <= %5e \n",norm(dJ),    tolG*(1+abs(JOld)) )
      Logging.info(s)
    end
    if( printActiveStoppingCirteria & STOP[4] )
	     s = @sprintf("STOPPING CRITERION:\n")
       Logging.info(s)
       s = @sprintf("   4. ||dJ|| = %5e <= %5e \n",norm(dJ), 1e6*eps())
       Logging.info(s)
    end

   if( printActiveStoppingCirteria & STOP[5] )

	       s = @sprintf("STOPPING CRITERION:\n")
         Logging.info(s)
         s = @sprintf("   5. |(JOld-J)/(JRef-J)|= %5e <= %5e \n",abs((JOld-J)/(JRef-J)), tolQ)
         Logging.info(s)
     end

    return all(STOP[1:3]) | STOP[4] | STOP[5]

end

end
