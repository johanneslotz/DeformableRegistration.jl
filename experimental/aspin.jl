function optimizeGaussNewtonASPIN(Jfunc::Function, S::Array{Float64,2}, # objective Function, subspaces
                             y::Array{Float64,1}, yInitial::Array{Float64,1}, options::regOptions)

    # Initilization of JRef: reference value of objective Function
    JRef = Jfunc(yInitial)[1]
    JOld = NaN

    for iter = 1:options.maxIterGaussNewton

        # get derivatives of objective function

        # solve system of linear equations d2J*dy=-dJ
          # (dy: change of the variables / search direction)
        # number of parameters < 10 => use backslash operator
        # else (large systems) => use conjugate gradient (cg) method

        cgIterations=0

        yp = repmat(y',size(S,1))
        print(size(yp))
        # S: each row one index vector
        for iterPatch = 1:size(S,1)
            for iterPatchGN = 1:3
                Jp, dJp, d2Jp = Jfunc(yp[iterPatch,:] ,doDerivative=true,doHessian=true)
                SiXY = repmat(S[iterPatch,:],2)
                dJp = dJp[SiXY .== 1]
                # Hi(x) = repmat(S[iterPatch,:],2) .* d2Jp(repmat(S[iterPatch,:],2) .* x)
                function embedInZeros(x,sz, idx)
                    y = zeros(sz)
                    y[idx[:] .==  1] = x
                    return y
                end
                Hi(x) = d2Jp( embedInZeros(x, size(SiXY), SiXY) )[SiXY .== 1]
                dy,flag,resvec,cgIterations = KrylovMethods.cg(Hi,dJp,maxIter=options.maxIterCG, tol=1e-5)[1:4]
                if( (dJp'*dy)[1] > 0)
                    dy = -dy
                    warn("Changing sign of computed descent direction. This is a suspicious move.")
                end
                # stepLength,LSiter,LSfailed = ArmijoLineSearch(Jfunc,Jp,dJp,yp[iterPatch,SiXY .== 1],dy)
                # if(LSfailed)
                #     info("Line search failed.")
                #     break
                # end
                stepLength = 0.1
                LSiter = -1

                # output
                if(cgIterations==0)
                  s = @sprintf("PATCH: %3d: J %8.4e     LSiter: %2d    J/Jref: %1.2f \n",iterPatchGN, Jp[1], LSiter, Jp[1]/JRef[1])
                  info(s)
                else
                  s = @sprintf("PATCH: %3d: J %8.4e     LSiter: %2d     CGiter: %3d     J/Jref: %1.2f \n",iterPatchGN, Jp[1], LSiter,cgIterations, Jp[1]/JRef[1])
                  info(s)
                end

                # update parameter y
                yp[iterPatch,SiXY .== 1] = yp[iterPatch,SiXY .== 1] + stepLength.*dy

                # stopping criteria
                if (iter>1)
                    if checkStoppingCriteria(Jp, JOld, JRef, dJp, yp[iterPatch,SiXY .== 1],  stepLength.*dy,
                            tolJ = options.stopping["tolJ"],    # tolerance: change of the objective function
                            tolY = options.stopping["tolY"],    # tolerance: change of the variables
                            tolG = options.stopping["tolG"],    # tolerance: change of the gradient
                            tolQ = options.stopping["tolQ"])    # tolerance: change of quotient)
                        break
                    end
                end
                # JOld = J[1]

            end
        end

        y[:] = 0 * y[:]
        for i = 1:size(yp,1)
            y[:] = y[:] + yp[i, :][:]
        end

        J, dJ, d2J = Jfunc(y,doDerivative=true,doHessian=true)

        println("TODO: finish ASPIN implementation.")

        # function aspind2J(x)
        #     d2(x) = d2J(x)
        #     for i = 1:size(yp,1)
        #         d2(x) = d2(x)  repmat(S[i,:],2)
        #     end
        # end

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
            info("Line search failed.")
            break
        end

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
