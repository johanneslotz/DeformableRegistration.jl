using DeformableRegistration
using Images
# using Logging
import DeformableRegistration.Distance.ngfDistance


# Normalized gradient field distance
function ngfDistance(referenceImage::regImage,templateImage::regImage,
                     transformedGrid::Array{Float64,1};
                     doDerivative::Bool=false,doHessian::Bool=false,options::regOptions=regOptions(), centeredGrid::Array{Float64,1}=zeros(1)
                     )

    opt = options

    edgeParameterR=options.edgeParameterR
    edgeParameterT=options.edgeParameterT
    debug = false
    useEdgeParameterInNumerator=options.useEdgeParameterInNumerator
    parametricOnly=options.parametricOnly

    if(size(centeredGrid)[1]==1)
      centeredGrid = getCellCenteredGrid(referenceImage).data
    end

    if(checkForOddNumberOfGridPoints(referenceImage.data,transformedGrid))
      error("StaggeredGrids are not supported. Please transform to cell centered first.")
    end


    #c,dTtuple = linearInter2D(templateImage.data, ΩT, mT, deformedGrid,doDerivative=doDerivative)
    if doDerivative
      transformedImage, dX_transformedImage, dY_transformedImage =
          interpolateImage(templateImage,transformedGrid,doDerivative=true)
      dT = spdiagm((dX_transformedImage, dY_transformedImage),[0,prod(size(transformedImage))])

    else
      transformedImage, dX_transformedImage, dY_transformedImage =
          interpolateImage(templateImage,transformedGrid,doDerivative=false)
    end

    Logging.debug("\n\n edgeParameterR: %5e    edgeParameterT %5e \n", edgeParameterR, edgeParameterT)

    m = size(referenceImage.data)
    h = referenceImage.voxelsize
    Rc = Array(referenceImage.data)
    shift = referenceImage.shift

    shortFiniteDiffX = spdiagm((-ones(m[1]-1,1),ones(m[1]-2,1)),[0, 1],m[1]-1,m[1])/(2*h[1])
    #shortFiniteDiffX[1] = shortFiniteDiffX[2]
    #shortFiniteDiffX[end] = shortFiniteDiffX[end-1]

    shortFiniteDiffY = spdiagm((-ones(m[2]-1,1),ones(m[2]-2,1)),[0, 1],m[2]-1,m[2])/(2*h[2])
    #shortFiniteDiffY[1] = shortFiniteDiffY[2]
    #shortFiniteDiffY[end] = shortFiniteDiffY[end-1]

    averageX = spdiagm((0.5*ones(m[1]-2,1),0.5*ones(m[1]-1,1)),[-1, 0],m[1],m[1]-1)
    averageY = spdiagm((ones(m[2]-2,1)*0.5,ones(m[2]-1,1)*0.5),[-1, 0],m[2],m[2]-1)

    G1 = sparse(kron(speye(m[2]),shortFiniteDiffX));
    G2 = sparse(kron(shortFiniteDiffY,speye(m[1])));

    AvgX = sparse(kron(speye(m[2]),averageX));
    AvgY = sparse(kron(averageY,speye(m[1])));

    gradTx  = G1*transformedImage[:]
    gradTy  = G2*transformedImage[:]

    gradRx  = G1*Rc[:]
    gradRy =  G2*Rc[:]

    lengthGT = sqrt.(sum(AvgX * (gradTx.^2) + AvgY * (gradTy.^2), 2) + edgeParameterT^2) #epsilon-norm
    lengthGR = sqrt.(sum(AvgX * (gradRx.^2) + AvgY * (gradRy.^2), 2) + edgeParameterR^2)

    if useEdgeParameterInNumerator
      r1  = sum(AvgX * (gradRx.*gradTx) +
		AvgY * (gradRy.*gradTy) ,2) .+ edgeParameterR*edgeParameterT    # numerator
    else
      r1  = sum(AvgX * (gradRx.*gradTx) +
		AvgY * (gradRy.*gradTy) ,2)    # numerator
    end
    r2  = 1./(lengthGT.*lengthGR)             #  ... and denominator
    rc  = r1 .* r2                            # combine things and finalize

    dFunctionValue = 0
    drc = 0

    if doDerivative
      dr1 = (AvgX * spdiagm(gradRx[:]) * G1) + (AvgY * spdiagm(gradRy[:]) * G2) #(AvgX * gradRx) * AvgX * G1 + (AvgY * gradRy) * AvgX * G2
      dr2 = ( spdiagm(-r2.^2)
              * spdiagm(lengthGR)
              * (spdiagm(1./(2 .* lengthGT)))
              * ((2 * AvgX * spdiagm(gradTx) * G1) + (2 * AvgY * spdiagm(gradTy) * G2))
            )

      if(parametricOnly)
          N = prod(m)
          Q = sparse([centeredGrid[1:N] centeredGrid[N+1:end] ones(N)])
          Q = [Q spzeros(size(Q)[1],size(Q)[2]); spzeros(size(Q)[1],size(Q)[2]) Q]
          dT = dT * Q
      end
      #dr2Partial = reshape(  r1 .* -1 ./ (lengthGR .* lengthGT.^3)  ,(size(r1)))

      #drc = (    spdiag((r2 .* (AvgX * gradRx) + dr2Partial .* (AvgX * gradTx) )[:]) * AvgX * G1
      #         + spdiag((r2 .* (AvgY * gradRy) + dr2Partial .* (AvgY * gradTy) )[:]) * AvgY * G2 )
      drc = spdiagm(r1)*dr2 + spdiagm(r2) * dr1
      drc_times_dT = drc*dT;
      dFunctionValue  = Array(-2.0*prod(h)*rc'*drc_times_dT)

      if doHessian
        d2FunctionValue(x) =  (2*prod(h) .* drc_times_dT' * drc_times_dT) * x   # note the missing minus sign is not a bug!
			else
			  d2FunctionValue = 0
		  end
	else
		d2FunctionValue = 0
	end

	functionValue  = (length(rc)*prod(h) - prod(h) * (rc'*rc))[1]

	if(ndims(dFunctionValue)>1)
		dFunctionValue = vec(dFunctionValue)
	end
	#return functionValue[1], dFunctionValue', d2FunctionValue, r1, 1.-(rc.^2)
	return [functionValue,dFunctionValue,d2FunctionValue] #, drc

end
