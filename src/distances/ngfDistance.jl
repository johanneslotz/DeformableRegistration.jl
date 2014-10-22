
function ngfDistance(referenceImage::Image,templateImage::Image,
                     transformedGrid::Array{Float64,1};
		                 doDerivative=false,doHessian=false,
		                 edgeParameterR=1,edgeParameterT=1,
		                 debug=false,useEdgeParameterInNumerator=true)


    if(checkStaggered(referenceImage,transformedGrid))
      error("StaggeredGrids are not supported. Please transform to cell centered first.")
    end


    #c,dTtuple = linearInter2D(templateImage.data, ΩT, mT, deformedGrid,doDerivative=doDerivative)
    if doDerivative
      transformedImage, dY_transformedImage, dX_transformedImage =
        linearImageInterpolationAtGridWithDerivative(templateImage,transformedGrid)
      dT = spdiagm((dY_transformedImage, dX_transformedImage),[0,prod(size(transformedImage))])

    else
      transformedImage =
        linearImageInterpolationAtGrid(templateImage,transformedGrid)
    end



    if debug==true
      @printf("\n\n edgeParameterR: %5e    edgeParameterT %5e \n", edgeParameterR, edgeParameterT)
    end

    m = size(referenceImage)
    h = referenceImage.properties["pixelspacing"]
    Rc = referenceImage.data
    ΩR = referenceImage.properties["spatialdomain"]

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

    lengthGT = sqrt(sum(AvgX * (gradTx.^2) + AvgY * (gradTy.^2), 2) + edgeParameterT^2) #epsilon-norm
    lengthGR = sqrt(sum(AvgX * (gradRx.^2) + AvgY * (gradRy.^2), 2) + edgeParameterR^2)

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
    d2FunctionValue = 0
    drc = 0

    if doDerivative
      dr1 = (AvgX * spdiagm(gradRx[:]) * G1) + (AvgY * spdiagm(gradRy[:]) * G2) #(AvgX * gradRx) * AvgX * G1 + (AvgY * gradRy) * AvgX * G2
      dr2 = ( spdiagm(-r2.^2)
              * spdiagm(lengthGR)
              * (spdiagm(1./(2 .* lengthGT)))
              * ((2 * AvgX * spdiagm(gradTx) * G1) + (2 * AvgY * spdiagm(gradTy) * G2))
            )

      #dr2Partial = reshape(  r1 .* -1 ./ (lengthGR .* lengthGT.^3)  ,(size(r1)))

      #drc = (    spdiag((r2 .* (AvgX * gradRx) + dr2Partial .* (AvgX * gradTx) )[:]) * AvgX * G1
      #         + spdiag((r2 .* (AvgY * gradRy) + dr2Partial .* (AvgY * gradTy) )[:]) * AvgY * G2 )
      drc = spdiagm(r1)*dr2 + spdiagm(r2) * dr1
      dFunctionValue  = -2*prod(h)*rc'*drc*dT;
      if doHessian
        d2FunctionValue =  2*prod(h) .* dT' * drc' * drc * dT   # note the missing minus sign is not a bug!
      end
    end

    functionValue  = prod(ΩR[2:2:end]-ΩR[1:2:end]) - prod(h) * (rc'*rc)

    #return functionValue[1], dFunctionValue', d2FunctionValue, r1, 1.-(rc.^2)
    return functionValue,dFunctionValue',d2FunctionValue, drc


end
