using ImageRegistration.ImageProcessing
using ImageRegistration.Transformation
using ImageRegistration.Distance
using Images
using Stats


function estimateNGFEpsilon(referenceImage::Image)

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

    gradTx  = G1*referenceImage[:]
    gradTy  = G2*referenceImage[:]

    quantileX = percentile(AvgX * (gradTx.^2),95)
    quantileY = percentile(AvgY * (gradTy.^2),95)

    estimatedEps = sqrt(0.5*(quantileX + quantileY))

    return estimatedEps


end

