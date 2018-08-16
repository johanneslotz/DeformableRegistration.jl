using Stats
using Printf
# using Logging

function estimateNGFEpsilon(referenceImage::regImage;cutoffPercent = 80)
  m = size(referenceImage.data)
  h = referenceImage.voxelsize
  Rc = referenceImage.data


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

  gradTx  = G1*referenceImage.data[:]
  gradTy  = G2*referenceImage.data[:]

  quantileX = percentile(AvgX * (gradTx.^2),cutoffPercent)
  quantileY = percentile(AvgY * (gradTy.^2),cutoffPercent)

  estimatedEps = sqrt(0.5*(quantileX + quantileY))

  lengthGT = sqrt(sum(AvgX * (gradTx.^2) + AvgY * (gradTy.^2), 2) + estimatedEps^2) #epsilon-norm

  vectorlength = reshape(lengthGT,(size(referenceImage)))

  try
    figure(51)
    clf()
    showImage(vectorlength)
    s = @sprintf("estimated epsilon: %3.3e", estimatedEps)
    title(s)
    colorbar()
    figure(52)
    clf()
    plt.hist(lengthGT,bins=100)
  catch err
    warn("Error while plotting. This is normal if it happens in the automated test framework.")
    println(err)
  end

  return (estimatedEps,vectorlength)


end
