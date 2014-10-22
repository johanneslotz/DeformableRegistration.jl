module Regularizer

export regularizer, createElasticOperatorCentered, createElasticOperatorStaggered, createDiffusiveOperatorCentered

function regularizer(deformationField::Array{Float64,1},
                            operator;
                            doDerivative=false,doHessian=false)

    d2functionValue = operator'*operator
    dfunctionValue  = d2functionValue*deformationField
    functionValue   = (0.5 .* deformationField'*dfunctionValue)[1]
    return functionValue, dfunctionValue, d2functionValue

end


function createElasticOperatorCentered(h::Array{Float64,1},m::Array{Int64,1}; mu::Number = 1, lambda::Number = 0)

  m = m.-1

  a   = sqrt(mu)
  b   = sqrt(mu+lambda)
  dx(k) = spdiagm((-1.*ones(m[k],1),ones(m[k],1)),[0:1],m[k],m[k]+1)/h[k]
  av(k) = spdiagm((ones(m[k],1)/2,ones(m[k],1)/2),[0:1],m[k],m[k]+1)
  D1   = kron(speye(m[2]+1),dx(1))
  D2   = kron(dx(2),speye(m[1]+1))
  A1   = kron( av(2) , speye(m[1]))
  A2   = kron( speye(m[2]),  av(1))
  p1,p2 = size(D1)
  p3,p4 = size(D2)

  B = [ a*D1 spzeros(p1,p2);
       a*D2 spzeros(p3,p4);
       spzeros(p1,p2) a*D1;
       spzeros(p3,p4) a*D2;
       b*A1*D1 b*A2*D2 ]

  return B

end

function createElasticOperatorStaggered(h::Array{Float64,1},m::Array{Int64,1}; mu::Number = 1, lambda::Number = 0)

    a   = sqrt(mu)
    b   = sqrt(mu+lambda)

    dx(h,m) = spdiagm((-ones(m,1),ones(m,1)),[0,1],m,m+1)/h

    D11 = kron(dx(h[2],m[2]),speye(m[1]))
    D21 = kron(speye(m[2]+1),dx(h[1],m[1]-1))
    D12 = kron(dx(h[2],m[2]-1),speye(m[1]+1))
    D22 = kron(speye(m[2]),dx(h[1],m[1]))

    B = spzeros(2*size(D11,1)+size(D21,1)+size(D12,1)+size(D22,1),size(D11,2)+size(D12,2))
    B[1:size(D11,1), 1:size(D11,2)]                         = a*D11
    B[size(D11,1)+1:size(D11,1)+size(D21,1), 1:size(D11,2)] = a*D21
    B[size(D11,1)+size(D21,1)+1:size(D11,1)+size(D21,1)+size(D12,1), size(D11,2)+1:end]                         = a*D12
    B[size(D11,1)+size(D21,1)+size(D12,1)+1:size(D11,1)+size(D21,1)+size(D12,1)+size(D22,1), size(D11,2)+1:end] = a*D22
    B[size(D11,1)+size(D21,1)+size(D12,1)+size(D22,1)+1:end, 1:size(D11,2)]     = b*D11
    B[size(D11,1)+size(D21,1)+size(D12,1)+size(D22,1)+1:end, size(D11,2)+1:end] = b*D22

    return B

end


function createDiffusiveOperatorCentered(h::Array{Float64,1},m::Array{Int64,1})

  d(k) = spdiagm((-ones(m[k]-1,1),ones(m[k]-1,1)),[0 1],m[k]-1,m[k])/(h[k])
  dx = d(1)
  dy = d(2)
  D1 = kron(speye(m[2]),dx)
  D2 = kron(dy,speye(m[1]))
  p1,p2 = size(D1)
  p3,p4 = size(D2)
  B = [ D1 spzeros(p1,p2);
        D2 spzeros(p3,p4);
        spzeros(p1,p2) D1;
        spzeros(p3,p4) D2
      ]
end


end

