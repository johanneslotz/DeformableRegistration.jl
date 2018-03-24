h=[1.0,1.0]
m=(1000,100)
import PyPlot


α = 10000



function createDiffusiveOperatorDirichlet1D(h::Array{Float64,1},m::Array{Int64,1})

  d(k) = spdiagm((-ones(m[k],1),ones(m[k],1)),[-1 0],m[k]+1,m[k])/(h[k])
  dx = d(1)
  return dx
end

function createDiffusiveOperatorCentered1D(h::Array{Float64,1},m::Array{Int64,1})

  d(k) = spdiagm((-ones(m[k]-1,1),ones(m[k]-1,1)),[0 1],m[k]-1,m[k])/(h[k])
  dx = d(1)
  return dx
end

function compute1DRegCoarseFine(dfine, dcoarse; nIter=10)
    u = 0*x
    uc = u

    T1(u) = cos.(x*f1 *2*π + u - dcoarse * 2*π)
    T2(u) = amp*cos.(x*f2 *2*π + u - dfine * 2*π) .* mask
    T(u) = (T1(u) + T2(u))
    dT1(u) = -sin.(x*f1 *2*π + u - dcoarse * 2*π)
    dT2(u) = -amp*sin.(x*f2 *2*π + u - dfine * 2*π) .* mask
    dT(u) = dT1(u) + dT2(u)

    R1 = cos.(x*4 *2*π )
    R2 = amp*cos.(x*f2 *2*π) .* mask
    R = R1 + R2
    ∇D = (T(u)-R).*dT(u)

    ## update
    for j=1:nIter
        ∇D = (T(u)-R).*dT(u)
        ∇S = B'*B* u

        H = α * B'*B + spdiagm(dT(u).*dT(u))
        s = H \ (∇D + ∇S)
        u = u - s

        Hc = α * B'*B + spdiagm(dT1(uc).*dT1(uc))
        sc = Hc \ ((T1(uc)-R1).*dT1(uc) + B'*B* uc)
        uc = uc - sc

    end
    # PyPlot.figure(2)
    # PyPlot.plot(x,u-uc, label = "u-uc", color="red", alpha=1/i)
    # PyPlot.axhline(y=0,color="black")
    # PyPlot.axvline(x=x[401], color="black", linestyle=":")
    # PyPlot.axvline(x=x[600], color="black", linestyle=":")
    return uc, u, T1, T2, T, dT1, dT2, dT, R1, R2, R, ∇D
end

B = createDiffusiveOperatorCentered1D(h,[m[1],m[2]])
#B = createDiffusiveOperatorDirichlet1D(h,[m[1],m[2]])
#B = B'*B

x = linspace(0,1,m[1])
mask = zeros(size(x))
mask[401:600]=1

f1 = 4
f2 = 40
amp = 1
dcoarse = 1/4

## deformation error vs. true deformation

emax = zeros(80)
for i=1:length(emax)
    dfine = dcoarse + 1/f2*(i-1)
    uc, u, T1, T2, T, dT1, dT2, dT, R1, R2, R, ∇D  = compute1DRegCoarseFine(dfine, dcoarse, nIter=100)
    emax[i] = maximum(abs.(u-uc))
end

PyPlot.figure(1)
PyPlot.clf()
wavelength = 1/f2
PyPlot.plot(2*π  *1/f2*(0:length(emax)-1) , emax)
PyPlot.xlabel("additional deformation in high frequency image patch")
PyPlot.ylabel("maximum error in low frequency solution")
PyPlot.tight_layout()
## Model Problem with two sine waves

dcoarse = 1/4
dfine = dcoarse + 1/20

uc, u, T1, T2, T, dT1, dT2, dT, R1, R2, R, ∇D = compute1DRegCoarseFine(dfine, dcoarse, nIter=100)

PyPlot.figure(2)
PyPlot.clf()

PyPlot.subplot(411)
PyPlot.title("data")
#PyPlot.plot(x,T1(0), label = "T1")
#PyPlot.plot(x,T2(0), label = "T2")
PyPlot.plot(x,T1(0)+T2(0), label = "T1 + T2")
PyPlot.plot(x,R, label = "R1 + R2")

PyPlot.axvline(x=x[401], color="black", linestyle=":")
PyPlot.axvline(x=x[600], color="black", linestyle=":")
PyPlot.legend()

PyPlot.subplot(412)
PyPlot.title("coarse + fine")
PyPlot.plot(x,T(u), label = "T")
PyPlot.plot(x,R, label = "R")
PyPlot.plot(x,∇D, label = "∇D")
PyPlot.plot(x,u, label = "u")
PyPlot.plot([x[401] x[600]][:],2*pi*[dfine dfine][:], color="black", linestyle="--")
PyPlot.axhline(y= dcoarse*2π, color="black", linestyle="--", label="true deformation")
PyPlot.axvline(x=x[401], color="black", linestyle=":")
PyPlot.axvline(x=x[600], color="black", linestyle=":")
PyPlot.legend()

PyPlot.subplot(413)
PyPlot.title("coarse")

PyPlot.plot(x,T1(uc), label = "T1")
PyPlot.plot(x,R1, label = "R1")
PyPlot.plot(x,(T1(uc)-R1).*dT1(uc), label = "∇D1")
PyPlot.plot(x,uc, label = "u")
PyPlot.axhline(y= dcoarse*2π, color="black", linestyle="--")
PyPlot.axvline(x=x[401], color="black", linestyle=":")
PyPlot.axvline(x=x[600], color="black", linestyle=":")

PyPlot.legend()


PyPlot.subplot(414)
PyPlot.title("difference")
PyPlot.plot(x,(u-uc)/(2*π), label = "u-uc", color="red")
PyPlot.axhline(y=0,color="black")
PyPlot.axvline(x=x[401], color="black", linestyle=":")
PyPlot.axvline(x=x[600], color="black", linestyle=":")

PyPlot.legend()
PyPlot.xlabel("deformation difference inside: $(abs(dcoarse-dfine))")
PyPlot.tight_layout()
