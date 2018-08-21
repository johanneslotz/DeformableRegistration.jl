function checkDerivative(fullF::Function,x;doPlot::Bool=false, smoothOffset::Bool=false)
    df = fullF(x)[2]'
    f(x) = fullF(x)[1]
    return checkDerivative(f,df,x,doPlot=doPlot, smoothOffset=smoothOffset)
end
function checkDerivative(f::Function,df,x;doPlot::Bool=false, smoothOffset::Bool=false)
    N = 17
    h = exp10.(2:-1:-N)
    if smoothOffset
        v = smoothArray(randn(size(x)).-0.5,5,5.0)
    else
        v = randn(size(x)).-0.5
    end
    fx = f(x)
    errquad = zeros(length(h))
    errlin  = zeros(length(h))
    for i=1:length(h)
        errlin[i]  = norm(fx             - f(x+h[i]*v))
        errquad[i] = norm(fx + h[i]*df*v - f(x+h[i]*v))
        s = @sprintf(" h: %1.1e  ||f(x+h*v)||: %2.2e   elin: %2.2e   equad: %2.2e \n",h[i], norm(f(x+h[i]*v)), errlin[i],errquad[i])
        @info s
    end
    #
    if(doPlot)
      PyPlot.rc("legend",fontsize=10)
      PyPlot.rc("xtick", labelsize=10, color="black", direction="in")
      PyPlot.rc("ytick", labelsize=10, color="black", direction="in")
      #fig = PyPlot.figure(figsize=(6,5),facecolor="white")
      PyPlot.clf()
      PyPlot.plot(h,errlin,label="||f(x) - f(x+hv)||")
      PyPlot.plot(h,errquad,".-",label="||f(x) + h nabla f v - f(x+hv)||")
      PyPlot.yscale("log")
      PyPlot.xscale("log")
      PyPlot.xlabel("h")
      PyPlot.ylabel("error")
      PyPlot.legend(loc="upper left")
    end
    return errlin, errquad
end

function checkErrorDecay(errquad::Vector)
  errorReduction = log10.(errquad[1:end-1]./errquad[2:end])
  errorReduction[isinf.(errorReduction)] .= 1
  @debug "# errorReduction.>1.7 = ", sum(errorReduction.>1.7)
  @info "sum(errorReduction...) = ", sum(errorReduction[errorReduction.>1.7]), "   (needs >6)"
  return sum(errorReduction[errorReduction.>1.7])>6
end
