using Logging

# try
#   using PyPlot
# catch e
#   println("caught: ", e)
# end

function checkDerivative(fullF::Function,x;doPlot::Bool=false)
    df = fullF(x)[2]'
    f(x) = fullF(x)[1]
    return checkDerivative(f,df,x,doPlot=doPlot)
end
function checkDerivative(f,df,x;doPlot::Bool=false)
    N = 17
    h = exp10.(2:-1:-N)
    v = randn(size(x)).-0.5
    fx = f(x)
    errquad = zeros(length(h))
    errlin  = zeros(length(h))
    for i=1:length(h)
        errlin[i]  = norm(fx             - f(x+h[i]*v))
        errquad[i] = norm(fx + h[i]*df*v - f(x+h[i]*v))
        s = @sprintf(" h: %5e  ||f(x+h*v)||: %5e   elin: %5e   equad: %5e \n",h[i], norm(f(x+h[i]*v)), errlin[i],errquad[i])
        Logging.info(s)
    end
    #
    # if(doPlot)
    #   PyPlot.rc("legend",fontsize=10)
    #   PyPlot.rc("xtick", labelsize=10, color="black", direction="in")
    #   PyPlot.rc("ytick", labelsize=10, color="black", direction="in")
    #   fig = figure(figsize=(6,5),facecolor="white")
    #   plot(h,errlin,label=L"$||f(x) - f(x+hv)||$")
    #   plot(h,errquad,".-",label=L"$||f(x) + h \nabla f v - f(x+hv)||$")
    #   yscale("log")
    #   xscale("log")
    #   xlabel(L"$h$")
    #   ylabel(L"$error$")
    #   legend(loc="upper left")
    # end
    return errlin, errquad
end

function checkErrorDecay(errquad::Vector)
  errorReduction = log10.(errquad[1:end-1]./errquad[2:end])
  errorReduction[isinf.(errorReduction)]=1
  return sum(errorReduction[errorReduction.>1.7])>6
end
