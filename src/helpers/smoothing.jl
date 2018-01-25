function getGaussianKernel(s::Int, width::Float64)
    σ = width / 3 # full width at tenth maximum = 4.3 * σ
    k = 1./(σ .* sqrt(2*π)).*exp.(-0.5.*(   (0.5:s)   .-s/2).^2/(σ.^2))
    #k = g(0.5:s)
    k = k / sum(k)
    return k
end

function smoothArray(a::Array{Float64,1}, kernelSize::Int, width::Float64;
        kernel::Array{Float64,1}=getGaussianKernel(kernelSize, width))
    a_smooth = zeros(size(a))
    halfKernelSize = Int((kernelSize+1)/2)
    for i = 1:halfKernelSize
        a_smooth[i] = sum(a[1:i])/i
        a_smooth[end-i+1] = sum(a[end-i+1:end])/i
    end
    for i=halfKernelSize:(length(a)-halfKernelSize)
        a_smooth[i] = (a[(i-halfKernelSize+1):(i+halfKernelSize-1)])'*kernel
    end
    return a_smooth
end

function smoothArray(a::Array{Float64,2}, kernelSize::Int, width::Float64)
    kernel = getGaussianKernel(kernelSize, width)
    a_smooth = zeros(size(a))
    for i = 1:size(a,1)
        a_smooth[i,:] = smoothArray(a[i,:], kernelSize, width, kernel=kernel)
    end
    for i = 1:size(a,2)
        a_smooth[:,i] = smoothArray(a_smooth[:,i], kernelSize, width, kernel=kernel)
    end
    return a_smooth
end
