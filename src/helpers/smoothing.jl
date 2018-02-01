function getGaussianKernel(s::Int, width::Float64)
    σ = width / 4.3 # full width at tenth maximum = 4.3 * σ
    k = 1./(σ .* sqrt(2*π)).*exp.(-0.5.*(   (0.5:s)   .-s/2).^2/(σ.^2))
    #k = g(0.5:s)
    k = k ./ sum(k)
    return k
end

function smoothArray(a::Array{Float64,1}, kernelSize::Int, width::Float64;
        kernel::Array{Float64,1}=getGaussianKernel(kernelSize, width))
    a_smooth = zeros(size(a))
    halfKernelSize = Int((kernelSize-1)/2)
    for i = 1:halfKernelSize
        a_smooth[i] = a[i]
        a_smooth[end-i+1] = a[end-i+1]
    end
    for i=halfKernelSize+1:(length(a)-halfKernelSize)
        a_smooth[i] = (a[(i-halfKernelSize):(i+halfKernelSize)])'*kernel
    end
    return a_smooth
end

function smoothArray(a::Array{Float64,2}, kernelSize::Int, width::Float64)
    kernel = getGaussianKernel(kernelSize, width)
    #kernel=[1/6, 2/3, 1/6]
    a_smooth = zeros(size(a))
    for i = 1:size(a,1)
        a_smooth[i,:] = smoothArray(a[i,:], kernelSize, width, kernel=kernel)
    end
    for i = 1:size(a,2)
        a_smooth[:,i] = smoothArray(a_smooth[:,i], kernelSize, width, kernel=kernel)
    end
    return a_smooth
end

function smoothArray(a::scaledArray, kernelSize::Int, width::Float64)
    x = reshape(a.data[1:Int(end/2)], a.dimensions)
    xs = smoothArray(x, kernelSize, width)[:]
    y = reshape(a.data[1+Int(end/2):end], a.dimensions)
    ys = smoothArray(y, kernelSize, width)[:]

    return scaledArray(vcat(xs,ys), a.dimensions, a.voxelsize, a.shift)
end
