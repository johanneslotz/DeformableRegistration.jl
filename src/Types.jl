module Types
    import Base.+
    import Base.-
    import Base.*
    import Base.size
    export scaledArray, regImage, + , - , *, t2a

    struct scaledArray
               data::Array{Float64, 1}
               dimensions::Tuple{Vararg{Int64}}
               voxelsize::Array{Float64, 1}
               shift::Array{Float64, 1}
    end

    struct regImage
               data::Array{Float64,2}
               voxelsize::Array{Float64, 1}
               shift::Array{Float64, 1}
    end


    function +(a::scaledArray, b::scaledArray)
        @assert a.voxelsize == b.voxelsize "Array's world matrices must match"
        @assert a.shift == b.shift "Array's world matrices must match"
        return(scaledArray(a.data+b.data, a.dimensions, a.voxelsize, a.shift))
    end

    function +(a::scaledArray, b::Array{Float64,1}) # for derivative check
        return(scaledArray(a.data+b, a.dimensions, a.voxelsize, a.shift))
    end

    function +(b::Number, a::scaledArray ) # for derivative check
        return(scaledArray(a.data+b, a.dimensions, a.voxelsize, a.shift))
    end

    function -(a::scaledArray, b::scaledArray)
        return(a + (-1 * b))
    end

    function *(s::Number, a::scaledArray)
        b = deepcopy(a)
        b.data[:] = s * b.data[:]
        return b
    end

    function size(a::scaledArray)
        return size(a.data)
    end

    function size(a::regImage)
        return size(a.data)
    end



    # helpers to make addition/multiplication of distance and regularizer possible
    function +(a::Function, b::SparseMatrixCSC)
    	f(x) = a(x) + b*x
    	return f
    end
    function +(b::Function, a::Function)
    	f(x) = a(x) + b(x)
    	return f
    end
    function *(a::Number, b::Function)
    	f(x) = a * b(x)
    	return f
    end


    function t2a(t::Tuple)
    	return [a for a in t]
    end


end
