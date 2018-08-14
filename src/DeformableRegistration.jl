module DeformableRegistration



module Types
    using Images
    struct scaledArray
               data::Array{Float64, 1}
               dimensions::Tuple{Vararg{Int64}}
               voxelsize::Array{Float64, 1}
               shift::Array{Float64, 1}
    end

    struct regImage
               data::ImageMeta
               voxelsize::Array{Float64, 1}
               shift::Array{Float64, 1}
    end



    import Base.+
    import Base.-
    import Base.*
    import Base.size

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

    export scaledArray, regImage, + , - , *

end

# helping functions
include("helpers/matrixFreeHelpers.jl")
include("helpers/registrationOptions.jl");
#include("helpers/objectiveFunctionCreation.jl")

# submodules
include("ImageProcessing.jl")
include("Transformation.jl")
include("Interpolation.jl")

include("Distance.jl")
include("Regularizer.jl")
include("Optimization.jl")

include("Examples.jl")

try
    include("Visualization.jl")
catch err
    println("caught: ", err)
    println("This is 'normal' during automated testing: ")
end

end
