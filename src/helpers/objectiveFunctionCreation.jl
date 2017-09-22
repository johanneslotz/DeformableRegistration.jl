# helpers to make addition/multiplication of distance and regularizer possible

function *(factor::Number,Regularizer::Tuple{Number,Vector,SparseMatrixCSC{Float64,Int64}})
    return factor.*Regularizer[1], factor.*Regularizer[2], factor.*Regularizer[3]
end

function +(Distance::Tuple{Number,Number,Number,Tuple{Vector,Vector}},Regularizer::Tuple{Number,Vector,SparseMatrixCSC{Float64,Int64}})
    return Distance[1].+Regularizer[1], 0, 0
end

function +(Distance::Tuple{Number,Number,Number},Regularizer::Tuple{Number,Vector,SparseMatrixCSC{Float64,Int64}})
    return Distance[1].+Regularizer[1], 0, 0
end

function +(Distance::Tuple{Number,Vector,Function,Tuple{Vector,Vector}},Regularizer::Tuple{Number,Vector,SparseMatrixCSC{Float64,Int64}})
    d2J(grid) = Distance[3](grid) + Regularizer[3]*grid
    return Distance[1].+Regularizer[1], Distance[2].+Regularizer[2], d2J
end

function +(Distance::Tuple{Number,Vector,Function},Regularizer::Tuple{Number,Vector,SparseMatrixCSC{Float64,Int64}})
    d2J(grid) = Distance[3](grid) + Regularizer[3]*grid
    return Distance[1].+Regularizer[1], Distance[2].+Regularizer[2], d2J
end

function +(Distance::Tuple{Number,Vector,SparseMatrixCSC{Float64,Int64},Tuple{Vector,Vector}},Regularizer::Tuple{Number,Vector,SparseMatrixCSC{Float64,Int64}})
    return Distance[1].+Regularizer[1], Distance[2].+Regularizer[2], Distance[3].+Regularizer[3]
end
