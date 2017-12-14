# helpers to make addition/multiplication of distance and regularizer possible

import Base.+
function +(a::Function, b::SparseMatrixCSC)
	f(x) = a(x) + b*x
	return f
end
function +(b::Function, a::Function)
	f(x) = a(x) + b(x)
	return f
end
