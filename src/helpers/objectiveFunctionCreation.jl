# helpers to make addition/multiplication of distance and regularizer possible

import Base.+
import Base.*
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
