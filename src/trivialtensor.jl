@inline space(A::AbstractArray) = Tuple(TrivialVectorSpace(n) for n in size(A))
@inline space(A::AbstractArray, l::Int) = TrivialVectorSpace(size(A, l))
