struct AntiHermitianToeplitzStorage{T, scalar} <: AbstractStorage{T, scalar}
  data::Array{T,3}
  norb::Int
  N::Int

  function AntiHermitianToeplitzStorage(data::AbstractArray{T,3}, scalar=false) where T
    norb = size(data,1)
    N = size(data,3)
    @assert size(data) == (norb, norb, N)
    @assert !(scalar && norb > 1)
    new{T,scalar}(data, norb, N)
  end
end

function AntiHermitianToeplitzStorage(::Type{T}, N::Integer, norb=1, scalar=false) where T <: Number
  data = zeros(T, norb, norb, N)
  return AntiHermitianToeplitzStorage(data, scalar)
end

AntiHermitianToeplitzStorage(N::Integer, norb=1, scalar=false) = TimeInvariantAntiHermitian(ComplexF64, N, norb,  scalar)

function Base.getindex(X::AntiHermitianToeplitzStorage{T,scalar}, i, j) where {T, scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)
  i < j && return copy(-adjoint(X[j,i]))
  k = i-j+1
  return scalar ? X.data[1,1,k] : X.data[:,:,k]
end

function Base.setindex!(X::AntiHermitianToeplitzStorage{T,scalar}, v, i, j) where {T, scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)
  if i == j
    @assert iszero(real(v))
  end
  i < j && return X[j,i] = -adjoint(v)
  k = i-j+1
  if scalar
    return X.data[1,1,k] = v
  else
    return X.data[:,:,k] = v
  end
end

function Base.:+(X::AntiHermitianToeplitzStorage{T,scalar}, Y::AntiHermitianToeplitzStorage{T,scalar}) where {T,scalar}
  @assert size(X) == size(Y)
  return AntiHermitianToeplitzStorage(X.data .+ Y.data, scalar)
end

function Base.:-(X::AntiHermitianToeplitzStorage{T,scalar}, Y::AntiHermitianToeplitzStorage{T,scalar}) where {T,scalar}
  @assert size(X) == size(Y)
  return AntiHermitianToeplitzStorage(X.data .- Y.data, scalar)
end

function Base.:*(X::AntiHermitianToeplitzStorage{T,scalar}, α::Number) where {T,scalar}
  return AntiHermitianToeplitzStorage(X.data * α, scalar)
end

Base.:*(α::Number, X::AntiHermitianToeplitzStorage) = X * α
