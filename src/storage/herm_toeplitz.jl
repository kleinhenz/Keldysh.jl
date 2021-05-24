struct AntiHermitianToeplitzStorage{T, scalar} <: AbstractStorage{T, scalar}
  data::Array{T,3}
  norb::Int
  N::Int
  M::Int

  function AntiHermitianToeplitzStorage{T,scalar}(data::AbstractArray{T,3}) where {T, scalar}
    norb = size(data,1)
    N = size(data,3)
    @assert size(data) == (norb, norb, N)
    @assert !(scalar && norb > 1)
    new(data, norb, N, N)
  end
end

function AntiHermitianToeplitzStorage(data::AbstractArray{T,3}, scalar=false) where T
  AntiHermitianToeplitzStorage{T,scalar}(data)
end

function AntiHermitianToeplitzStorage(data::AbstractArray{T,4}, scalar=false) where T
  norb = size(data,1)
  N = size(data,3)
  @assert size(data) == (norb, norb, N, N)
  col = zeros(T, norb, norb, N)
  for i in 1:N
    col[:,:,i] = data[:,:,i,1]
  end
  AntiHermitianToeplitzStorage(col, scalar)
end

function AntiHermitianToeplitzStorage(::Type{T}, N::Integer, norb=1, scalar=false) where T <: Number
  data = zeros(T, norb, norb, N)
  return AntiHermitianToeplitzStorage(data, scalar)
end

AntiHermitianToeplitzStorage(N::Integer, norb=1, scalar=false) = TimeInvariantAntiHermitian(ComplexF64, N, norb,  scalar)

function Base.getindex(X::AntiHermitianToeplitzStorage{T,scalar}, i, j) where {T, scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)
  i < j && return copy(-adjoint(X[j,i]))
  n = i-j+1
  return scalar ? X.data[1,1,n] : X.data[:,:,n]
end

function Base.getindex(X::AntiHermitianToeplitzStorage{T,scalar}, k, l, i, j) where {T, scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N) && (1 <= k <= X.norb) && (1 <= l <= X.norb)
  i < j && return -conj(X[l,k,j,i])
  n = i-j+1
  return X.data[k,l,n]
end

function Base.setindex!(X::AntiHermitianToeplitzStorage{T,scalar}, v, i, j) where {T, scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)
  i < j && return X[j,i] = -adjoint(v)
  n = i-j+1
  if scalar
    return X.data[1,1,n] = v
  else
    return X.data[:,:,n] = v
  end
end

function Base.setindex!(X::AntiHermitianToeplitzStorage{T,scalar}, v, k, l, i, j) where {T, scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N) && (1 <= k <= X.norb) && (1 <= l <= X.norb)
  i < j && return X[l,k,j,i] = -conj(v)
  n = i-j+1
  X.data[k,l,n] = v
end
