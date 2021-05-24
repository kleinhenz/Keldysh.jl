struct AntiHermitianStorage{T,scalar} <: AbstractStorage{T, scalar}
  data::Array{T,3}
  norb::Int
  N::Int

  function AntiHermitianStorage{T, scalar}(data::AbstractArray{T,3}) where {T, scalar}
    norb = size(data,1)
    M = size(data,3)
    @assert norb == size(data,1) == size(data,2)
    @assert !(scalar && norb > 1)
    N = div(-1 + isqrt(1 + 8M), 2)
    @assert div(N*(N+1),2) == M
    new(data,norb,N)
  end
end

function AntiHermitianStorage(data::AbstractArray{T,3}, scalar=false) where T
  AntiHermitianStorage{T,scalar}(data)
end

function AntiHermitianStorage(data::AbstractArray{T,4}, scalar=false) where T
  norb = size(data,1)
  N = size(data,3)
  @assert size(data) == (norb, norb, N, N)
  M = div(N*(N+1), 2)

  lo = zeros(T, norb, norb, M)
  k = 1
  for i in 1:N
    for j in 1:i
      lo[:,:,k] .= data[:,:,i,j]
      k += 1
    end
  end

  AntiHermitianStorage(lo, scalar)
end

function AntiHermitianStorage(::Type{T}, N::Integer, norb=1, scalar=false) where T <: Number
  data = zeros(T, norb, norb, div(N*(N+1), 2))
  return AntiHermitianStorage(data, scalar)
end

AntiHermitianStorage(N::Integer, norb=1, scalar=false) = AntiHermitianStorage(ComplexF64, N, norb, scalar)

function Base.getindex(X::AntiHermitianStorage{T,scalar}, i, j) where {T, scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)
  i < j && return copy(-adjoint(X[j,i]))
  n = div((i-1)*i,2)
  n += j
  return scalar ? X.data[1,1,n] : X.data[:,:,n]
end

function Base.getindex(X::AntiHermitianStorage{T,scalar}, k, l, i, j) where {T, scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N) && (1 <= k <= X.norb) && (1 <= l <= X.norb)
  i < j && return -conj(X[l,k,j,i])
  n = div((i-1)*i,2)
  n += j
  return X.data[k,l,n]
end

function Base.setindex!(X::AntiHermitianStorage{T,scalar}, v, i, j) where {T, scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)
  i < j && return X[j,i] = -adjoint(v)
  n = div((i-1)*i,2)
  n += j
  if scalar
    return X.data[1,1,n] = v
  else
    return X.data[:,:,n] = v
  end
end

function Base.setindex!(X::AntiHermitianStorage{T,scalar}, v, k, l, i, j) where {T, scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N) && (1 <= k <= X.norb) && (1 <= l <= X.norb)
  i < j && return X[l,k,j,i] = -conj(v)
  n = div((i-1)*i,2)
  n += j
  X.data[k,l,n] = v
end

function Base.size(X::AntiHermitianStorage)
  return (X.norb, X.norb, X.N, X.N)
end

function Base.similar(X::T) where T <: AntiHermitianStorage
  T(similar(X.data))
end

function Base.zero(X::T) where T <: AntiHermitianStorage
  T(zero(X.data))
end

function Base.:+(X::AntiHermitianStorage{T,scalar}, Y::AntiHermitianStorage{T,scalar}) where {T,scalar}
  @assert size(X) == size(Y)
  return AntiHermitianStorage(X.data .+ Y.data, scalar)
end

function Base.:-(X::AntiHermitianStorage{T,scalar}, Y::AntiHermitianStorage{T,scalar}) where {T,scalar}
  @assert size(X) == size(Y)
  return AntiHermitianStorage(X.data .- Y.data, scalar)
end

function Base.:*(X::AntiHermitianStorage{T,scalar}, α::Number) where {T,scalar}
  return AntiHermitianStorage(X.data * α, scalar)
end

Base.:*(α::Number, X::AntiHermitianStorage) = X * α
