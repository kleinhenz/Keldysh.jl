struct AntiHermitianStorage{T,scalar} <: AbstractStorage{T, scalar}
  data::Array{T,3}
  norb::Int
  N::Int

  function AntiHermitianStorage(data::AbstractArray{T,3}, scalar=false) where T
    norb = size(data,1)
    M = size(data,3)
    @assert norb == size(data,1) == size(data,2)
    @assert !(scalar && norb > 1)

    N = div(-1 + isqrt(1 + 8M), 2)
    @assert div(N*(N+1),2) == M
    new{T,scalar}(data,norb,N)
  end
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
  i < j && return -conj(X[j,i])
  k = div((i-1)*i,2)
  k += j
  return scalar ? X.data[1,1,k] : X.data[:,:,k]
end

function Base.setindex!(X::AntiHermitianStorage{T,scalar}, v, i, j) where {T, scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)
#  if i == j
#    @assert iszero(real(v))
#  end
  i < j && return X[j,i] = -conj(v)
  k = div((i-1)*i,2)
  k += j
  if scalar
    return X.data[1,1,k] = v
  else
    return X.data[:,:,k] = v
  end
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
