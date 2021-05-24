struct AntiHermitianStorage{T,scalar} <: AbstractStorage{T, scalar}
  data::Array{T,3}
  norb::Int
  N::Int
  M::Int

  function AntiHermitianStorage{T, scalar}(data::AbstractArray{T,3}) where {T, scalar}
    norb = size(data,1)
    m = size(data,3)
    @assert norb == size(data,1) == size(data,2)
    @assert !(scalar && norb > 1)
    N = div(-1 + isqrt(1 + 8m), 2)
    @assert div(N*(N+1),2) == m
    new(data,norb,N,N)
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

function Base.getindex(X::AntiHermitianStorage{T,scalar}, k, l, i::Int, j::Int) where {T, scalar}
  i < j && return copy(-adjoint(X[l,k,j,i])) #FIXME k,l must BOTH be integer or colons for this to make sense
  n = div((i-1)*i,2)
  n += j
  return X.data[k,l,n]
end

function Base.setindex!(X::AntiHermitianStorage{T,scalar}, v, k, l, i::Int, j::Int) where {T, scalar}
  i < j && return X[l,k,j,i] = -adjoint(v)
  n = div((i-1)*i,2)
  n += j
  X.data[k,l,n] = v
end
