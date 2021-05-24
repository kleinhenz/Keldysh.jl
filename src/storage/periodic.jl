"""
    PeriodicStorage{T,scalar} <: AbstractStorage{T, scalar}

Storage for a (discontinuous) shift invariant, periodic function ``G(τ1, τ2) = G(τ1 - τ2)``, ``G(τ) = G(τ + β)``. \\
Represents data on the uniform grid ``{0, Δτ, 2Δτ, ..., (N-1)Δτ}`` where ``Δτ = β/(N-1)``. \\
Note don't generically have ``G(0) = G(β)`` so we need to store both endpoints. \\
Because of this, there is a discontinuity along the diagonal and the storage is not quite circulant. \\
Diagonal discontinuity is handled by `greater` argument to `getindex` method.
"""
struct PeriodicStorage{T,scalar} <: AbstractStorage{T, scalar}
  data::Array{T,3}
  norb::Int
  N::Int
  M::Int

  function PeriodicStorage{T,scalar}(data::AbstractArray{T,3}) where {T, scalar}
    norb = size(data,1)
    N = size(data,3)
    @assert size(data) == (norb, norb, N)
    @assert !(scalar && norb > 1)
    new(data,norb,N,N)
  end
end

function PeriodicStorage(data::AbstractArray{T,3}, scalar=false) where T
  PeriodicStorage{T,scalar}(data)
end

function PeriodicStorage(::Type{T}, N::Integer, norb=1, scalar=false) where T <: Number
  data = zeros(T, norb, norb, N)
  return PeriodicStorage(data, scalar)
end

PeriodicStorage(N::Integer, norb=1, scalar=false) = PeriodicStorage(ComplexF64, N, norb, scalar)

function Base.getindex(X::PeriodicStorage{T,scalar}, i, j, greater=true) where {T,scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)

  greater = i == j ? greater : i > j

  if scalar
    return greater ? X.data[1,1,i-j+1] : X.data[1,1,i-j+X.N]
  else
    return greater ? X.data[:,:,i-j+1] : X.data[:,:,i-j+X.N]
  end
end

function Base.getindex(X::PeriodicStorage{T,scalar}, k, l, i, j, greater=true) where {T,scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N) && (1 <= k <= X.norb) && (1 <= l <= X.norb)

  greater = i == j ? greater : i > j
  return greater ? X.data[k,l,i-j+1] : X.data[k,l,i-j+X.N]
end

function Base.setindex!(X::PeriodicStorage{T,scalar}, v, i, j) where {T,scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)

  if scalar
    if i >= j
      return X.data[1,1,i-j+1] = v
    else
      return X.data[1,1,i-j+X.N] = v
    end
  else
    if i >= j
      return X.data[:,:,i-j+1] = v
    else
      return X.data[:,:,i-j+X.N] = v
    end
  end
end

function Base.setindex!(X::PeriodicStorage{T,scalar}, v, k, l, i, j) where {T,scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N) && (1 <= k <= X.norb) && (1 <= l <= X.norb)
  greater = i >= j
  if greater
    return X.data[k,l,i-j+1] = v
  else
    return X.data[k,l,i-j+X.N] = v
  end
end
