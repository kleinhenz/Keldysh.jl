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

function Base.getindex(X::PeriodicStorage{T,scalar}, k, l, i::Int, j::Int, greater=true) where {T,scalar}
  greater = i == j ? greater : i > j
  return greater ? X.data[k,l,i-j+1] : X.data[k,l,i-j+X.N]
end
# need to define these methods to forward greater argument
Base.getindex(X::PeriodicStorage{T,true}, i, j, greater=true) where {T} = X[1,1,i,j, greater]
Base.getindex(X::PeriodicStorage{T,false}, i, j, greater=true) where {T} = X[:,:,i,j, greater]

function Base.setindex!(X::PeriodicStorage{T,scalar}, v, k, l, i::Int, j::Int) where {T,scalar}
  greater = i >= j
  if greater
    return X.data[k,l,i-j+1] = v
  else
    return X.data[k,l,i-j+X.N] = v
  end
end
