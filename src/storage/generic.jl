struct GenericStorage{T,scalar} <: AbstractStorage{T, scalar}
  data::Array{T,4}
  norb::Int
  N::Int
  M::Int

  function GenericStorage{T, scalar}(data::AbstractArray{T,4}) where {T, scalar}
    norb = size(data,1)
    N = size(data,3)
    M = size(data,4)
    @assert size(data) == (norb, norb, N, M)
    @assert !(scalar && norb > 1)
    return new(data, norb, N, M)
  end
end

function GenericStorage(data::AbstractArray{T,4}, scalar=false) where T
  GenericStorage{T,scalar}(data)
end

function GenericStorage(::Type{T}, N::Integer, M::Integer, norb=1, scalar=false) where T <: Number
  data = zeros(T, norb, norb, N, M)
  return GenericStorage(data, scalar)
end

GenericStorage(N::Integer, M::Integer, norb=1, scalar=false) = GenericStorage(ComplexF64, N, M, norb, scalar)

@inline function Base.getindex(X::GenericStorage{T,scalar}, k, l, i, j) where {T, scalar}
  return X.data[k,l,i,j]
end

@inline function Base.setindex!(X::GenericStorage{T,scalar}, v, k, l, i, j) where {T, scalar}
  return X.data[k,l,i,j] = v
end
