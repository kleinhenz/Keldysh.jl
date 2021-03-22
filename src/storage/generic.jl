struct GenericStorage{T,scalar} <: AbstractStorage{T, scalar}
  data::Array{T,4}
  norb::Int
  N::Int
  M::Int

  function GenericStorage(data::AbstractArray{T,4}, scalar=false) where T
    norb = size(data,1)
    N = size(data,3)
    M = size(data,4)
    @assert size(data) == (norb, norb, N, M)
    @assert !(scalar && norb > 1)
    new{T,scalar}(data,norb,N,M)
  end
end

function GenericStorage(::Type{T}, N::Integer, M::Integer, norb=1, scalar=false) where T <: Number
  data = zeros(T, norb, norb, N, M)
  return GenericStorage(data, scalar)
end

GenericStorage(N::Integer, M::Integer, norb=1, scalar=false) = GenericStorage(ComplexF64, N, M, norb, scalar)

function Base.size(X::GenericStorage)
  return (X.norb, X.norb, X.N, X.M)
end

function Base.getindex(X::GenericStorage{T,scalar}, i, j) where {T, scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.M)
  return scalar ? X.data[1,1,i,j] : X.data[:,:,i,j]
end

function Base.setindex!(X::GenericStorage{T,scalar}, v, i, j) where {T, scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.M)
  if scalar
    return X.data[1,1,i,j] = v
  else
    return X.data[:,:,i,j] = v
  end
end

function Base.setindex!(X::GenericStorage{T,true}, v, i, j) where T
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.M)
  return X.data[1,1,i,j] = v
end

function Base.:+(X::GenericStorage{T,scalar}, Y::GenericStorage{T,scalar}) where {T,scalar}
  @assert size(X) == size(Y)
  return GenericStorage(X.data .+ Y.data, scalar)
end

function Base.:-(X::GenericStorage{T,scalar}, Y::GenericStorage{T,scalar}) where {T,scalar}
  @assert size(X) == size(Y)
  return GenericStorage(X.data .- Y.data, scalar)
end

function Base.:*(X::GenericStorage{T,scalar}, α::Number) where {T,scalar}
  return GenericStorage(X.data * α, scalar)
end

Base.:*(α::Number, X::GenericStorage) = X * α
