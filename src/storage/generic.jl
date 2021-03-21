struct GenericStorage{T,norb} <: AbstractStorage{T, norb}
  data::Array{T,4}
  N::Int
  M::Int

  function GenericStorage(data::AbstractArray{T,4}) where T
    norb = size(data,1)
    N = size(data,3)
    M = size(data,4)
    @assert size(data) == (norb, norb, N, M)
    new{T,norb}(data,N,M)
  end
end

function GenericStorage(::Type{T}, N::Integer, M::Integer, norb=1) where T <: Number
  data = zeros(T, norb, norb, N, M)
  return GenericStorage(data)
end

GenericStorage(N::Integer, M::Integer, norb=1) = GenericStorage(ComplexF64, N, M, norb)

function Base.size(X::GenericStorage)
  norb = norbitals(X)
  return (norb, norb, X.N, X.M)
end


function Base.getindex(X::GenericStorage, i, j)
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.M)

  if norbitals(X) == 1
    return X.data[1,1,i,j]
  else
    return X.data[:,:,i,j]
  end
end

function Base.setindex!(X::GenericStorage, v, i, j)
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.M)

  if norbitals(X) == 1
    return X.data[1,1,i,j] = v
  else
    return X.data[:,:,i,j] = v
  end
end

function Base.:+(X::GenericStorage{T,norb}, Y::GenericStorage{T,norb}) where {T,norb}
  @assert size(X) == size(Y)
  return GenericStorage{T,norb}(X.data .+ Y.data, X.N)
end

function Base.:-(X::GenericStorage{T,norb}, Y::GenericStorage{T,norb}) where {T,norb}
  @assert size(X) == size(Y)
  return GenericStorage{T,norb}(X.data .- Y.data, X.N)
end

function Base.:*(X::GenericStorage{T,norb}, α::Number) where {T,norb}
  return GenericStorage{T,norb}(X.data * α, X.N)
end

Base.:*(α::Number, X::GenericStorage) = X * α
