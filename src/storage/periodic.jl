struct PeriodicStorage{T,scalar} <: AbstractStorage{T, scalar}
  data::Array{T,3}
  norb::Int
  N::Int

  function PeriodicStorage(data::AbstractArray{T,3}, scalar=false) where T
    norb = size(data,1)
    N = size(data,3)
    @assert size(data) == (norb, norb, N)
    @assert !(scalar && norb > 1)
    new{T,scalar}(data,norb,N)
  end
end

function PeriodicStorage(::Type{T}, N::Integer, norb=1, scalar=false) where T <: Number
  data = zeros(T, norb, norb, N)
  return PeriodicStorage(data, scalar)
end

PeriodicStorage(N::Integer, norb=1, scalar=false) = PeriodicStorage(ComplexF64, N, norb, scalar)

function Base.getindex(X::PeriodicStorage{T,scalar}, i, j) where {T,scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)

  if scalar
    return i >= j ? X.data[1,1,i-j+1] : X.data[1,1,i-j+X.N]
  else
    return i >= j ? X.data[:,:,i-j+1] : X.data[:,:,i-j+X.N]
  end
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

function Base.:+(X::PeriodicStorage{T,scalar}, Y::PeriodicStorage{T,scalar}) where {T,scalar}
  @assert size(X) == size(Y)
  return PeriodicStorage(X.data .+ Y.data, scalar)
end

function Base.:-(X::PeriodicStorage{T,scalar}, Y::PeriodicStorage{T,scalar}) where {T,scalar}
  @assert size(X) == size(Y)
  return PeriodicStorage(X.data .- Y.data, scalar)
end

function Base.:*(X::PeriodicStorage{T,scalar}, α::Number) where {T,scalar}
  return PeriodicStorage(X.data * α, scalar)
end

Base.:*(α::Number, X::PeriodicStorage) = X * α
