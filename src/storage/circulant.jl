struct CirculantStorage{T,scalar} <: AbstractStorage{T, scalar}
  data::Array{T,3}
  norb::Int
  N::Int

  function CirculantStorage(data::AbstractArray{T,3}, scalar=false) where T
    norb = size(data,1)
    N = size(data,3)
    @assert size(data) == (norb, norb, N)
    @assert !(scalar && norb > 1)
    new{T,scalar}(data,norb,N)
  end
end

function CirculantStorage(::Type{T}, N::Integer, norb=1, scalar=false) where T <: Number
  data = zeros(T, norb, norb, N)
  return CirculantStorage(data, scalar)
end

CirculantStorage(N::Integer, norb=1, scalar=false) = CirculantStorage(ComplexF64, N, norb, scalar)

function Base.getindex(X::CirculantStorage{T,scalar}, i, j, greater=true) where {T,scalar}
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)

  greater = i == j ? greater : i > j

  if scalar
    return greater ? X.data[1,1,i-j+1] : X.data[1,1,i-j+X.N]
  else
    return greater ? X.data[:,:,i-j+1] : X.data[:,:,i-j+X.N]
  end
end

function Base.setindex!(X::CirculantStorage{T,scalar}, v, i, j) where {T,scalar}
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

function Base.:+(X::CirculantStorage{T,scalar}, Y::CirculantStorage{T,scalar}) where {T,scalar}
  @assert size(X) == size(Y)
  return CirculantStorage(X.data .+ Y.data, scalar)
end

function Base.:-(X::CirculantStorage{T,scalar}, Y::CirculantStorage{T,scalar}) where {T,scalar}
  @assert size(X) == size(Y)
  return CirculantStorage(X.data .- Y.data, scalar)
end

function Base.:*(X::CirculantStorage{T,scalar}, α::Number) where {T,scalar}
  return CirculantStorage(X.data * α, scalar)
end

Base.:*(α::Number, X::CirculantStorage) = X * α
