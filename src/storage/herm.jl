struct AntiHermitianStorage{T,norb} <: AbstractStorage{T, norb}
  data::Array{T,3}
  N::Int

  function AntiHermitianStorage(data::AbstractArray{T,3}) where T
    norb = size(data,1)
    M = size(data,3)
    @assert norb == size(data,1) == size(data,2)

    N = div(-1 + isqrt(1 + 8M), 2)
    @assert div(N*(N+1),2) == M
    new{T,norb}(data,N)
  end
end

function AntiHermitianStorage(data::AbstractArray{T,4}) where T
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

  AntiHermitianStorage(lo)
end

function AntiHermitianStorage(::Type{T}, N::Integer, norb=1) where T <: Number
  data = zeros(T, norb, norb, div(N*(N+1), 2))
  return AntiHermitianStorage(data)
end

AntiHermitianStorage(N::Integer, norb=1) = AntiHermitianStorage(ComplexF64, N, norb)

function Base.getindex(X::AntiHermitianStorage, i, j)
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)

  i < j && return -conj(X[j,i])

  k = div((i-1)*i,2)
  k += j

  if norbitals(X) == 1
    return X.data[1,1,k]
  else
    return X.data[:,:,k]
  end
end

function Base.setindex!(X::AntiHermitianStorage, v, i, j)
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)

  if i == j
    @assert iszero(real(v))
  end

  i < j && return X[j,i] = -conj(v)

  k = div((i-1)*i,2)
  k += j

  if norbitals(X) == 1
    return X.data[1,1,k] = v
  else
    return X.data[:,:,k] = v
  end
end

function Base.:+(X::AntiHermitianStorage{T,norb}, Y::AntiHermitianStorage{T,norb}) where {T,norb}
  @assert size(X) == size(Y)
  return AntiHermitianStorage{T,norb}(X.data .+ Y.data, X.N)
end

function Base.:-(X::AntiHermitianStorage{T,norb}, Y::AntiHermitianStorage{T,norb}) where {T,norb}
  @assert size(X) == size(Y)
  return AntiHermitianStorage{T,norb}(X.data .- Y.data, X.N)
end

function Base.:*(X::AntiHermitianStorage{T,norb}, α::Number) where {T,norb}
  return AntiHermitianStorage{T,norb}(X.data * α, X.N)
end

Base.:*(α::Number, X::AntiHermitianStorage) = X * α
