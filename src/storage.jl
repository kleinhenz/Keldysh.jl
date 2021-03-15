struct CompactAntiHermitian{T,norb}
  data::Array{T,3}
  N::Int
end

function CompactAntiHermitian(::Type{T}, N, norb=1) where T <: Number
  data = zeros(norb, norb, div(N*(N+1), 2))
  return CompactAntiHermitian{T, norb}(data, N)
end

CompactAntiHermitian(N, norb=1) = CompactAntiHermitian(ComplexF64, N, norb)

norbitals(::Type{<:CompactAntiHermitian{T,norb}}) where {T, norb} = norb
norbitals(X::CompactAntiHermitian) = norbitals(typeof(X))

Base.eltype(::Type{<:CompactAntiHermitian{T,norb}}) where {T, norb} = T
Base.eltype(X::CompactAntiHermitian) = eltype(typeof(X))

function Base.size(X::CompactAntiHermitian)
  norb = norbitals(X)
  return (norb, norb, X.N, X.N)
end

function Base.getindex(X::CompactAntiHermitian, i, j)
  @boundscheck begin
    @assert 1 <= i <= X.N
    @assert 1 <= j <= X.N
  end

  i < j && return -conj(X[j,i])

  k = div((i-1)*i,2)
  k += j

  if norbitals(X) == 1
    return X.data[1,1,k]
  else
    return X.data[:,:,k]
  end
end

function Base.setindex!(X::CompactAntiHermitian, v, i, j)
  @boundscheck begin
    @assert 1 <= i <= X.N
    @assert 1 <= j <= X.N
  end

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

function Base.:+(X::CompactAntiHermitian{T,norb}, Y::CompactAntiHermitian{T,norb}) where {T,norb}
  @assert X.N == Y.N
  return CompactAntiHermitian{T,norb}(X.data .+ Y.data, X.N)
end

function Base.:-(X::CompactAntiHermitian{T,norb}, Y::CompactAntiHermitian{T,norb}) where {T,norb}
  @assert X.N == Y.N
  return CompactAntiHermitian{T,norb}(X.data .- Y.data, X.N)
end

function Base.:*(X::CompactAntiHermitian{T,norb}, α::Number) where {T,norb}
  return CompactAntiHermitian{T,norb}(X.data * α, X.N)
end

Base.:*(α::Number, X::CompactAntiHermitian) = X * α
