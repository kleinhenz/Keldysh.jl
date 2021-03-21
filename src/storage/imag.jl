struct ImaginaryTimeStorage{T,norb} <: AbstractStorage{T, norb}
  data::Array{T,3}
  N::Int

  function ImaginaryTimeStorage(data::AbstractArray{T,3}) where T
    norb = size(data,1)
    N = size(data,3)
    @assert size(data) == (norb, norb, N)
    new{T,norb}(data,N)
  end
end

function ImaginaryTimeStorage(::Type{T}, N::Integer, norb=1) where T <: Number
  data = zeros(T, norb, norb, N)
  return ImaginaryTimeStorage(data)
end

ImaginaryTimeStorage(N::Integer, norb=1) = ImaginaryTimeStorage(ComplexF64, N, norb)

function Base.getindex(X::ImaginaryTimeStorage, i, j)
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)

  if norbitals(X) == 1
    return i >= j ? X.data[1,1,i-j+1] : -X.data[1,1,i-j+X.N]
  else
    return i >= j ? X.data[:,:,i-j+1] : -X.data[:,:,i-j+X.N]
  end
end

function Base.setindex!(X::ImaginaryTimeStorage, v, i, j)
  @boundscheck @assert (1 <= i <= X.N) && (1 <= j <= X.N)

  if norbitals(X) == 1
    if i >= j
      return X.data[1,1,i-j+1] = v
    else
      return X.data[1,1,i-j+X.N] = -v
    end
  else
    if i >= j
      return X.data[:,:,i-j+1] = v
    else
      return X.data[:,:,i-j+X.N] = -v
    end
  end
end

function Base.:+(X::ImaginaryTimeStorage{T,norb}, Y::ImaginaryTimeStorage{T,norb}) where {T,norb}
  @assert size(X) == size(Y)
  return ImaginaryTimeStorage{T,norb}(X.data .+ Y.data, X.N)
end

function Base.:-(X::ImaginaryTimeStorage{T,norb}, Y::ImaginaryTimeStorage{T,norb}) where {T,norb}
  @assert size(X) == size(Y)
  return ImaginaryTimeStorage{T,norb}(X.data .- Y.data, X.N)
end

function Base.:*(X::ImaginaryTimeStorage{T,norb}, α::Number) where {T,norb}
  return ImaginaryTimeStorage{T,norb}(X.data * α, X.N)
end

Base.:*(α::Number, X::ImaginaryTimeStorage) = X * α
