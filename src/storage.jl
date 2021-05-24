#TODO should this be <: AbstractArray{T,4}?
"""
    AbstractStorage{T,scalar}

Abstract supertype for storage of different Keldysh components of nonequilibrium Green's functions taking into account symmetries.
Represents an array with dimensions (norb,norb,N,M).
Supports basic arithmetic and indexing.
"""
abstract type AbstractStorage{T, scalar} end

Base.eltype(::Type{<:AbstractStorage{T,scalar}}) where {T, scalar} = T
Base.eltype(X::AbstractStorage) = eltype(typeof(X))

function Base.size(X::AbstractStorage)
  return (X.norb, X.norb, X.N, X.M)
end

function Base.similar(X::T) where T <: AbstractStorage
  T(similar(X.data))
end

function Base.zero(X::T) where T <: AbstractStorage
  T(zero(X.data))
end

function Base.:+(X::T, Y::T) where T <: AbstractStorage
  @assert size(X) == size(Y)
  return T(X.data .+ Y.data)
end

function Base.:-(X::T, Y::T) where T <: AbstractStorage
  @assert size(X) == size(Y)
  return T(X.data .- Y.data)
end

function Base.:*(X::T, α::Number) where T <: AbstractStorage
  return T(X.data * α)
end

function Base.:*(α::Number, X::T) where T <: AbstractStorage
  return X * α
end

include("storage/generic.jl")
include("storage/herm.jl")
include("storage/herm_toeplitz.jl")
include("storage/periodic.jl")
