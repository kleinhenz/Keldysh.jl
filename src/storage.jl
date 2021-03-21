abstract type AbstractStorage{T, norb} end

norbitals(::Type{<:AbstractStorage{T,norb}}) where {T, norb} = norb
norbitals(X::AbstractStorage) = norbitals(typeof(X))

Base.eltype(::Type{<:AbstractStorage{T,norb}}) where {T, norb} = T
Base.eltype(X::AbstractStorage) = eltype(typeof(X))

include("storage/generic.jl")
include("storage/herm.jl")
include("storage/imag.jl")
