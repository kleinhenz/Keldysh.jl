abstract type AbstractStorage{T, scalar} end

Base.eltype(::Type{<:AbstractStorage{T,scalar}}) where {T, scalar} = T
Base.eltype(X::AbstractStorage) = eltype(typeof(X))

include("storage/generic.jl")
include("storage/herm.jl")
include("storage/time_invariant_herm.jl")
include("storage/periodic.jl")
