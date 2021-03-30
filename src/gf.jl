abstract type AbstractTimeGF{T} end
#TODO decide on full interface that should be supported by subtypes

@enum GFSignEnum fermionic=-1 bosonic=1

Base.eltype(::Type{<:AbstractTimeGF{T}}) where {T} = T
Base.eltype(X::AbstractTimeGF) = eltype(typeof(X))

function Base.getindex(G::AbstractTimeGF, b1::BranchEnum, b2::BranchEnum)
  grid = G.grid
  @assert b1 ∈ grid.contour && b2 ∈ grid.contour
  x = b1 == backward_branch ? reverse(grid[b1]) : grid[b1]
  y = b2 == backward_branch ? reverse(grid[b2]) : grid[b2]

  norb = norbitals(G)
  n = length(x)
  m = length(y)

  out = zeros(eltype(G), norb, norb, n, m)

  for (i,t1) in enumerate(x)
    for (j,t2) in enumerate(y)
      out[:,:,i,j] .= G[t1, t2]
    end
  end

  return norb == 1 ? out[1,1,:,:] : out
end

function Base.getindex(G::AbstractTimeGF, component::Symbol)
  grid = G.grid
  norb = norbitals(G)

  if component == :greater
    return G[backward_branch, forward_branch]
  elseif component == :lesser
    return G[forward_branch, backward_branch]
  elseif component == :matsubara
    X = real(-im .* G[imaginary_branch, imaginary_branch])
    return norb == 1 ? X[:,1] : X[:,:,:,1]
  elseif component == :retarded
    ret = G[:greater] - G[:lesser]
    return ret
  elseif component == :advanced
    adv = G[:lesser] - G[:greater]
    return adv
  elseif component == :leftmixing
    G[backward_branch, imaginary_branch]
  else
    throw(ArgumentError("component $component not recognized"))
  end
end

include("gf/generic.jl")
include("gf/full.jl")
include("gf/time_invariant_full.jl")
include("gf/keldysh.jl")
include("gf/time_invariant_keldysh.jl")
include("gf/imag.jl")
