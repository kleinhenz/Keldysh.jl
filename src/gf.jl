# TODO decide on full interface that should be supported by subtypes
# An AbstractTimeGF G should support at least:
# * G[k,l,t1,t2] where k,l are orbital indices and t1,t2 are TimeGridPoints
# * TimeDomain(G)
# * similar(G)
# * zero(G)
# * norbitals(G)
# * G.grid
"""
    AbstractTimeGF{T, scalar}

Abstract supertype for Keldysh Green's function objects.\\
An `AbstractTimeGF` object `G` supports indexing with `TimeGridPoint`s and interpolation with `BranchPoint`s.
"""
abstract type AbstractTimeGF{T, scalar} end

@enum GFSignEnum fermionic=-1 bosonic=1

Base.eltype(::Type{<:AbstractTimeGF{T, scalar}}) where {T, scalar} = T
Base.eltype(X::AbstractTimeGF) = eltype(typeof(X))

is_scalar(::Type{<:AbstractTimeGF{T, scalar}}) where {T, scalar} = scalar
is_scalar(X::AbstractTimeGF) = is_scalar(typeof(X))

Base.getindex(X::AbstractTimeGF{T,true}, i, j, greater=true) where {T} = X[1,1,i,j,greater]
Base.getindex(X::AbstractTimeGF{T,false}, i, j, greater=true) where {T} = X[:,:,i,j,greater]
Base.setindex!(X::AbstractTimeGF{T,true}, v, i, j) where {T} = X[1,1,i,j] = v
Base.setindex!(X::AbstractTimeGF{T,false}, v, i, j) where {T} = X[:,:,i,j] = v

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

function interpolate(G::AbstractTimeGF{T, true}, t1::BranchPoint, t2::BranchPoint) where T
  grid = G.grid

  t1_l = find_lower(grid, t1)
  t1_r = grid[t1_l.idx + 1]

  t2_l = find_lower(grid, t2)
  t2_r = grid[t2_l.idx + 1]

  has_cut = t1.domain == t2.domain
  greater = heaviside(t1, t2)

  do_triangular_interp = has_cut && (t1_l.idx == t2_l.idx && t1_r.idx == t2_r.idx)

  if t1 == t2
    return _line_interp(t1.val,
                        t1_l.bpoint.val,
                        t1_r.bpoint.val,
                        G[t1_l, t1_l, greater],
                        G[t1_r, t1_r, greater])
  end

  if !do_triangular_interp
    return _square_interp(t1.val,
                          t1_l.bpoint.val,
                          t1_r.bpoint.val,
                          t2.val,
                          t2_l.bpoint.val,
                          t2_r.bpoint.val,
                          G[t1_l, t2_l, greater],
                          G[t1_l, t2_r, greater],
                          G[t1_r, t2_l, greater],
                          G[t1_r, t2_r, greater])
  else
    if greater # t1 >= t2
      return _tri_interp(t1.val,
                         t1_l.bpoint.val,
                         t1_r.bpoint.val,
                         t2.val,
                         t2_l.bpoint.val,
                         t2_r.bpoint.val,
                         G[t1_l, t2_l, greater],
                         G[t1_r, t2_l, greater],
                         G[t1_r, t2_r, greater])
    else
      return _tri_interp(t2.val,
                         t2_l.bpoint.val,
                         t2_r.bpoint.val,
                         t1.val,
                         t1_l.bpoint.val,
                         t1_r.bpoint.val,
                         G[t1_l, t2_l, greater],
                         G[t1_l, t2_r, greater],
                         G[t1_r, t2_r, greater])
    end
  end
end

function interpolate!(x, G::AbstractTimeGF{T, false}, t1::BranchPoint, t2::BranchPoint) where T
  grid = G.grid

  t1_l = find_lower(grid, t1)
  t1_r = grid[t1_l.idx + 1]

  t2_l = find_lower(grid, t2)
  t2_r = grid[t2_l.idx + 1]

  has_cut = t1.domain == t2.domain
  greater = heaviside(t1, t2)

  do_triangular_interp = has_cut && (t1_l.idx == t2_l.idx && t1_r.idx == t2_r.idx)

  if t1 == t2
    return _line_interp!(x,
                        t1.val,
                        t1_l.bpoint.val,
                        t1_r.bpoint.val,
                        G[t1_l, t1_l, greater],
                        G[t1_r, t1_r, greater])
  end

  if !do_triangular_interp
    return _square_interp!(x,
                          t1.val,
                          t1_l.bpoint.val,
                          t1_r.bpoint.val,
                          t2.val,
                          t2_l.bpoint.val,
                          t2_r.bpoint.val,
                          G[t1_l, t2_l, greater],
                          G[t1_l, t2_r, greater],
                          G[t1_r, t2_l, greater],
                          G[t1_r, t2_r, greater])
  else
    if greater # t1 >= t2
      return _tri_interp!(x,
                         t1.val,
                         t1_l.bpoint.val,
                         t1_r.bpoint.val,
                         t2.val,
                         t2_l.bpoint.val,
                         t2_r.bpoint.val,
                         G[t1_l, t2_l, greater],
                         G[t1_r, t2_l, greater],
                         G[t1_r, t2_r, greater])
    else
      return _tri_interp!(x,
                         t2.val,
                         t2_l.bpoint.val,
                         t2_r.bpoint.val,
                         t1.val,
                         t1_l.bpoint.val,
                         t1_r.bpoint.val,
                         G[t1_l, t2_l, greater],
                         G[t1_l, t2_r, greater],
                         G[t1_r, t2_r, greater])
    end
  end
end

function (G::AbstractTimeGF{T, true})(t1::BranchPoint, t2::BranchPoint) where T
  return interpolate(G, t1, t2)
end

function (G::AbstractTimeGF{T, false})(t1::BranchPoint, t2::BranchPoint) where T
  norb = norbitals(G)
  x = zeros(T, norb, norb)
  return interpolate!(x, G, t1, t2)
end

include("gf/generic.jl")
include("gf/full.jl")
include("gf/time_invariant_full.jl")
include("gf/keldysh.jl")
include("gf/time_invariant_keldysh.jl")
include("gf/imag.jl")
