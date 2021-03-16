using LinearAlgebra

abstract type AbstractTimeGF{T, norb} end

norbitals(::Type{<:AbstractTimeGF{T, norb}}) where {T,norb} = norb
norbitals(x::AbstractTimeGF) = norbitals(typeof(x))

function Base.getindex(G::AbstractTimeGF, b1::BranchEnum, b2::BranchEnum)
  grid = G.grid
  @assert b1 ∈ grid.contour && b2 ∈ grid.contour
  x = b1 == backward_branch ? reverse(grid[b1]) : grid[b1]
  y = b2 == backward_branch ? reverse(grid[b2]) : grid[b2]

  n = length(x)
  m = length(y)
  norb = size(G.data,1)

  out = zeros(eltype(G.data), norb, norb, n, m)

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

struct TimeGF{T, norb} <: AbstractTimeGF{T, norb}
  grid::TimeGrid
  data::Array{T, 4}

  function TimeGF(grid::TimeGrid, data::Array{T,4}) where T
    @assert size(data,1) == size(data,2)
    norb = size(data,1)
    return new{T,norb}(grid, data)
  end
end

function TimeGF(::Type{T}, grid::TimeGrid, norb = 1) where T <: Number
  N = length(grid)
  data = zeros(T, norb, norb, N, N)
  TimeGF(grid, data)
end
TimeGF(grid::TimeGrid, norb = 1) = TimeGF(ComplexF64, grid, norb)

function TimeGF(f::Function, ::Type{T}, grid::TimeGrid, norb = 1; lower = false) where T <: Number
  N = length(grid)
  G = TimeGF(T, grid, norb)

  for t1 in grid
    for t2 in grid
      lower && t1.idx < t2.idx && continue
      G[t1, t2] = f(t1, t2)
    end
  end
  return G
end
TimeGF(f::Function, grid::TimeGrid, norb = 1; lower = false) = TimeGF(f, ComplexF64, grid, norb; lower)

AbstractArray4 = AbstractArray{T,4} where T
AbstractArray3 = AbstractArray{T,3} where T

function TimeGF(les::AbstractArray4,
                ret::AbstractArray4,
                tv::AbstractArray4,
                mat::AbstractArray3,
                grid::TimeGrid)

  @assert grid.contour.domain == full_contour

  t = realtimes(grid)
  tau = imagtimes(grid)

  nt = length(t)
  ntau = length(tau)
  norb = size(les, 1)

  @assert size(les) == (norb, norb, nt, nt)
  @assert size(ret) == (norb, norb, nt, nt)
  @assert size(tv) == (norb, norb, nt, ntau)
  @assert size(mat) == (norb, norb, ntau)

  G = TimeGF(grid, norb) do t1, t2
    greater = heaviside(t1.val, t2.val)
    i = t1.ridx
    j = t2.ridx

    x =
    if ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) &&
        (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
      greater ? ret[:,:,i,j] .+ les[:,:,i,j] : les[:,:,i,j]
    elseif (t1.val.domain == imaginary_branch && (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
      conj(tv[:,:,j, ntau+1-i]) # akoi 19c
    elseif ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) && t2.val.domain == imaginary_branch)
      tv[:,:,i,j]
    else
      greater ? 1.0im * mat[:,:,i - j + 1] : -1.0im * mat[:,:,i - j + ntau]
    end

    return norb == 1 ? x[] : x
  end

  return G
end

function TimeGF(les::AbstractArray4,
                ret::AbstractArray4,
                grid::TimeGrid)

  @assert grid.contour.domain == keldysh_contour

  t = realtimes(grid)

  nt = length(t)
  norb = size(les, 1)

  @assert size(les) == (norb,norb,nt,nt)
  @assert size(ret) == (norb,norb,nt,nt)

  G = TimeGF(grid,norb) do t1, t2
    greater = heaviside(t1.val, t2.val)
    i = t1.ridx
    j = t2.ridx
    x = greater ? ret[:,:,i,j] .+ les[:,:,i,j] : les[:,:,i,j]
    return norb == 1 ? x[] : x
  end

  return G
end

#### Indexing ###
Base.@propagate_inbounds function Base.getindex(G::TimeGF, i::Int, j::Int)
  if norbitals(G) == 1
    return G.data[1,1,i,j]
  else
    return G.data[:,:,i,j]
  end
end

Base.@propagate_inbounds function Base.setindex!(G::TimeGF, v, i::Int, j::Int)
  if norbitals(G) == 1
    G.data[1, 1, i, j] = v
  else
    G.data[:,:,i,j] = v
  end
end

# indexing with TimeGridPoint
function Base.getindex(G::TimeGF, t1::TimeGridPoint, t2::TimeGridPoint, gtr=true)
  val = @inbounds G[t1.idx, t2.idx]
  (!gtr && t1.idx == t2.idx) && (val += jump(G))
  return val
end

function Base.setindex!(G::TimeGF, v, t1::TimeGridPoint, t2::TimeGridPoint)
  @inbounds G[t1.idx, t2.idx] = v
end

function jump(G::TimeGF)
  t0_plus = branch_bounds(G.grid, forward_branch)[1]
  t0_minus = branch_bounds(G.grid, backward_branch)[2]
  return G[t0_plus, t0_minus] - G[t0_plus, t0_plus]
end
