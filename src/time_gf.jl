using LinearAlgebra

abstract type AbstractTimeGF end

struct TimeGF <: AbstractTimeGF
  data::Array{ComplexF64, 4}
  grid::TimeGrid
end

function TimeGF(grid::TimeGrid, norb = 1)
  N = length(grid)
  data = zeros(ComplexF64, norb, norb, N, N)
  TimeGF(data, grid)
end

function TimeGF(f::Function, grid::TimeGrid, norb = 1; lower = false)
  N = length(grid)
  gf = TimeGF(grid, norb)

  for t1 in grid
    for t2 in grid
      lower && t1.idx < t2.idx && continue
      gf[t1, t2] = f(t1, t2)
    end
  end
  return gf
end

function TimeGF(les::AbstractArray{ComplexF64,4},
                ret::AbstractArray{ComplexF64,4},
                tv::AbstractArray{ComplexF64,4},
                mat::AbstractArray{Float64,3},
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
    if ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) &&
        (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
      i = findfirst(isapprox(real(t1.val.val)), t)
      j = findfirst(isapprox(real(t2.val.val)), t)
      greater ? ret[:,:,i,j] .+ les[:,:,i,j] : les[:,:,i,j]
    elseif (t1.val.domain == imaginary_branch && (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
      i = findfirst(isapprox(imag(t1.val.val)), -tau)
      j = findfirst(isapprox(real(t2.val.val)), t)
      conj(tv[:,:,j, ntau+1-i]) # akoi 19c
    elseif ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) && t2.val.domain == imaginary_branch)
      i = findfirst(isapprox(real(t1.val.val)), t)
      j = findfirst(isapprox(imag(t2.val.val)), -tau)
      tv[:,:,i,j]
    else
      i = findfirst(isapprox(imag(t1.val.val)), -tau)
      j = findfirst(isapprox(imag(t2.val.val)), -tau)
      greater ? 1.0im * mat[:,:,i - j + 1] : -1.0im * mat[:,:,i - j + ntau]
    end
  end

  return G
end

function TimeGF(les::AbstractArray{ComplexF64,4},
                ret::AbstractArray{ComplexF64,4},
                grid::TimeGrid)

  @assert grid.contour.domain == keldysh_contour

  t = realtimes(grid)

  nt = length(t)

  @assert size(les) == (nt,nt)
  @assert size(ret) == (nt,nt)

  G = TimeGF(grid) do t1, t2
    greater = heaviside(t1.val, t2.val)
    i = findfirst(isapprox(real(t1.val.val)), t)
    j = findfirst(isapprox(real(t2.val.val)), t)
    greater ? ret[:,:,i,j] .+ les[:,:,i,j] : les[:,:,i,j]
  end

  return G
end

#### Indexing ###
Base.@propagate_inbounds function Base.getindex(gf::TimeGF, i::Int, j::Int)
  gf.data[:, :, i, j]
end

Base.@propagate_inbounds function Base.setindex!(gf::TimeGF, v::AbstractMatrix, i::Int, j::Int)
  gf.data[:, :, i, j] = v
end

# indexing with TimeGridPoint
Base.@propagate_inbounds function getindex(gf::TimeGF, t1::TimeGridPoint, t2::TimeGridPoint, gtr=true)
  val = gf[t1.idx, t2.idx]
  (!gtr && t1.idx == t2.idx) && (val += jump(gf))
  return val
end
Base.@propagate_inbounds function setindex!(gf::TimeGF, v::AbstractMatrix, t1::TimeGridPoint, t2::TimeGridPoint)
  gf[t1.idx, t2.idx] = v
end

function jump(gf::TimeGF)
  t0_plus = branch_bounds(gf.grid, forward_branch)[1]
  t0_minus = branch_bounds(gf.grid, backward_branch)[2]
  return gf[t0_plus, t0_minus] - gf[t0_plus, t0_plus]
end

function getindex(G::TimeGF, b1::BranchEnum, b2::BranchEnum)
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

  return out
end

function getindex(G::TimeGF, component::Symbol)
  grid = G.grid
  if component == :greater
    return G[backward_branch, forward_branch]
  elseif component == :lesser
    return G[forward_branch, backward_branch]
  elseif component == :matsubara
    real(-im .* G[imaginary_branch, imaginary_branch][:,:,:,1])
  elseif component == :retarded
    ret = G[:greater] - G[:lesser]
    return ret
  elseif component == :advanced
    adv = G[:lesser] - G[:greater]
    return adv
  elseif component == :leftmixing
    gf[backward_branch, imaginary_branch]
  else
    throw(ArgumentError("component $component not recognized"))
  end
end
