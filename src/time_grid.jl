struct TimeGridPoint
  idx::Int64 # contour index
  ridx::Int64 # real/imag time index
  bpoint::BranchPoint
end

abstract type AbstractTimeGrid <: AbstractVector{TimeGridPoint} end

function _make_grid_helper(branches, npts)
  points = TimeGridPoint[]
  branch_bounds_ = Pair{TimeGridPoint, TimeGridPoint}[]

  for (b, n) in zip(branches, npts)
    @assert n > 1 "every branch must have at least two points (b = $b, npts = $n)"

    indices = (1:n) .+ length(points)

    r_indices = b.domain == backward_branch ? reverse(1:n) : 1:n

    branch_points = b.(range(0, 1, length=n))
    time_grid_points = TimeGridPoint.(indices, r_indices, branch_points)

    append!(points, time_grid_points)
    push!(branch_bounds_, time_grid_points[1]=>time_grid_points[end])
  end

  branch_bounds = ntuple(i -> branch_bounds_[i], length(branches))
  return points, branch_bounds
end

struct FullTimeGrid <: AbstractTimeGrid
  contour::FullContour
  points::Vector{TimeGridPoint}
  branch_bounds::NTuple{3, Pair{TimeGridPoint, TimeGridPoint}}
  nt::Int
  ntau::Int

  function FullTimeGrid(c::FullContour, nt, ntau)
    npts = ntuple(i -> c.branches[i].domain == imaginary_branch ? ntau : nt, 3)
    points, branch_bounds = _make_grid_helper(c.branches, npts)
    return new(c, points, branch_bounds, nt, ntau)
  end
end

struct KeldyshTimeGrid <: AbstractTimeGrid
  contour::KeldyshContour
  points::Vector{TimeGridPoint}
  branch_bounds::NTuple{2, Pair{TimeGridPoint, TimeGridPoint}}
  nt::Int

  function KeldyshTimeGrid(c::KeldyshContour, nt)
    npts = (nt, nt)
    points, branch_bounds = _make_grid_helper(c.branches, npts)
    return new(c, points, branch_bounds, nt)
  end
end

struct ImaginaryTimeGrid <: AbstractTimeGrid
  contour::ImaginaryContour
  points::Vector{TimeGridPoint}
  branch_bounds::NTuple{1, Pair{TimeGridPoint, TimeGridPoint}}
  ntau::Int

  function ImaginaryTimeGrid(c::ImaginaryContour, ntau)
    npts = (ntau,)
    points, branch_bounds = _make_grid_helper(c.branches, npts)
    return new(c, points, branch_bounds, ntau)
  end
end

### AbstractArray Interface ###
IndexStyle(::Type{<:AbstractTimeGrid}) = IndexLinear()
Base.size(grid::AbstractTimeGrid) = size(grid.points)
Base.getindex(grid::AbstractTimeGrid, i::Int) = grid.points[i]
Base.setindex!(grid::AbstractTimeGrid, v::TimeGridPoint, i::Int) = grid.points[i] = v

function Base.length(grid::AbstractTimeGrid, b::BranchEnum)
  bounds = branch_bounds(grid, b)
  return bounds[2].idx - bounds[1].idx + 1
end

function branch_bounds(grid::AbstractTimeGrid, b::BranchEnum)
  @assert b ∈ grid.contour
  idx = findfirst(x -> x.domain == b, grid.contour.branches)
  return grid.branch_bounds[idx]
end

function Base.getindex(grid::AbstractTimeGrid, b::BranchEnum)
  bounds = branch_bounds(grid, b)
  return view(grid, bounds[1].idx:bounds[2].idx)
end

function Base.step(grid::AbstractTimeGrid, b::BranchEnum)::ComplexF64
  @assert b ∈ grid.contour
  if b == forward_branch
    return grid.contour.tmax / (grid.nt - 1)
  elseif b == backward_branch
    return -grid.contour.tmax / (grid.nt - 1)
  else
    return -im * grid.contour.β / (grid.ntau - 1)
  end
end

function integrate(f, grid::AbstractTimeGrid, t1::TimeGridPoint, t2::TimeGridPoint, init=0.0im)
  # wrap around
  (t1.idx < t2.idx) && return integrate(f, grid, t1, grid[1]) + integrate(f, grid, grid[end], t2)

  integral = init

  first_branch_idx = findfirst(x -> x.domain == t2.bpoint.domain, grid.contour.branches)
  last_branch_idx = findfirst(x -> x.domain == t1.bpoint.domain, grid.contour.branches)

  first = t2
  for b in first_branch_idx:last_branch_idx
    Δt = step(grid, grid.contour.branches[b].domain)
    last = (b == last_branch_idx ? t1 : grid.branch_bounds[b][2])
    @assert first.bpoint.domain == last.bpoint.domain
    if (last.idx != first.idx) #trapezoid rule
      branch_integral = 0.5 * (f(grid[first.idx]) + f(grid[last.idx]))
      for i in (first.idx+1):(last.idx-1)
        t = @inbounds grid[i]
        branch_integral += f(t)
      end
      integral += (Δt * branch_integral)
    end
    b == last_branch_idx || (first = grid[last.idx + 1])
  end
  return integral
end

integrate(f, grid::AbstractTimeGrid, init=0.0im) = integrate(f, grid, grid[end], grid[1], init)

realtimes(grid::FullTimeGrid) = range(0.0, grid.contour.tmax, length=grid.nt)
realtimes(grid::KeldyshTimeGrid) = range(0.0, grid.contour.tmax, length=grid.nt)

imagtimes(grid::FullTimeGrid) = range(0.0, grid.contour.β, length=grid.ntau)
imagtimes(grid::ImaginaryTimeGrid) = range(0.0, grid.contour.β, length=grid.ntau)

"""
    struct TimeDomain <: Any

Container for a vector of pairs of `TimeGridPoint`s.
Represents the domain over which a Green's function is defined.

Fields
======

`points :: Vector{Tuple{TimeGridPoint, TimeGridPoint}}`
"""
struct TimeDomain
  points::Vector{Tuple{TimeGridPoint, TimeGridPoint}}
end

function find_lower(grid::AbstractTimeGrid, t::BranchPoint)
  bounds = branch_bounds(grid, t.domain)
  if t.ref == 1.0 # don't hit the bound
    idx = bounds.second.idx - 1
  else
    npts = t.domain == imaginary_branch ? grid.ntau : grid.nt
    step = 1 / (npts - 1)
    idx = bounds.first.idx + Int(floor(t.ref / step))
  end
  return grid[idx]
end
