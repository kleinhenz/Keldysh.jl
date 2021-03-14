struct TimeGridPoint
  idx::Int64 # contour index
  ridx::Int64 # real/imag time index
  val::BranchPoint
end

struct TimeGrid <: AbstractVector{TimeGridPoint}
  contour::Contour
  points::Vector{TimeGridPoint}
  step::Vector{ComplexF64}
  branch_bounds::Vector{Tuple{TimeGridPoint, TimeGridPoint}}

  function TimeGrid(c::Contour; npts_real = 0, npts_imag = 0)
    points = TimeGridPoint[]
    step = ComplexF64[]
    branch_bounds = Tuple{TimeGridPoint, TimeGridPoint}[]

    for b in c.branches
      npts = (b.domain == imaginary_branch) ? npts_imag : npts_real
      @assert npts > 1 "every branch must have at least two points (b = $b, npts = $npts)"

      indices = (1:npts) .+ length(points)

      r_indices = b.domain == backward_branch ? reverse(1:npts) : 1:npts

      branch_points = b.(range(0, 1, length=npts))
      time_grid_points = TimeGridPoint.(indices, r_indices, branch_points)

      append!(points, time_grid_points)
      append!(step, (b.max_val - b.min_val) / (npts - 1))
      push!(branch_bounds, (time_grid_points[1], time_grid_points[end]))
    end
    return new(c, points, step, branch_bounds)
  end
end


### AbstractArray Interface ###
IndexStyle(::Type{<:TimeGrid}) = IndexLinear()
Base.size(grid::TimeGrid) = size(grid.points)
Base.getindex(grid::TimeGrid, i::Int) = grid.points[i]
Base.setindex!(grid::TimeGrid, v::TimeGridPoint, i::Int) = grid.points[i] = v

"get ith point on branch b of grid"
function Base.getindex(grid::TimeGrid, b::BranchEnum, i::Int)
  bounds = branch_bounds(grid, b)
  return grid[bounds[1].idx + i - 1]
end

"get all points on branch b of grid"
function Base.getindex(grid::TimeGrid, b::BranchEnum)
  bounds = branch_bounds(grid, b)
  return grid[bounds[1].idx:bounds[2].idx]
end

function branch_bounds(grid::TimeGrid, b::BranchEnum)
  return grid.branch_bounds[grid.contour.branch_indices[Int(b)]]
end

function Base.length(grid::TimeGrid, b::BranchEnum)
  bounds = branch_bounds(grid, b)
  return bounds[2].idx - bounds[1].idx + 1
end

function Base.step(grid::TimeGrid, b::BranchEnum)
  return grid.step[grid.contour.branch_indices[Int(b)]]
end

function integrate(f, grid::TimeGrid, t1::TimeGridPoint, t2::TimeGridPoint)
  # wrap around
  (t1.idx < t2.idx) && return integrate(f, grid, t1, grid[1]) + integrate(f, grid, grid[end], t2)

  integral = zero(fieldtype(BranchPoint, :val)(0.0) * f(t2)) # extra evaluation just to get the type

  first_branch_idx = grid.contour.branch_indices[Int(t2.val.domain)]
  last_branch_idx =  grid.contour.branch_indices[Int(t1.val.domain)]

  first = t2
  for b in first_branch_idx:last_branch_idx
    Δt = grid.step[b]
    last = (b == last_branch_idx ? t1 : grid.branch_bounds[b][2])
    @assert first.val.domain == last.val.domain
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

integrate(f, grid::TimeGrid) = integrate(f, grid, grid[end], grid[1])

realtimes(grid::TimeGrid) = map(t -> real(t.val.val), grid[forward_branch])
imagtimes(grid::TimeGrid) = map(t -> real(1.0im * t.val.val), grid[imaginary_branch]) # τ = it
