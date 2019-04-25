struct TimeGridPoint
  idx::Int64
  val::BranchPoint
end

struct TimeGrid
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
      branch_points = get_point.(Ref(b), range(0, 1, length=npts))
      time_grid_points = TimeGridPoint.(indices, branch_points)

      append!(points, time_grid_points)
      append!(step, (b.max_val - b.min_val) / (npts - 1))
      push!(branch_bounds, (time_grid_points[1], time_grid_points[end]))
    end
    return new(c, points, step, branch_bounds)
  end
end

function integrate(f, grid::TimeGrid)
  integrate(f, grid, grid.points[end], grid.points[1])
end

function integrate(f, grid::TimeGrid, t1::TimeGridPoint, t2::TimeGridPoint)
  # wrap around
  if (t1.idx < t2.idx)
    return integrate(f, grid, t1, grid.points[1]) + integrate(f, grid, grid.points[end], t2)
  end

  T = typeof(fieldtype(BranchPoint, :val)(0.0) * f(t2))
  s = zero(T)

  first_branch_idx = grid.contour.branch_indices[Int(t2.val.domain)]
  last_branch_idx =  grid.contour.branch_indices[Int(t1.val.domain)]

  first = t2
  for b in first_branch_idx:last_branch_idx
    Δt = grid.step[b]
    last = (b == last_branch_idx ? t1 : grid.branch_bounds[b][2])
    m = last.idx - first.idx + 1
    @assert first.val.domain == last.val.domain
    if (m >= 2) #trapezoid rule
      sb = zero(T)
      sb += 0.5 * (f(grid.points[first.idx]) + f(grid.points[last.idx]))
      for t in grid.points[(first.idx+1):(last.idx-1)]
        sb += f(t)
      end
      s += (Δt * sb)
    end
    b == last_branch_idx || (first = grid.points[last.idx + 1])
  end
  return s
end
