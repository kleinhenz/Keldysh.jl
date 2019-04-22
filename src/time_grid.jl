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
