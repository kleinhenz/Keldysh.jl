struct TimeGrid
  contour::Contour
  points::Vector{BranchPoint}
  step::Vector{ComplexF64}

  function TimeGrid(c::Contour; npts_real = 0, npts_imag = 0)
    points = BranchPoint[]
    step = ComplexF64[]
    for b in c.branches
      npts = (b.domain == imaginary_branch) ? npts_imag : npts_real
      @assert npts > 1 "every branch must have at least two points (b = $b, npts = $npts)"
      append!(points, get_point.(Ref(b), range(0, 1, length=npts)))
      append!(step, (b.max_val - b.min_val) / (npts - 1))
    end
    return new(c, points, step)
  end
end
