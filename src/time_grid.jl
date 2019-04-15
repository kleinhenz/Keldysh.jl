struct TimeGrid
  contour::Contour
  points::Vector{BranchPoint}
  function TimeGrid(c::Contour; npts_real = 0, npts_imag = 0)
    p = BranchPoint[]
    for b in c.branches
      npts = (b.domain == imaginary_branch) ? npts_imag : npts_real
      @assert npts > 1 "every branch must have at least two points (b = $b, npts = $npts)"
      append!(p, get_point.(Ref(b), range(0, 1, length=npts)))
    end
    return new(c, p)
  end
end

