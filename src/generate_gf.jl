function gf_1level(grid::TimeGrid; epsilon, beta=0.0)
  im_b = get_branch(grid.contour, imaginary_branch)
  if im_b !== nothing
    beta = length(im_b)
  end

  f = (t1, t2) -> -1.0im * (θ(t1, t2) - fermi(epsilon, beta)) * exp(-1.0im * (t1.val - t2.val) * epsilon)
  return f.(grid.points, permutedims(grid.points))
end

