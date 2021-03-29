struct GenericTimeGF{T, scalar, U <: AbstractTimeGrid} <: AbstractTimeGF{T}
  grid::U
  data::GenericStorage{T, scalar}

  function GenericTimeGF(grid::U, data::GenericStorage{T,scalar}) where {T, U <: AbstractTimeGrid, scalar}
    return new{T,scalar, U}(grid, data)
  end
end
norbitals(G::GenericTimeGF) = G.data.norb

function GenericTimeGF(::Type{T}, grid::AbstractTimeGrid, norb=1, scalar=false) where T <: Number
  N = length(grid)
  data = GenericStorage(T, N, N, norb, scalar)
  GenericTimeGF(grid, data)
end
GenericTimeGF(grid::AbstractTimeGrid, norb=1, scalar=false) = GenericTimeGF(ComplexF64, grid, norb, scalar)

function GenericTimeGF(f::Function, ::Type{T}, grid::AbstractTimeGrid, norb=1, scalar=false) where T <: Number
  N = length(grid)
  G = GenericTimeGF(T, grid, norb, scalar)

  for t1 in grid
    for t2 in grid
      G[t1, t2] = f(t1, t2)
    end
  end
  return G
end
GenericTimeGF(f::Function, grid::AbstractTimeGrid, norb=1, scalar=false) = GenericTimeGF(f, ComplexF64, grid, norb, scalar)

@inline function Base.getindex(G::GenericTimeGF, t1::TimeGridPoint, t2::TimeGridPoint, greater=true)
  val = G.data[t1.idx, t2.idx]
  (!greater && t1.idx == t2.idx) && (val += jump(G))
  return val
end

function Base.setindex!(G::GenericTimeGF, v, t1::TimeGridPoint, t2::TimeGridPoint)
  G.data[t1.idx, t2.idx] = v
end

function jump(G::GenericTimeGF)
  t0_plus = branch_bounds(G.grid, forward_branch)[1]
  t0_minus = branch_bounds(G.grid, backward_branch)[2]
  return G[t0_plus, t0_minus] - G[t0_plus, t0_plus]
end

@inline function _line_interp_weights(x, x1, x2)
  lx = x2 - x1
  dx = (x - x1)/lx
  return (1 - dx, dx)
end

@inline function _tri_interp_weights(x, x1, x2, y, y1, y2)
  lx = (x2 - x1)
  dx = (x - x1)/lx

  ly = (y2 - y1)
  dy = (y - y1)/ly

  return (1 - dx, -dy + dx, dy)
end

@inline function _square_interp_weights(x, x1, x2, y, y1, y2)
  lx = (x2 - x1)
  dx = (x - x1)/lx

  ly = (y2 - y1)
  dy = (y - y1)/ly

  w = ((1 - dx) * (1 - dy),
       (1 - dx) * dy,
       dx * (1 - dy),
       dx * dy)

  return w
end

@inline function _line_interp(x, x1, x2, f1, f2)
  w = _line_interp_weights(x, x1, x2)
  return w[1] * f1 + w[2] * f2
end

@inline function _line_interp!(f, x, x1, x2, f1, f2)
  w = _line_interp_weights(x, x1, x2)
  return f .= w[1] .* f1 .+ w[2] .* f2
end

@inline function _tri_interp(x, x1, x2, y, y1, y2, v1, v2, v3)
  w = _tri_interp_weights(x, x1, x2, y, y1, y2)
  return w[1] * v1 + w[2] * v2 + w[3] * v3
end

@inline function _tri_interp!(f, x, x1, x2, y, y1, y2, v1, v2, v3)
  w = _tri_interp_weights(x, x1, x2, y, y1, y2)
  return f .= w[1] .* v1 .+ w[2] .* v2 .+ w[3] .* v3
end

@inline function _square_interp(x, x1, x2, y, y1, y2, f11, f12, f21, f22)
  w = _square_interp_weights(x, x1, x2, y, y1, y2)
  return w[1] * f11 + w[2] * f12 + w[3] * f21 + w[4] * f22
end

@inline function _square_interp!(f, x, x1, x2, y, y1, y2, f11, f12, f21, f22)
  w = _square_interp_weights(x, x1, x2, y, y1, y2)
  return f .= w[1] .* f11 .+ w[2] .* f12 .+ w[3] .* f21 .+ w[4] .* f22
end

function interpolate(G::GenericTimeGF{T, true}, t1::BranchPoint, t2::BranchPoint) where T
  grid = G.grid
  bounds1 = branch_bounds(grid, t1.domain)
  bounds2 = branch_bounds(grid, t2.domain)

  t1_l = find_lower(grid, t1)
  t1_r = grid[t1_l.idx + 1]

  t2_l = find_lower(grid, t2)
  t2_r = grid[t2_l.idx + 1]

  has_cut = t1.domain == t2.domain
  greater = heaviside(t1, t2)

  do_triangular_interp = has_cut && (t1_l.idx == t2_l.idx && t1_r.idx == t2_r.idx)

  if t1 == t2
    return _line_interp(t1.val,
                        t1_l.val.val,
                        t1_r.val.val,
                        G[t1_l, t1_l, greater],
                        G[t1_r, t1_r, greater])
  end

  if !do_triangular_interp
    return _square_interp(t1.val,
                          t1_l.val.val,
                          t1_r.val.val,
                          t2.val,
                          t2_l.val.val,
                          t2_r.val.val,
                          G[t1_l, t2_l, greater],
                          G[t1_l, t2_r, greater],
                          G[t1_r, t2_l, greater],
                          G[t1_r, t2_r, greater])
  else
    if greater # t1 >= t2
      return _tri_interp(t1.val,
                         t1_l.val.val,
                         t1_r.val.val,
                         t2.val,
                         t2_l.val.val,
                         t2_r.val.val,
                         G[t1_l, t2_l, greater],
                         G[t1_r, t2_l, greater],
                         G[t1_r, t2_r, greater])
    else
      return _tri_interp(t2.val,
                         t2_l.val.val,
                         t2_r.val.val,
                         t1.val,
                         t1_l.val.val,
                         t1_r.val.val,
                         G[t1_l, t2_l, greater],
                         G[t1_l, t2_r, greater],
                         G[t1_r, t2_r, greater])
    end
  end
end

function interpolate!(x, G::GenericTimeGF{T, false}, t1::BranchPoint, t2::BranchPoint) where T
  grid = G.grid
  bounds1 = branch_bounds(grid, t1.domain)
  bounds2 = branch_bounds(grid, t2.domain)

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
                        t1_l.val.val,
                        t1_r.val.val,
                        G[t1_l, t1_l, greater],
                        G[t1_r, t1_r, greater])
  end

  if !do_triangular_interp
    return _square_interp!(x,
                          t1.val,
                          t1_l.val.val,
                          t1_r.val.val,
                          t2.val,
                          t2_l.val.val,
                          t2_r.val.val,
                          G[t1_l, t2_l, greater],
                          G[t1_l, t2_r, greater],
                          G[t1_r, t2_l, greater],
                          G[t1_r, t2_r, greater])
  else
    if greater # t1 >= t2
      return _tri_interp!(x,
                         t1.val,
                         t1_l.val.val,
                         t1_r.val.val,
                         t2.val,
                         t2_l.val.val,
                         t2_r.val.val,
                         G[t1_l, t2_l, greater],
                         G[t1_r, t2_l, greater],
                         G[t1_r, t2_r, greater])
    else
      return _tri_interp!(x,
                         t2.val,
                         t2_l.val.val,
                         t2_r.val.val,
                         t1.val,
                         t1_l.val.val,
                         t1_r.val.val,
                         G[t1_l, t2_l, greater],
                         G[t1_l, t2_r, greater],
                         G[t1_r, t2_r, greater])
    end
  end
end

function (G::GenericTimeGF{T, true})(t1::BranchPoint, t2::BranchPoint) where T
  return interpolate(G, t1, t2)
end

function (G::GenericTimeGF{T, false})(t1::BranchPoint, t2::BranchPoint) where T
  norb = norbitals(G)
  x = zeros(T, norb, norb)
  return interpolate!(x, G, t1, t2)
end
