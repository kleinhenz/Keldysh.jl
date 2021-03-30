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

