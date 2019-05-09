@enum BranchEnum forward_branch=1 backward_branch=2 imaginary_branch=3

struct BranchPoint
  val::ComplexF64
  ref::Float64
  domain::BranchEnum
end

struct Branch
  domain::BranchEnum
  min_val::ComplexF64
  max_val::ComplexF64
  function Branch(d::BranchEnum, len::Real)
    @assert len > 0.0 "branch must have non-zero length (d = $d, len = $len)"
    if d == forward_branch
      return new(d, 0.0, len)
    elseif d == backward_branch
      return new(d, len, 0.0)
    else
      return new(d, 0.0, -1.0im * len)
    end
  end
end

"""
maps from the unit interval [0, 1] to a point on a branch
"""
function get_point(b::Branch, ref::Real)
  @assert 0.0 <= ref <= 1.0
  return BranchPoint(ref * b.max_val + (1 - ref) * b.min_val, ref, b.domain)
end

(b::Branch)(ref::Real) = get_point(b, ref)

function length(b::Branch)
  return abs(b.max_val - b.min_val)
end
