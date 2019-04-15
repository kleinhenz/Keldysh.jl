import Base.length

@enum BranchEnum forward_branch backward_branch imaginary_branch

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
  return BranchPoint((b.max_val - b.min_val) * ref, ref, b.domain)
end

function length(b::Branch)
  return abs(b.max_val - b.min_val)
end
