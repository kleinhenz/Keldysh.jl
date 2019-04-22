@enum ContourEnum full_contour=1 keldysh_contour=2 imaginary_contour=3

function branch_set(c::ContourEnum)
  if c == full_contour
    return Set([forward_branch, backward_branch, imaginary_branch])
  elseif c == keldysh_contour
    return Set([forward_branch, backward_branch])
  else
    return Set([imaginary_branch])
  end
end

function get_contour_enum(branches::Set{BranchEnum})
  for c in instances(ContourEnum)
    branch_set(c) == branches && return c
  end
  @assert false "invalid branch set"
end

function get_contour_enum(branches::AbstractVector{Branch})
  branch_set = Set([b.domain for b in branches])
  @assert length(branch_set) == length(branches) "branches must be unique"
  return get_contour_enum(branch_set)
end

struct Contour
  domain::ContourEnum
  branches::Vector{Branch}

  function Contour(branches::AbstractVector{Branch})
    contour_enum = get_contour_enum(branches)
    new(contour_enum, branches)
  end

end

function Contour(d::ContourEnum; tmax=0.0, β=0.0)
  if d == full_contour
    branches = [Branch(forward_branch, tmax), Branch(backward_branch, tmax), Branch(imaginary_branch, β)]
  elseif d == keldysh_contour
    branches = [Branch(forward_branch, tmax), Branch(backward_branch, tmax)]
  else
    branches = [Branch(imaginary_branch, β)]
  end
  return Contour(branches)
end

function twist(c::Contour)
  return Contour(circshift(c.branches, -1))
end

function get_branch(c::Contour, d::BranchEnum)
  for b in c.branches
    b.domain == d && return b
  end
  return nothing
end

"""
θ(t1, t2) = t1 >= t2 ? 1.0 : 0.0
"""
function heaviside(c::Contour, t1::BranchPoint, t2::BranchPoint)
  if t1.domain == t2.domain
    return t1.ref >= t2.ref
  else
    for b in c.branches
      b.domain == t1.domain && return true
      b.domain == t2.domain && return false
    end
  end
end

"""
Heaviside function on standard keldysh contour
"""
function heaviside(t1::BranchPoint, t2::BranchPoint)
  if t1.domain == t2.domain
    return t1.ref >= t2.ref
  else
    return t1.domain > t2.domain
  end
end

const θ = heaviside
