@enum ContourEnum full_contour=1 keldysh_contour=2 imaginary_contour=3

struct Contour
  domain::ContourEnum
  branches::Vector{Branch}
  branch_indices::NTuple{3, Int} # branches[branch_indices[Int(BranchEnum)]].domain == BranchEnum

  function Contour(branches::AbstractVector{Branch})
    contour_enum = get_contour_enum(branches)
    branch_indices::Vector{Int} = [-1, -1, -1]
    for (i, b) in enumerate(branches)
      branch_indices[Int(b.domain)] = i
    end
    new(contour_enum, branches, Tuple(branch_indices))
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

function nbranches(c::ContourEnum)
  if c == full_contour
    return 3
  elseif c == keldysh_contour
    return 2
  else
    return 1
  end
end

function nbranches(c::Contour)
  return nbranches(c.domain)
end

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
  branch_enums = Set([b.domain for b in branches])
  @assert length(branch_enums) == length(branches) "branches must be unique"
  return get_contour_enum(branch_enums)
end

function twist(c::Contour)
  return Contour(circshift(c.branches, -1))
end

function get_branch(c::Contour, d::BranchEnum)
  idx = c.branch_indices[Int(d)]
  return idx == -1 ? nothing : c.branches[idx]
end

Base.getindex(c::Contour, d::BranchEnum) = get_branch(c, d)

"""
θ(t1, t2) = t1 >= t2 ? 1.0 : 0.0
"""
function heaviside(c::Contour, t1::BranchPoint, t2::BranchPoint)
  t1.domain == t2.domain ? t1.ref >= t2.ref : c.branch_indices[Int(t1.domain)] > c.branch_indices[Int(t2.domain)]
end

"""
Heaviside function on standard keldysh contour
"""
function heaviside(t1::BranchPoint, t2::BranchPoint)
  t1.domain == t2.domain ?  t1.ref >= t2.ref : t1.domain > t2.domain
end

const θ = heaviside

function Base.in(b::BranchEnum, c::ContourEnum)
  if c == keldysh_contour
    return b == imaginary_branch ? false : true
  elseif c == imaginary_branch
    return b == imaginary ? true : false
  else
    return true
  end
end

function Base.in(b::BranchEnum, c::Contour)
  return in(b, c.domain)
end
