abstract type AbstractContour end

struct ImaginaryContour <: AbstractContour
  branches::NTuple{1, Branch}
  β::Float64
end
function ImaginaryContour(;β)
  branches = (Branch(imaginary_branch, β),)
  return ImaginaryContour(branches, β)
end

struct KeldyshContour <: AbstractContour
  branches::NTuple{2, Branch}
  tmax::Float64
end
function KeldyshContour(; tmax)
  branches = (Branch(forward_branch, tmax), Branch(backward_branch, tmax))
  return KeldyshContour(branches, tmax)
end

struct FullContour <: AbstractContour
  branches::NTuple{3, Branch}
  tmax::Float64
  β::Float64
end
function FullContour(; tmax, β)
  branches = (Branch(forward_branch, tmax), Branch(backward_branch, tmax), Branch(imaginary_branch, β))
  return FullContour(branches, tmax, β)
end

nbranches(::Type{FullContour}) = 3
nbranches(::Type{KeldyshContour}) = 2
nbranches(::Type{ImaginaryContour}) = 1
nbranches(c::AbstractContour) = nbranches(typeof(c))

Base.length(c::ImaginaryContour) = c.β
Base.length(c::KeldyshContour) = 2*c.tmax
Base.length(c::FullContour) = 2*c.tmax + c.β

function get_point(c::AbstractContour, ref)
  @assert ref >= 0.0 && ref <= length(c)

  for b in c.branches
    lb = length(b)
    if ref <= lb
      return b(ref / lb)
    else
      ref -= lb
    end
  end
end

(c::AbstractContour)(ref::Real) = get_point(c, ref)

function twist(c::FullContour)
  n = nbranches(c)
  FullContour(ntuple(i -> c.branches[mod1(i+1,n)], n), c.tmax, c.β)
end

function twist(c::KeldyshContour)
  n = nbranches(c)
  KeldyshContour(ntuple(i -> c.branches[mod1(i+1,n)], n), c.tmax)
end

"""
θ(t1, t2) = t1 >= t2 ? 1.0 : 0.0
"""
function heaviside(c::AbstractContour, t1::BranchPoint, t2::BranchPoint)
  if t1.domain == t2.domain
    return t1.ref >= t2.ref
  else
    return findfirst(b -> b.domain == t1.domain, c.branches) >
           findfirst(b -> b.domain == t2.domain, c.branches)
  end
end

"""
Heaviside function on standard keldysh contour
"""
function heaviside(t1::BranchPoint, t2::BranchPoint)
  t1.domain == t2.domain ?  t1.ref >= t2.ref : t1.domain > t2.domain
end

function Base.getindex(c::AbstractContour, d::BranchEnum)
  idx = findfirst(b -> b.domain == d, c.branches)
  return isnothing(idx) ? nothing : c.branches[idx]
end

Base.in(b::BranchEnum, ::Type{FullContour}) = true
Base.in(b::BranchEnum, ::Type{KeldyshContour}) = b == imaginary_branch ? false : true
Base.in(b::BranchEnum, ::Type{ImaginaryContour}) = b == imaginary_branch ? true : false
Base.in(b::BranchEnum, c::AbstractContour) = in(b, typeof(c))

#function Base.in(b::BranchEnum, c::ContourEnum)
#  if c == keldysh_contour
#    return b == imaginary_branch ? false : true
#  elseif c == imaginary_branch
#    return b == imaginary ? true : false
#  else
#    return true
#  end
#end

#struct Contour
#  domain::ContourEnum
#  branches::Vector{Branch}
#  branch_indices::NTuple{3, Int} # branches[branch_indices[Int(BranchEnum)]].domain == BranchEnum
#
#  function Contour(branches::AbstractVector{Branch})
#    contour_enum = get_contour_enum(branches)
#    branch_indices::Vector{Int} = [-1, -1, -1]
#    for (i, b) in enumerate(branches)
#      branch_indices[Int(b.domain)] = i
#    end
#    new(contour_enum, branches, Tuple(branch_indices))
#  end
#
#end
#
#function Contour(d::ContourEnum; tmax=0.0, β=0.0)
#  if d == full_contour
#    branches = [Branch(forward_branch, tmax), Branch(backward_branch, tmax), Branch(imaginary_branch, β)]
#  elseif d == keldysh_contour
#    branches = [Branch(forward_branch, tmax), Branch(backward_branch, tmax)]
#  else
#    branches = [Branch(imaginary_branch, β)]
#  end
#  return Contour(branches)
#end
#
#function nbranches(c::ContourEnum)
#  if c == full_contour
#    return 3
#  elseif c == keldysh_contour
#    return 2
#  else
#    return 1
#  end
#end
#
#function nbranches(c::Contour)
#  return nbranches(c.domain)
#end
#
#function branch_set(c::ContourEnum)
#  if c == full_contour
#    return Set([forward_branch, backward_branch, imaginary_branch])
#  elseif c == keldysh_contour
#    return Set([forward_branch, backward_branch])
#  else
#    return Set([imaginary_branch])
#  end
#end
#
#function get_contour_enum(branches::Set{BranchEnum})
#  for c in instances(ContourEnum)
#    branch_set(c) == branches && return c
#  end
#  @assert false "invalid branch set"
#end
#
#function get_contour_enum(branches::AbstractVector{Branch})
#  branch_enums = Set([b.domain for b in branches])
#  @assert length(branch_enums) == length(branches) "branches must be unique"
#  return get_contour_enum(branch_enums)
#end
#
#function twist(c::Contour)
#  return Contour(circshift(c.branches, -1))
#end
#
#function get_branch(c::Contour, d::BranchEnum)
#  idx = c.branch_indices[Int(d)]
#  return idx == -1 ? nothing : c.branches[idx]
#end
#
#Base.getindex(c::Contour, d::BranchEnum) = get_branch(c, d)
#
#"""
#θ(t1, t2) = t1 >= t2 ? 1.0 : 0.0
#"""
#function heaviside(c::Contour, t1::BranchPoint, t2::BranchPoint)
#  t1.domain == t2.domain ? t1.ref >= t2.ref : c.branch_indices[Int(t1.domain)] > c.branch_indices[Int(t2.domain)]
#end
#
#"""
#Heaviside function on standard keldysh contour
#"""
#function heaviside(t1::BranchPoint, t2::BranchPoint)
#  t1.domain == t2.domain ?  t1.ref >= t2.ref : t1.domain > t2.domain
#end
#
#const θ = heaviside
#
#function Base.in(b::BranchEnum, c::ContourEnum)
#  if c == keldysh_contour
#    return b == imaginary_branch ? false : true
#  elseif c == imaginary_branch
#    return b == imaginary ? true : false
#  else
#    return true
#  end
#end
#
#function Base.in(b::BranchEnum, c::Contour)
#  return in(b, c.domain)
#end
