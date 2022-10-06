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

function get_ref(c::AbstractContour, t::BranchPoint)
  ref = 0
  for b in c.branches
      lb = length(b)
      if t.domain == b.domain
          return ref + (t.ref * lb)
      else
          ref += lb
      end
  end
  @assert false
end

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
