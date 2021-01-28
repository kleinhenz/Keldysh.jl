using LinearAlgebra

abstract type TimeGF <: AbstractArray{ComplexF64, 2} end

struct TimeScalarGF <: TimeGF
  data::Array{ComplexF64, 2}
  grid::TimeGrid
end

struct TimeMatrixGF <: TimeGF
  data::Array{ComplexF64, 4}
  grid::TimeGrid
end

function TimeGF(data::Array{ComplexF64, N}, grid::TimeGrid) where {N}
    if N == 2
        return TimeScalarGF(data, grid)
    elseif N == 4
        return TimeMatrixGF(data, grid)
    else
        throw(ArgumentError())
    end
end

function TimeGF(grid::TimeGrid, Norb = nothing)
  N = length(grid)
  if Norb == nothing
      data = zeros(ComplexF64, N, N)
  else
      data = zeros(ComplexF64, Norb, Norb, N, N)
  end
  TimeGF(data, grid)
end

function TimeGF(f::Function, grid::TimeGrid; Norb = nothing, lower = false)
  N = length(grid)
  gf = TimeGF(grid, Norb)
    
  for t1 in grid
    for t2 in grid
      lower && t1.idx < t2.idx && continue
      gf[t1, t2] = f(t1, t2)
    end
  end
  return gf
end

function TimeGF(les::AbstractMatrix, ret::AbstractMatrix, tv::AbstractMatrix, mat::AbstractVector, grid::TimeGrid)
  @assert grid.contour.domain == full_contour

  t = realtimes(grid)
  tau = imagtimes(grid)

  nt = length(t)
  ntau = length(tau)

  @assert size(les) == (nt,nt)
  @assert size(ret) == (nt,nt)
  @assert size(tv) == (nt,ntau)
  @assert size(mat) == (ntau,)

  G = TimeGF(grid) do t1, t2
    greater = heaviside(t1.val, t2.val)
    if ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) &&
        (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
      i = findfirst(isapprox(real(t1.val.val)), t)
      j = findfirst(isapprox(real(t2.val.val)), t)
      greater ? ret[i,j] + les[i,j] : les[i,j]
    elseif (t1.val.domain == imaginary_branch && (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
      i = findfirst(isapprox(imag(t1.val.val)), -tau)
      j = findfirst(isapprox(real(t2.val.val)), t)
      conj(tv[j, ntau+1-i]) # akoi 19c
    elseif ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) && t2.val.domain == imaginary_branch)
      i = findfirst(isapprox(real(t1.val.val)), t)
      j = findfirst(isapprox(imag(t2.val.val)), -tau)
      tv[i,j]
    else
      i = findfirst(isapprox(imag(t1.val.val)), -tau)
      j = findfirst(isapprox(imag(t2.val.val)), -tau)
      greater ? 1.0im * mat[i - j + 1] : -1.0im * mat[i - j + ntau]
    end
  end

  return G
end

function TimeGF(les::AbstractMatrix, ret::AbstractMatrix, grid::TimeGrid)
  @assert grid.contour.domain == keldysh_contour

  t = realtimes(grid)

  nt = length(t)

  @assert size(les) == (nt,nt)
  @assert size(ret) == (nt,nt)

  G = TimeGF(grid) do t1, t2
    greater = heaviside(t1.val, t2.val)
    i = findfirst(isapprox(real(t1.val.val)), t)
    j = findfirst(isapprox(real(t2.val.val)), t)
    greater ? ret[i,j] + les[i,j] : les[i,j]
  end

  return G
end

### AbstractArray Interface ###
IndexStyle(::Type{<:TimeGF}) = IndexLinear()
size(gf::TimeGF) = size(gf.data)
Base.@propagate_inbounds getindex(gf::TimeGF, i::Int) = gf.data[i]
Base.@propagate_inbounds setindex!(gf::TimeGF, v, i::Int) = gf.data[i] = v

Base.@propagate_inbounds Base.getindex(gf::TimeMatrixGF, i::Int, j::Int) = view(gf.data, :, :, i, j)
Base.@propagate_inbounds Base.setindex!(gf::TimeMatrixGF, v, i::Int, j::Int) = gf.data[:, :, i, j] = v


function similar(gf::TimeGF, ::Type{S}) where S
  data = similar(gf.data, S)
  TimeGF(data, gf.grid) # TODO copy grid?
end

### BroadCasting ###
Base.BroadcastStyle(::Type{<:TimeGF}) = Broadcast.ArrayStyle{TimeGF}()
function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{TimeGF}}, ::Type{ElType}) where ElType
    # Scan the inputs for first TimeGF
    gf = find_time_gf(bc)
    # copy grid to output
    TimeGF(similar(Array{ElType}, axes(bc)), gf.grid) # TODO copy grid?
end

"`A = find_time_gf(As)` returns the first TimeGF among the arguments."
find_time_gf(bc::Base.Broadcast.Broadcasted) = find_time_gf(bc.args)
find_time_gf(args::Tuple) = find_time_gf(find_time_gf(args[1]), Base.tail(args))
find_time_gf(x) = x
find_time_gf(a::TimeGF, rest) = a
find_time_gf(::Any, rest) = find_time_gf(rest)

# indexing with TimeGridPoint
Base.@propagate_inbounds function getindex(gf::TimeGF, t1::TimeGridPoint, t2::TimeGridPoint, gtr=true)
  val = gf[t1.idx, t2.idx]
  (!gtr && t1.idx == t2.idx) && (val += jump(gf))
  return val
end
Base.@propagate_inbounds setindex!(gf::TimeGF, v, t1::TimeGridPoint, t2::TimeGridPoint) = gf[t1.idx, t2.idx] = v

function jump(gf::TimeGF)
  t0_plus = branch_bounds(gf.grid, forward_branch)[1]
  t0_minus = branch_bounds(gf.grid, backward_branch)[2]
  return gf[t0_plus, t0_minus] - gf[t0_plus, t0_plus]
end

"""
Transposed view of gf taking into account diagonal discontinuity
"""
struct TimeGFTranspose <: AbstractArray{ComplexF64, 2}
  gf::TimeGF
end

adjoint(gf::TimeGF) = TimeGFTranspose(gf)

### AbstractArray Interface ###
IndexStyle(::Type{<:TimeGFTranspose}) = IndexCartesian()
size(gfa::TimeGFTranspose) = size(gfa.gf)
Base.@propagate_inbounds getindex(gfa::TimeGFTranspose, i::Int, j::Int) = i == j ? gfa.gf[i,i] + jump(gfa.gf) : gfa.gf[j,i]
Base.@propagate_inbounds setindex!(gfa::TimeGFTranspose, v, i::Int, j::Int) = gfa.gf[j,i] = v

# indexing with TimeGridPoint
Base.@propagate_inbounds getindex(gfa::TimeGFTranspose, t1::TimeGridPoint, t2::TimeGridPoint) = gfa[t1.idx, t2.idx]
Base.@propagate_inbounds setindex!(gfa::TimeGFTranspose, v::ComplexF64, t1::TimeGridPoint, t2::TimeGridPoint) = gfa[t1.idx, t2.idx] = v

function getindex(gf::TimeGF, b1::BranchEnum, b2::BranchEnum)
  grid = gf.grid
  @assert b1 ∈ grid.contour && b2 ∈ grid.contour
  x = b1 == backward_branch ? reverse(grid[b1]) : grid[b1]
  y = b2 == backward_branch ? reverse(grid[b2]) : grid[b2]
  return [gf[t1, t2] for t1 in x, t2 in y]
end

function getindex(G::TimeGF, component::Symbol)
  grid = G.grid
  if component == :greater
    return G[backward_branch, forward_branch]
  elseif component == :lesser
    return G[forward_branch, backward_branch]
  elseif component == :matsubara
    real(-im .* G[imaginary_branch, imaginary_branch][:,1])
  elseif component == :retarded
    ret = G[:greater] - G[:lesser]
    θ = LowerTriangular(ones(size(ret)...))
    return θ .* ret
  elseif component == :advanced
    adv = G[:lesser] - G[:greater]
    θ = UpperTriangular(ones(size(adv)...))
    return θ .* adv
  elseif component == :leftmixing
    gf[backward_branch, imaginary_branch]
  else
    throw(ArgumentError("component $component not recognized"))
  end
end

"""
    density(gf::TimeGF)

n(t) = <c†(t)c(t)> = -i G<(t,t)
"""
function density(gf::TimeGF)
  return -1.0im * diag(gf[:lesser])
end

"""
    equilibrium_spectrum(gf::TimeGF, ω)

A(ω) = Im -1/π ∫dt Gʳ(t, 0) exp(iωt)
"""
function equilibrium_spectrum(gf::TimeGF, ω)
  grid = gf.grid
  t = realtimes(grid)
  trapz = ft -> sum((ft[i] + ft[i+1]) * (t[i+1] - t[i]) / 2 for i in 1:(length(t) - 1))
  gfʳ = gf[:retarded][:,1]
  Aω = [imag((-1.0 / π) * trapz(gfʳ .* exp.(1.0im .* ωi .* t))) for ωi in ω]
end

function aux_spectrum(gf::TimeGF, ω::Number)
  grid = gf.grid

  delta_f = (t1, t2) -> -1.0im *(θ(t1.val, t2.val) - 1) * exp(-1.0im * (t1.val.val - t2.val.val) * ω)
  delta_e = (t1, t2) -> -1.0im *(θ(t1.val, t2.val) - 0) * exp(-1.0im * (t1.val.val - t2.val.val) * ω)

  At = map(zip(grid[forward_branch], reverse(grid[backward_branch]))) do (t⁺, t⁻)
    Iω_f = -2.0 * integrate(t -> @inbounds(gf[t⁺, t] * delta_f(t, t⁻)), grid, t⁺, t⁻)
    Iω_e = 2.0 * integrate(t -> @inbounds(delta_e(t⁺, t) * gf[t, t⁻]), grid, t⁺, t⁻)
    (1 / (2π)) * (Iω_e - Iω_f)
  end
end

"""
    aux_spectrum(gf::TimeGF, ω)

Compute auxiliary current spectral function of gf for frequencies ω
See Eq. 24 from [1]

[1] Cohen, Guy, David R. Reichman, Andrew J. Millis, and Emanuel Gull. “Green’s
Functions from Real-Time Bold-Line Monte Carlo.” Physical Review B 89, no. 11
(March 31, 2014): 115139. https://doi.org/10.1103/PhysRevB.89.115139.
"""
function aux_spectrum(gf::TimeGF, ω)
  reduce(vcat, (aux_spectrum(gf, ωi)' for ωi in ω))
end
