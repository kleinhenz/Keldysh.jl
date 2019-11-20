using LinearAlgebra

struct TimeGF <: AbstractArray{ComplexF64, 2}
  data::Array{ComplexF64, 2}
  grid::TimeGrid
end

function TimeGF(grid::TimeGrid)
  N = length(grid)
  TimeGF(zeros(ComplexF64, N, N), grid)
end

function TimeGF(f::Function, grid::TimeGrid; lower = false, time_invariant = false)
  if time_invariant
    δ = minimum(abs.(grid.step)) / 10

    make_key = (t1, t2) -> begin
      theta = θ(t1.val, t2.val)
      Δt = t1.val.val - t2.val.val
      return (Complex{Int}(round(Δt / δ)), theta)
    end

    T = typeof(f(grid[1], grid[1])) # extra evaluation to get type for cache

    cache = Dict{Tuple{Complex{Int}, Bool}, T}()

    g = (t1, t2) -> begin
      key = make_key(t1, t2)
      key ∈ keys(cache) ? cache[key] : cache[key] = f(t1, t2)
    end

    return TimeGF(g, grid, lower=lower, time_invariant=false)
  else
    N = length(grid)
    gf = TimeGF(grid)

    for t1 in grid
      for t2 in grid
        lower && t1.idx < t2.idx && continue
        gf[t1, t2] = f(t1, t2)
      end
    end
    return gf
  end
end

### AbstractArray Interface ###
IndexStyle(::Type{<:TimeGF}) = IndexLinear()
size(gf::TimeGF) = size(gf.data)
Base.@propagate_inbounds getindex(gf::TimeGF, i::Int) = gf.data[i]
Base.@propagate_inbounds setindex!(gf::TimeGF, v, i::Int) = gf.data[i] = v
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

function getindex(gf::TimeGF, component::Symbol)
  grid = gf.grid
  if component == :greater
    return gf[backward_branch, forward_branch]
  elseif component == :lesser
    return gf[forward_branch, backward_branch]
  elseif component == :matsubara
    gf[imaginary_branch, imaginary_branch][:,1]
  elseif component == :retarded
    gfʳ = gf[:greater] - gf[:lesser]
    θ = LowerTriangular(ones(size(gfʳ)...))
    return θ .* gfʳ
  else
    throw(ArgumentError("component $component not recognized"))
  end
end

"""
n(t) = <c†(t)c(t)> = -i G<(t,t)
"""
function density(gf::TimeGF)
  return -1.0im * diag(gf[:lesser])
end

"""
A(ω) Im -1/π ∫dt Gʳ(t, 0) exp(iωt)
"""
function equilibrium_spectrum(gf::TimeGF, ω)
  grid = gf.grid
  t = map(t -> real(t.val.val), grid[forward_branch])
  trapz = ft -> sum((ft[i] + ft[i+1]) * (t[i+1] - t[i]) / 2 for i in 1:(length(t) - 1))
  gfʳ = gf[:retarded][:,1]
  Aω = [imag((-1.0 / π) * trapz(gfʳ .* exp.(1.0im .* ωi .* t))) for ωi in ω]
end
