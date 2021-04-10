using LinearAlgebra

using Keldysh, HDF5

parse_param(::Type{T}, s::AbstractString) where T = parse(T, s)
parse_param(::Type{String}, s::AbstractString) = string(s)
function parse_params(args, param_def)
  N = length(param_def)

  names = ntuple(N) do i
    param_def[i][2]
  end

  vals = map(param_def) do (type, name, default)
    i = findfirst(isequal("--$name"), args)
    return isnothing(i) ? default : parse_param(type, args[i+1])
  end

  return NamedTuple{names}(vals)
end

@enum SpinEnum spin_up = UInt8(1) spin_down=UInt8(2)
Base.to_index(A, sp::SpinEnum) = Int(sp)
flip(sp::SpinEnum) = sp == spin_up ? spin_down : spin_up

# type representing state of anderson impurity
# occupation of up/down stored in first two bits
struct FockState
    state::UInt8
    function FockState(s)
        @assert s < 4
        new(s)
    end
end

# check whether ith spin is occupied
# note sp ∈ [1, 2] so don't need to do (1 << (sp - 1))
Base.getindex(st::FockState, sp::SpinEnum) = (st.state & UInt8(sp)) > 0

# flip occupation of ith spin
# note sp ∈ [1, 2] so don't need to do (sp << (i - 1))
flip(st::FockState, sp::SpinEnum) = FockState(st.state ⊻ UInt8(sp))

# convert internal state to index
Base.to_index(st::FockState) = Int(st.state + 1)

"""
Compute the populations (i.e. the diagonal components of the impurity density matrix) from a propagator
"""
function populations(P)
  nstates = length(P)
  ξ = [1.0, -1.0, -1.0, 1.0]
  P_lsr_diag = reduce(hcat, (1.0im * ξ[s] * diag(P[s][:lesser]) for s in 1:nstates))
  Zt = sum(P_lsr_diag, dims=2)
  ρt = P_lsr_diag ./ Zt
  return ρt, Zt
end

struct NCAParams
  dyson_rtol::Float64
  dyson_atol::Float64
  dyson_max_iter::Int
  max_order::Int
  function NCAParams(dyson_rtol, dyson_atol, dyson_max_iter, max_order)
    @assert 1 <= max_order <= 2
    new(dyson_rtol, dyson_atol, dyson_max_iter, max_order)
  end
end
NCAParams(; dyson_rtol = 1e-6, dyson_atol = 1e-10, dyson_max_iter = 100, max_order = 1) = NCAParams(dyson_rtol, dyson_atol, dyson_max_iter, max_order)

struct NCAData{T <: AbstractTimeGF, U <: AbstractTimeGrid}
  P0::Array{T,1} # bare propagator
  Δ::Array{T, 1} # hybridization function

  P::Array{T,1} # dressed propagator
  Σ::Array{T,1} # self-energy
  ΣxP::Array{T,1} # self-energy convolved with propagator
  G::Array{T,1} # green's function

  grid::U # time grid
  states::NTuple{4, FockState}
  spins::Tuple{SpinEnum, SpinEnum}
end

function NCAData(P0, Δ)
  states = ntuple(i -> FockState(i-1), 4)
  spins = instances(SpinEnum)

  statesize = length(states)
  indexsize = length(spins)

  @assert length(P0) == statesize
  @assert length(Δ) == indexsize

  grid = first(P0).grid

  X = P0[1]
  P = [zero(X) for _ in 1:statesize]
  Σ = [zero(X) for _ in 1:statesize]
  ΣxP = [zero(X) for _ in 1:statesize]
  G = [zero(X) for _ in 1:indexsize]

  NCAData(P0, Δ, P, Σ, ΣxP, G, grid, states, spins)
end

function Σnca(data::NCAData, t1::TimeGridPoint, t2::TimeGridPoint, st_sigma::FockState)
  sum(data.spins) do sp
    st_prop = flip(st_sigma, sp)
    1.0im * data.P[st_prop][t1, t2] * (st_sigma[sp] ? data.Δ[sp][t1, t2] : -data.Δ[sp][t2, t1, false])
  end
end

# collect all one-crossing terms
# define 4 points t1 > t2 > t3 > t4
#           _________
#           |        |
#      __________    |
# _____|    |    |___|___
#     t1   t2    t3  t4
# <------------------------
# t
# st_sigma <- st0 <- st1 <- st2 <- st_sigma
#
# delta[sp0](t1,t3) * delta[sp1](t2, t4)
function Σoca(data::NCAData, t1::TimeGridPoint, t4::TimeGridPoint, st_sigma::FockState)
  Δ, P, grid = data.Δ, data.P, data.grid

  sum(data.spins) do sp1
    sp0 = flip(sp1)

    st2 = flip(st_sigma, sp1)
    st1 = flip(st2, sp0)
    st0 = flip(st1, sp1)

    # delta[sp1](t2, t4)
    h1 = t2 -> 1.0im * (st_sigma[sp1] ? Δ[sp1][t2, t4] : -Δ[sp1][t4, t2, false])

    # delta[sp0](t1, t3)
    h0 = t3 -> 1.0im * (st_sigma[sp0] ? Δ[sp0][t1, t3] : -Δ[sp0][t3, t1, false])

    # integrate over t3
    f = t2 -> h1(t2) * integrate(t3 -> h0(t3) * P[st1][t2, t3] * P[st2][t3, t4], grid, t2, t4)

    # integrate over t2
    return integrate(t2 -> P[st0][t1, t2] * f(t2), grid, t1, t4)
  end
end

function dyson!(data::NCAData, t1::TimeGridPoint, t2::TimeGridPoint, params::NCAParams)
  @assert t1.idx >= t2.idx

  p_t1t2_cur = zeros(ComplexF64, length(data.states))
  p_t1t2_next = zeros(ComplexF64, length(data.states))

  for st in data.states
    p_t1t2_cur[st] = data.P0[st][t1,t2]
    data.P[st][t1,t2] = data.P0[st][t1,t2] # initial guess
  end

  ↻ = (A, B) -> integrate(t -> @inbounds(A[t1, t] * B[t, t2]), data.grid, t1, t2)

  done = false
  iter = 1
  diff = 0.0
  while iter <= params.dyson_max_iter && !done

    for st in data.states
      p_t1t2_next[st] = 0.0

      data.Σ[st][t1, t2] = 0.0
      params.max_order > 0 && (data.Σ[st][t1, t2] += Σnca(data, t1, t2, st))
      params.max_order > 1 && (data.Σ[st][t1, t2] += Σoca(data, t1, t2, st))

      # p = p₀ + p₀ ↻ Σ ↻ p
      data.ΣxP[st][t1, t2] = data.Σ[st] ↻ data.P[st]
      p_t1t2_next[st] += data.P0[st] ↻ data.ΣxP[st]

      p_t1t2_next[st] += data.P0[st][t1, t2]
    end

    diff = norm(p_t1t2_cur - p_t1t2_next)
    done = diff < max(params.dyson_atol, params.dyson_rtol * norm(p_t1t2_cur))
    for st in data.states
      data.P[st][t1,t2] = p_t1t2_next[st]
    end
    p_t1t2_cur .= p_t1t2_next
    iter += 1
  end
end

function nca!(data::NCAData, params::NCAParams)
  N = length(data.grid)
  for d in 0:(N-1) # solve diagonal by diagonal
    println("diagonal $(d+1)/$N")
    for j in 1:(N-d)
      i = j + d
      t1 = data.grid[i]
      t2 = data.grid[j]
      dyson!(data, t1, t2, params)
    end
  end

  return data
end

function make_bare_prop(grid::KeldyshTimeGrid, ρ, ϵ, U)
  E = [0.0, ϵ, ϵ, 2*ϵ + U]
  ξ = [1.0, -1.0, -1.0, 1.0]
  P = map(1:4) do s
    GenericTimeGF(grid, 1, true) do t1, t2
      t1.idx < t2.idx && return 0.0
      ϕ = integrate(t -> E[s], grid, t1, t2)
      heaviside(t1.val, t2.val) ? -im * exp(-im * ϕ) : -im * ξ[s] * ρ[s] * exp(-im * ϕ)
    end
  end
  return P
end

function make_bare_prop(grid::FullTimeGrid, ϵ, U)
  E = [0.0, ϵ, ϵ, 2*ϵ + U]
  ξ = [1.0, -1.0, -1.0, 1.0]
  P = map(1:4) do s
    GenericTimeGF(grid, 1, true) do t1, t2
      t1.idx < t2.idx && return 0.0
      ϕ = integrate(t -> E[s], grid, t1, t2)
      heaviside(t1.val, t2.val) ? -im * exp(-im * ϕ) : -im * ξ[s] * exp(-im * ϕ)
    end
  end
  return P
end

param_def = [(String, :contour, "keldysh"),
             (Float64, :tmax, 10.0),
             (Float64, :beta, 5.0),
             (Int, :nt, 201),
             (Int, :ntau, 101),
             (Float64, :rtol, 1e-6),
             (Float64, :atol, 1e-10),
             (Int, :max_iter, 200),
             (String, :mode, "nca"),
             (Float64, :D, 10.0),
             (Float64, :eps, -3.0),
             (Float64, :U, 8.0),
             (String, :output_file, "output.h5")]

function main()
  println("running nca...")

  p = parse_params(ARGS, param_def)

  for (k,v) in pairs(p)
    println("$k : $v")
  end

  data =
  if p.contour == "keldysh"
    c = twist(KeldyshContour(tmax=p.tmax))
    grid = KeldyshTimeGrid(c, p.nt)

    ρ = [1.0, 0.0, 0.0, 0.0]
    P0 = make_bare_prop(grid, ρ, p.eps, p.U)

    dos = Keldysh.flat_dos(ν=10.0, D=p.D)
    Δ = [GenericTimeGF(dos, p.beta, grid) for s in 1:2]

    NCAData(P0, Δ)
  elseif p.contour == "full"
    c = twist(FullContour(tmax=p.tmax, β=p.beta))
    grid = FullTimeGrid(c, p.nt, p.ntau)
    P0 = make_bare_prop(grid, p.eps, p.U)

    dos = Keldysh.flat_dos(ν=10.0, D=p.D)
    Δ = [GenericTimeGF(dos, grid) for s in 1:2]

    NCAData(P0, Δ)
  else
    error("unknown contour $(p.contour)")
  end

  max_order = Dict("nca"=>1, "oca"=>2)[p.mode]
  params = NCAParams(dyson_rtol = p.rtol, dyson_atol = p.atol, dyson_max_iter = p.max_iter, max_order = max_order)

  nca!(data, params)

  t = collect(realtimes(data.grid))
  ρt, Zt = populations(data.P)

  h5open(p.output_file, "w") do h5f
    for (k,v) in pairs(p)
      write(h5f, "/input/params/$k", v)
    end

    h5f["/output/obs/pop/rho"] = ρt
    h5f["/output/obs/pop/Z"] = Zt
    h5f["/output/obs/pop/t"] = collect(t)
  end
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
