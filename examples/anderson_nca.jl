using LinearAlgebra
import Base.getindex
import Base.to_index

using Keldysh, HDF5

@enum SpinEnum spin_up = UInt8(1) spin_down=UInt8(2)
to_index(A, sp::SpinEnum) = Int(sp)
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
getindex(st::FockState, sp::SpinEnum) = (st.state & UInt8(sp)) > 0

# flip occupation of ith spin
# note sp ∈ [1, 2] so don't need to do (sp << (i - 1))
flip(st::FockState, sp::SpinEnum) = FockState(st.state ⊻ UInt8(sp))

# convert internal state to index
to_index(st::FockState) = Int(st.state + 1)

# Propagation direction for dyson equation
# NOTE symmetric propagation maintains unitarity but is numerically unstable
@enum DysonDirectionEnum forward_prop backward_prop symmetric_prop

"""
Compute the populations (i.e. the diagonal components of the impurity density matrix) from a propagator
"""
function populations(p)
  grid = p[1].grid
  nstates = length(p)
  ξ = [1.0, -1.0, -1.0, 1.0]
  p_lsr_diag = reduce(hcat, (1.0im * ξ[s] * diag(p[s][:lesser]) for s in 1:nstates))
  Zt = sum(p_lsr_diag, dims=2)
  ρt = p_lsr_diag ./ Zt
  t = realtimes(grid)
  return t, ρt, Zt
end

struct nca_params
  dyson_tol::Float64
  dyson_dir::DysonDirectionEnum
  dyson_max_iter::Int
  max_order::Int
  function nca_params(dyson_tol, dyson_dir, dyson_max_iter, max_order)
    @assert 1 <= max_order <= 2
    new(dyson_tol, dyson_dir, dyson_max_iter, max_order)
  end
end
nca_params(; dyson_tol = 1e-6, dyson_dir = forward_prop, dyson_max_iter = 100, max_order = 1) = nca_params(dyson_tol, dyson_dir, dyson_max_iter, max_order)

struct nca_data
  p0::Array{TimeGF,1} # bare propagator
  Δ::Array{TimeGF, 1} # hybridization function

  p::Array{TimeGF,1} # dressed propagator
  Σ::Array{TimeGF,1} # self-energy
  Σxp::Array{TimeGF,1} # self-energy convolved with propagator
  pxΣ::Array{TimeGF,1} # propagator convolved with self-energy
  G::Array{TimeGF,1} # green's function

  grid::TimeGrid # time grid
  states::NTuple{4, FockState}
  spins::Tuple{SpinEnum, SpinEnum}

  function nca_data(p0, Δ)
    states = ntuple(i -> FockState(i-1), 4)
    spins = instances(SpinEnum)

    statesize = length(states)
    indexsize = length(spins)

    @assert length(p0) == statesize
    @assert length(Δ) == indexsize

    grid = first(p0).grid

    p = [TimeGF(grid) for _ in 1:statesize]
    Σ = [TimeGF(grid) for _ in 1:statesize]
    Σxp = [TimeGF(grid) for _ in 1:statesize]
    pxΣ = [TimeGF(grid) for _ in 1:statesize]
    G = [TimeGF(grid) for _ in 1:indexsize]

    new(p0, Δ, p, Σ, Σxp, pxΣ, G, grid, states, spins)
  end
end

function Σnca(data::nca_data, t1::TimeGridPoint, t2::TimeGridPoint, st_sigma::FockState)
  sum(data.spins) do sp
    st_prop = flip(st_sigma, sp)
    1.0im * data.p[st_prop][t1, t2] * (st_sigma[sp] ? data.Δ[sp][t1, t2] : -data.Δ[sp][t2, t1, false])
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
function Σoca(data::nca_data, t1::TimeGridPoint, t4::TimeGridPoint, st_sigma::FockState)
  Δ, p, grid = data.Δ, data.p, data.grid

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
    # FIXME workaround for failure to elide fancier indexing (seems to be problem in julia 1.5)
    f = t2 -> h1(t2) * integrate(t3 -> h0(t3) * p[st1].data[t2.idx, t3.idx] * p[st2].data[t3.idx, t4.idx], grid, t2, t4)
#    f = t2 -> h1(t2) * integrate(t3 -> h0(t3) * p[st1][t2, t3] * p[st2][t3, t4], grid, t2, t4)

    # integrate over t2
    return integrate(t2 -> p[st0][t1, t2] * f(t2), grid, t1, t4)
  end
end

function dyson(data::nca_data, t1::TimeGridPoint, t2::TimeGridPoint, params::nca_params)
  @assert t1.idx >= t2.idx

  p_t1t2_cur = zeros(eltype(data.p[1]), length(data.states))
  p_t1t2_next = zeros(eltype(data.p[1]), length(data.states))

  for st in data.states
    p_t1t2_cur[st] = data.p0[st][t1,t2]
    data.p[st][t1,t2] = data.p0[st][t1,t2] # initial guess
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
      if params.dyson_dir == forward_prop || params.dyson_dir == symmetric_prop
        data.Σxp[st][t1, t2] = data.Σ[st] ↻ data.p[st]
        p_t1t2_next[st] += data.p0[st] ↻ data.Σxp[st]
      end

      # p = p₀ + p ↻ Σ ↻ p₀
      if params.dyson_dir == backward_prop || params.dyson_dir == symmetric_prop
        data.pxΣ[st][t1, t2] = data.p[st] ↻ data.Σ[st]
        p_t1t2_next[st] += data.pxΣ[st] ↻ data.p0[st]
      end

      params.dyson_dir == symmetric_prop && (p_t1t2_next[st] *= 0.5)

      p_t1t2_next[st] += data.p0[st][t1, t2]
    end

    diff = norm(p_t1t2_cur - p_t1t2_next)
    done = diff < params.dyson_tol * norm(p_t1t2_cur)
    for st in data.states
      data.p[st][t1,t2] = p_t1t2_next[st]
    end
    p_t1t2_cur .= p_t1t2_next
    iter += 1
  end
end

function nca(p0, Δ, params::nca_params)
  data = nca_data(p0, Δ)
  N = length(data.grid)
  for d in 0:(N-1) # solve diagonal by diagonal
    for j in 1:(N-d)
      i = j + d
      t1 = data.grid[i]
      t2 = data.grid[j]
      dyson(data, t1, t2, params)
    end
  end

  return data
end

"""
    run_anderson_nca(;tmax=5.0, npts_real = 51, β = 1.0, dos = Keldysh.flat_dos(), U = 4.0, ρ0 = [1.0, 0.0, 0.0, 0.0], tol=1e-6)

Solve the Anderson impurity model on the two-branch Keldysh contour

# Arguments
- `tmax::Real=5.0`: the maximum time on the contour
- `npts_real::Int=51`: the number of points used to discretize the forward/backward branches of the contour
- `β::Real=1.0`: the inverse temperature of the bath
- `dos=Keldysh.flat_dos()`: the bath density of states
- `U::Real=4.0`: the coulomb interaction strength
- `ρ0::Array{Float,1}=[1.0, 0.0, 0.0, 0.0]`: the initial (diagonal) impurity density matrix
- `params::nca_params`: parameters controlling solution of nca equations
"""
function run_anderson_nca(;tmax=5.0, npts_real = 51, β = 1.0, dos = Keldysh.flat_dos(), U = 4.0, ρ0 = [1.0, 0.0, 0.0, 0.0], params::nca_params = nca_params())
  @assert length(ρ0) == 4

  c = twist(Contour(keldysh_contour, tmax=tmax))
  grid = TimeGrid(c, npts_real = npts_real)
  ϵ = [0.0, -U/2, -U/2, 0.0]

  ξ = [1.0, -1.0, -1.0, 1.0]

  p0 = map(1:4) do s
    TimeGF(grid, lower=true) do t1, t2
      val = -1.0im * exp(-1.0im * (t1.val.val - t2.val.val) * ϵ[s])
      heaviside(t1.val, t2.val) || (val *= ξ[s] * ρ0[s])
      return val
    end
  end

  Δup = dos2gf(dos, β, grid)
  Δ = [deepcopy(Δup) for _ in 1:2] # spin symmetric

  data = nca(p0, Δ, params)

  return data
end

function main()
  dos = Keldysh.flat_dos(ν=10.0, D=10.0)
  params = nca_params(max_order = 1)
  data = run_anderson_nca(β=1.0, tmax=5.0, npts_real = 101, U = 8.0, dos=dos, params=params)
  t, ρt, Zt = populations(data.p)
  h5write("output.h5", "output/rho", ρt)
  h5write("output.h5", "output/Z", ρt)
  h5write("output.h5", "output/t", t)
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
