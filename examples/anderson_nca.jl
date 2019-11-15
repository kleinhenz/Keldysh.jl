using LinearAlgebra

using Keldysh, HDF5

@enum DysonDirectionEnum forward backward symmetric

"""
Compute the populations (i.e. the diagonal components of the impurity density matrix) from a propagator
"""
function populations(p)
  grid = p[1].grid
  npts_real = length(grid, forward_branch)
  nstates = length(p)
  ξ = [1.0, -1.0, -1.0, 1.0]
  p_lsr_diag = reduce(hcat, (1.0im * ξ[s] * diag(p[s][:lesser]) for s in 1:nstates))
  Zt = sum(p_lsr_diag, dims=2)
  ρt = p_lsr_diag ./ Zt
  t = map(t -> real(t.val.val),  grid[forward_branch])
  return t, ρt, Zt
end

struct nca_params
  dyson_tol::Real
  dyson_dir::DysonDirectionEnum
  dyson_max_iter::Int
end
nca_params(; dyson_tol = 1e-6, dyson_dir = forward, dyson_max_iter = 100) = nca_params(dyson_tol, dyson_dir, dyson_max_iter)

struct nca_data
  p0 # bare propagator
  Δ # hybridization function

  p # dressed propagator
  Σ # self-energy
  Σxp # self-energy convolved with propagator
  pxΣ # self-energy convolved with propagator
  G # green's function

  grid # time grid
  statesize
  indexsize

  function nca_data(p0, Δ)
    statesize = length(p0)
    indexsize = length(Δ)
    @assert statesize == 2^indexsize

    grid = first(p0).grid

    p = [TimeGF(grid) for _ in 1:statesize]
    Σ = [TimeGF(grid) for _ in 1:statesize]
    Σxp = [TimeGF(grid) for _ in 1:statesize]
    pxΣ = [TimeGF(grid) for _ in 1:statesize]
    G = [TimeGF(grid) for _ in 1:indexsize]

    new(p0, Δ, p, Σ, Σxp, pxΣ, G, grid, statesize, indexsize)
  end
end

function Σnca(data::nca_data, t1::TimeGridPoint, t2::TimeGridPoint, st_sigma::UInt64)
  val = 0.0im
  for sp in 1:data.indexsize
    st_prop = ((st_sigma - 1) ⊻ sp) + 1
    gtr = st_sigma & (1 << (sp - 1)) > 0
    val += 1.0im * data.p[st_prop][t1, t2] * (gtr ? data.Δ[sp][t1, t2] : -data.Δ[sp][t2, t1, false])
  end
  return val
end

function dyson(data::nca_data, t1::TimeGridPoint, t2::TimeGridPoint, params::nca_params)
  @assert t1.idx >= t2.idx
  # NOTE symmetric propagation maintains unitarity but is numerically unstable

  x0 = zeros(eltype(data.p[1]), data.statesize)
  x1 = zeros(eltype(data.p[1]), data.statesize)

  for s in 1:data.statesize
    x0[s] = data.p0[s][t1,t2]
    data.p[s][t1,t2] = data.p0[s][t1,t2] # initial guess
  end

  convo_integr = (A, B) -> integrate(t -> A[t1, t] * B[t, t2], data.grid, t1, t2)

  done = false
  iter = 1
  diff = 0.0
  while iter <= params.dyson_max_iter && !done

    for s in 1:data.statesize
      x1[s] = 0.0

      data.Σ[s][t1, t2] = Σnca(data, t1, t2, UInt64(s))
      if params.dyson_dir == forward || params.dyson_dir == symmetric
        data.Σxp[s][t1, t2] = convo_integr(data.Σ[s], data.p[s])
        x1[s] += convo_integr(data.p0[s], data.Σxp[s])
      end

      if params.dyson_dir == backward || params.dyson_dir == symmetric
        data.pxΣ[s][t1, t2] = convo_integr(data.p[s], data.Σ[s])
        x1[s] += convo_integr(data.pxΣ[s], data.p0[s])
      end

      params.dyson_dir == symmetric && (x1[s] *= 0.5)

      x1[s] += data.p0[s][t1, t2]
    end

    diff = norm(x0 - x1)
    done = diff < params.dyson_tol * norm(x0)
    for s in 1:data.statesize
      data.p[s][t1,t2] = x1[s]
    end
    x0 .= x1
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
      θ(t1.val, t2.val) || (val *= ξ[s] * ρ0[s])
      return val
    end
  end

  Δup = dos2gf(dos, grid, β=β)
  Δ = [deepcopy(Δup) for _ in 1:2] # spin symmetric

  data = nca(p0, Δ, params)

  return data
end

function main()
  dos = Keldysh.flat_dos(ν=10.0, D=10.0)
  data = run_anderson_nca(β=1.0, tmax=5.0, npts_real = 101, U = 8.0, dos=dos)
  t, ρt, Zt = populations(data.p)
  h5write("output.h5", "output/rho", ρt)
  h5write("output.h5", "output/Z", ρt)
  h5write("output.h5", "output/t", t)
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
