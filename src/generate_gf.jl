using QuadGK

"""
check/get β from contour or explicitly given value
"""
function get_beta(grid::TimeGrid, β)
  im_b = get_branch(grid.contour, imaginary_branch)
  @assert (im_b === nothing && β !== nothing) || (im_b !== nothing && β === nothing)
  β === nothing ? length(im_b) : β
end

function _make_time_invariant_gf(f, grid)
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

    g.(grid, permutedims(grid))
end

function make_gf(f, grid; time_invariant = false)
  time_invariant && return _make_time_invariant_gf(f, grid)
  return f.(grid, permutedims(grid))
end

function gf_1level(t1::BranchPoint, t2::BranchPoint; ϵ, β)
    -1.0im * (θ(t1, t2) - fermi(ϵ, β)) * exp(-1.0im * (t1.val - t2.val) * ϵ)
end

function gf_1level(grid::TimeGrid; ϵ, β=nothing)
  β = get_beta(grid, β)
  make_gf(grid) do t1, t2
    gf_1level(t1.val, t2.val, ϵ=ϵ, β=β)
  end
end

function dos_integrator(f)
  integral, err = quadgk(f, -Inf, Inf, atol=1e-10, rtol=1e-10, maxevals=10^9, order=21)
  return integral
end

function dos2gf(dos, t1::BranchPoint, t2::BranchPoint; β, integrator=dos_integrator)
#    return integrator(ω -> dos(ω) * (θ(t1, t2) - fermi(ω, β)) * exp(-1.0im * (t1.val - t2.val) * ω))
    theta = θ(t1, t2)
    Δt = t1.val - t2.val
    integrand = ω -> dos(ω) * (ω > 0.0 ? exp(-1.0im * ω * (Δt - 1.0im * (1.0 - theta) * β)) / (exp(-β * ω) + 1) :
                                         exp(-1.0im * ω * (Δt + 1.0im * theta * β)) / (exp(β * ω) + 1))
    return -1.0im * (2 * theta - 1) * integrator(integrand)
end

function dos2gf(dos, grid::TimeGrid; β=nothing, integrator=dos_integrator)
  β = get_beta(grid, β)
  make_gf(grid, time_invariant=true) do t1, t2
    dos2gf(dos, t1.val, t2.val, β=β, integrator=integrator)
  end
end

"""
`flat_dos(;ν=1.0, D=5.0)`

return flat band dos with half-bandwith D and inverse cutoff width ν centered at zero
"""
flat_dos(ω; ν=1.0, D=5.0) = (1.0/π) / ((1 + exp(ν * (ω - D))) * (1 + exp(-ν * (ω + D))))

"""
`gaussian_dos(; ϵ=1.0, ν=1.0)`

return normalized gaussian dos centered at ϵ with width ν
"""
gaussian_dos(ω; ϵ=1.0, ν=1.0) = (1.0 / (2 * sqrt(π * ν))) * exp(-((ω - ϵ)^2)/(4ν))
