abstract type AbstractDOS end
abstract type AbstractDOSIntegrator end

struct DOS{F} <: AbstractDOS
  ωmin::Real
  ωmax::Real
  A::F
end
(dos::DOS)(ω) = dos.A(ω)

"""
Descriptor of an integrable singularity in DOS
"""
struct DOSSingularity
  # Position of the singularity, Ω_p.
  position::Real
  # Asymptotic behavior of the DOS near the singularity, S_p(ω).
  # To be chosen such that \lim_{ω \to Ω_p} [D(ω) - S_p(ω)] = 0
  asymptotics
  # Integral of 'asymptotics' over ω\in[ω_{min};ω_{max}].
  integral
end

"""
  SingularDOS represents densities of states D(ω) according to
  the following decomposition,

  D(ω) = R(ω) + \\sum_{p=1}^P S_p(ω).

  R(ω) is a smooth function on [ω_{min};ω_{max}], and the number P of
  singular contributions S_p(ω) equals the number of integrable singularities
  Ω_p of D(ω). Each term S_p(ω) is regular everywhere on
  ω\\in[ω_{min};ω_{max}] except for its corresponding ω = Ω_p.
"""
struct SingularDOS <: AbstractDOS
  # Support limits
  ωmin::Real
  ωmax::Real
  # Regular part R(ω)
  regular
  # List of singularities
  singularities::Vector{DOSSingularity}
end

SingularDOS(ωmin, ωmax, dos) = SingularDOS(ωmin, ωmax, dos, [])

"""Evaluate singular DOS at a given frequency"""
function (dos::SingularDOS)(ω)
  x = dos.regular(ω)
  for s in dos.singularities
    x += s.asymptotics(ω)
  end
  return x
end

struct DeltaDOS <: AbstractDOS
  ϵ::Vector{Float64}
  w::Vector{Float64}
end

DeltaDOS(ϵ::Number) = DeltaDOS([ϵ])
DeltaDOS(ϵ::Vector) = DeltaDOS(ϵ, ones(length(ϵ)) / length(ϵ))

"""
Support limits of a DOS object
"""
dos_support_limits(dos::AbstractDOS) = (-Inf, Inf)
dos_support_limits(dos::DOS) = (dos.ωmin, dos.ωmax)
dos_support_limits(dos::SingularDOS) = (dos.ωmin, dos.ωmax)

struct GaussKronrodDOSIntegrator <: AbstractDOSIntegrator
  atol::Float64
  rtol::Float64
  maxevals::Int
  order::Int
end

GaussKronrodDOSIntegrator(; atol=1e-10, rtol=1e-10, maxevals=10^9, order=21) = GaussKronrodDOSIntegrator(atol, rtol, maxevals, order)
function (integrator::GaussKronrodDOSIntegrator)(f, a, b)
  integral, err = quadgk(ω -> f(ω),
                         a,
                         b,
                         atol=integrator.atol,
                         rtol=integrator.rtol,
                         maxevals=integrator.maxevals,
                         order=integrator.order)
  return integral
end

function (integrator::GaussKronrodDOSIntegrator)(f, dos::AbstractDOS)
  limits = dos_support_limits(dos)
  return integrator(ω -> f(ω) * dos(ω), limits[1], limits[2])
end

"""
  Integrator for SingularDOS

  \\int_{ω_{min}}^{ω_{max}} dω D(ω) f(ω) =
  \\int_{ω_{min}}^{ω_{max}} dω R(ω) f(ω) +
  \\sum_{p=1}^P \\int_{ω_{min}}^{ω_{max}} dω S_p(ω) [f(ω) - f(Ω_p)] +
  \\sum_{p=1}^P f(Ω_p) \\int_{ω_{min}}^{ω_{max}} dω S_p(ω).

  The integrands in the first two terms of the RHS are smooth,
  and the value of the integral in the last term must be provided in the
  `integral` field of the corresponding `DOSSingularity` structure.
"""
function (integrator::GaussKronrodDOSIntegrator)(f, dos::SingularDOS)
  limits = dos_support_limits(dos)
  val = integrator(ω -> f(ω) * dos.regular(ω), limits[1], limits[2])
  for s in dos.singularities
    f_s = f(s.position)
    val += integrator(ω -> ω ≈ s.position ? .0 : s.asymptotics(ω) * (f(ω) - f_s), limits[1], limits[2])
    val += f_s * s.integral
  end
  return val
end

function (integrator::GaussKronrodDOSIntegrator)(f, dos::DeltaDOS)
  return sum(dos.w .* f.(dos.ϵ))
end

function dos2gf(dos, β, t1::BranchPoint, t2::BranchPoint, integrator = GaussKronrodDOSIntegrator())
    theta = heaviside(t1, t2)
    Δt = t1.val - t2.val
    f = ω -> (ω > 0.0 ? exp(-1.0im * ω * (Δt - 1.0im * (1.0 - theta) * β)) / (exp(-β * ω) + 1) :
                        exp(-1.0im * ω * (Δt + 1.0im * theta * β)) / (exp(β * ω) + 1))
    return -1.0im * (2 * theta - 1) * integrator(f, dos)
end

#
# DOS factory functions
#

"""
`flat_dos(;ν=1.0, D=5.0, μ=0.0)`

return flat band DOS centered at μ with half-bandwith D and inverse cutoff width ν
"""
flat_dos(; ν=1.0, D=5.0, μ=0.0) = DOS(-Inf, Inf, ω -> (1.0/π) * fermi(ν * (ω - μ - D)) * fermi(-ν * (ω - μ + D)))

"""
`gaussian_dos(; ϵ=1.0, ν=1.0)`

return normalized Gaussian DOS centered at ϵ with width ν
"""
gaussian_dos(; ϵ=1.0, ν=1.0) = DOS(-Inf, Inf, ω -> (1.0 / (2 * sqrt(π * ν))) * exp(-((ω - ϵ)^2)/(4ν)))

"""
`bethe_dos(; ϵ=0.0, t=1.0)`

return normalized DOS of a Bethe lattice with hopping constant t centered at ϵ
"""
bethe_dos(; ϵ=0.0, t=1.0) = SingularDOS(-2t + ϵ, 2t + ϵ,
  ω -> begin
    if ω == -2t + ϵ || ω == 2t + ϵ
      -2 / (π*t)
    else
      x = (ω - ϵ) / (2t)
      (sqrt(1 - x*x) - sqrt(2 * (1 - x)) - sqrt(2 * (1 + x))) / (π*t)
    end
  end,
  [
    DOSSingularity(-2t + ϵ, ω -> sqrt(2 * (1 + (ω - ϵ) / (2t))) / (π*t), 16 / (3*π)),
    DOSSingularity( 2t + ϵ, ω -> sqrt(2 * (1 - (ω - ϵ) / (2t))) / (π*t), 16 / (3*π))
  ]
)

"""
`chain_dos(; ϵ=0.0, t=1.0)`

return normalized DOS of a linear chain with hopping constant t centered at ϵ
"""
chain_dos(; ϵ=0.0, t=1.0) = SingularDOS(-2t + ϵ, 2t + ϵ,
  ω -> begin
    if ω == -2t + ϵ || ω == 2*t + ϵ
      -3 / (8*π*t)
    else
      x = (ω - ϵ) / (2t)
      rp = sqrt(2 * (1 - x))
      # The second term regularizes the derivative of near ω = 2t to ease integration
      sp = 1 / (2*π*t * rp) + rp / (16*π*t)

      rm = sqrt(2 * (1 + x))
      # The second term regularizes the derivative near ω = -2t to ease integration
      sm = 1 / (2*π*t * rm) + rm / (16*π*t)

      (1 / (2*π*t)) / sqrt(1 - x^2) - sp - sm
    end
  end,
  [
    DOSSingularity(-2t + ϵ, ω -> begin
                          x = (ω - ϵ) / (2t)
                          s = sqrt(2 * (1 + x))
                          # The second term regularizes the derivative
                          1 / (2*π*t*s) + s / (16*π*t)
                        end,
                    7 / (3*π)),
    DOSSingularity( 2t + ϵ, ω -> begin
                          x = (ω - ϵ) / (2t)
                          s = sqrt(2 * (1 - x))
                          # The second term regularizes the derivative
                          1 / (2*π*t*s) + s / (16*π*t)
                        end,
                    7 / (3*π))
  ]
)

"""
`square_dos(; ϵ=0.0, t=1.0)`

return normalized DOS of a 2D square lattice with hopping constant t centered at ϵ
"""
square_dos(; ϵ=0.0, t=1.0) = SingularDOS(-4t + ϵ, 4t + ϵ,
  ω -> begin
    if ω ≈ ϵ
      0
    else
      x = (ω - ϵ) / (4t)
      (1 / (2 * π^2 * t)) * (Elliptic.K(1 - x^2) + log(abs(x) / 4));
    end
  end,
  [
    DOSSingularity(ϵ, ω -> -1 / (2 * π^2 * t) * log(abs(ω - ϵ) / (16t)),
                   4 / (π^2) * (1 + 2 * log(2)))
  ]
)
