function gf_1level(t1::BranchPoint, t2::BranchPoint; ϵ, β)
    -1.0im * (θ(t1, t2) - fermi(ϵ, β)) * exp(-1.0im * (t1.val - t2.val) * ϵ)
end

function gf_1level(grid::TimeGrid; ϵ, β)
  TimeGF(grid) do t1, t2
    gf_1level(t1.val, t2.val; ϵ, β)
  end
end

function dos2gf(dos, t1::BranchPoint, t2::BranchPoint; β, integrator = dos_integrator)
    theta = θ(t1, t2)
    Δt = t1.val - t2.val
    f = ω -> (ω > 0.0 ? exp(-1.0im * ω * (Δt - 1.0im * (1.0 - theta) * β)) / (exp(-β * ω) + 1) :
                        exp(-1.0im * ω * (Δt + 1.0im * theta * β)) / (exp(β * ω) + 1))
    return -1.0im * (2 * theta - 1) * integrator(f, dos)
end

function _dos2gf_ret(dos, t::Real, integrator=dos_integrator)
  f = ω -> exp(-im * ω * t)
  return -im * integrator(f, dos)
end

function _dos2gf_les(dos, t::Real, β::Real, integrator=dos_integrator)
  f = ω -> exp(-im * ω * t) * fermi(ω, β)
  return im * integrator(f, dos)
end

function _dos2gf_tv(dos, t::Real, τ::Real, β::Real, integrator=dos_integrator)
#   f = ω -> exp(-im * ω * t) * exp(ω * τ) * fermi(ω, β)

  f = ω -> ω > 0.0 ?
    exp(-im * ω * t) * exp(-ω * (β - τ)) * fermi(-ω, β) :
    exp(-im * ω * t) * exp(ω * τ) * fermi(ω, β)

  return im * integrator(f, dos)
end

function _dos2gf_mat(dos, τ::Real, β::Real, integrator=dos_integrator)
#   f = ω -> exp(-ω * τ) * fermi(-ω, β)

  f = ω -> ω > 0.0 ?
    exp(-ω * τ) * fermi(-ω, β) :
    exp(ω * (β - τ)) * fermi(ω, β)

  return -1 * integrator(f, dos)
end

function _dos2gf_full_contour(dos, β, grid::TimeGrid)
  @assert grid.contour.domain == full_contour

  t = realtimes(grid)
  τ = imagtimes(grid)

  les_ = [_dos2gf_les(dos, ti, β) for ti in t]
  les = [i >= j ? les_[i - j + 1] : -conj(les_[j - i + 1]) for i in 1:length(t), j in 1:length(t)]

  ret_ = [_dos2gf_ret(dos, ti) for ti in t]
  ret = [i >= j ? ret_[i - j + 1] : -conj(ret_[j - i + 1]) for i in 1:length(t), j in 1:length(t)]

  tv = [_dos2gf_tv(dos, ti, τj, β) for ti in t, τj in τ]
  mat = [_dos2gf_mat(dos, τi, β) for τi in τ]

  return TimeGF(les, ret, tv, mat, grid)

#   return (;les, ret, tv, mat)
end

function _dos2gf_keldysh_contour(dos, β, grid::TimeGrid)
  @assert grid.contour.domain == keldysh_contour

  t = realtimes(grid)

  les_ = [_dos2gf_les(dos, ti, β) for ti in t]
  les = [i >= j ? les_[i - j + 1] : -conj(les_[j - i + 1]) for i in 1:length(t), j in 1:length(t)]

  ret_ = [_dos2gf_ret(dos, ti) for ti in t]
  ret = [i >= j ? ret_[i - j + 1] : -conj(ret_[j - i + 1]) for i in 1:length(t), j in 1:length(t)]

  return TimeGF(les, ret, grid)
end

function dos2gf(dos, β, grid::TimeGrid)
  if grid.contour.domain == full_contour
    _dos2gf_full_contour(dos, β, grid)
  elseif grid.contour.domain == keldysh_contour
    _dos2gf_keldysh_contour(dos, β, grid)
  else
    error("unsupported")
  end
end
