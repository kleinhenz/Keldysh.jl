
#function gf_1level(t1::BranchPoint, t2::BranchPoint; ϵ, β)
#    -1.0im * (heaviside(t1, t2) - fermi(ϵ, β)) * exp(-1.0im * (t1.val - t2.val) * ϵ)
#end
#
#function gf_1level(grid::TimeGrid; ϵ, β)
#  TimeGF(grid, 1, true) do t1, t2
#    gf_1level(t1.val, t2.val; ϵ, β)
#  end
#end

#function _dos2gf_ret(dos, t::Real, integrator=dos_integrator)
#  f = ω -> exp(-im * ω * t)
#  return -im * integrator(f, dos)
#end
#
#function _dos2gf_les(dos, t::Real, β::Real, integrator=dos_integrator)
#  f = ω -> exp(-im * ω * t) * fermi(ω, β)
#  return im * integrator(f, dos)
#end
#
#function _dos2gf_tv(dos, t::Real, τ::Real, β::Real, integrator=dos_integrator)
##   f = ω -> exp(-im * ω * t) * exp(ω * τ) * fermi(ω, β)
#
#  f = ω -> ω > 0.0 ?
#    exp(-im * ω * t) * exp(-ω * (β - τ)) * fermi(-ω, β) :
#    exp(-im * ω * t) * exp(ω * τ) * fermi(ω, β)
#
#  return im * integrator(f, dos)
#end
#
#function _dos2gf_mat(dos, τ::Real, β::Real, integrator=dos_integrator)
##   f = ω -> exp(-ω * τ) * fermi(-ω, β)
#
#  f = ω -> ω > 0.0 ?
#    exp(-ω * τ) * fermi(-ω, β) :
#    exp(ω * (β - τ)) * fermi(ω, β)
#
#  return -1 * integrator(f, dos)
#end
#
#function _dos2gf_full_contour(dos, β, grid::TimeGrid)
#  @assert grid.contour.domain == full_contour
#
#  t = realtimes(grid)
#  τ = imagtimes(grid)
#
#  nt = length(t)
#  ntau = length(τ)
#
#  les_ = [_dos2gf_les(dos, ti, β) for ti in t]
#  les = [i >= j ? les_[i - j + 1] : -conj(les_[j - i + 1]) for i in 1:length(t), j in 1:length(t)]
#  les = reshape(les, 1, 1, nt, nt)
#
#  ret_ = [_dos2gf_ret(dos, ti) for ti in t]
#  ret = [i >= j ? ret_[i - j + 1] : -conj(ret_[j - i + 1]) for i in 1:length(t), j in 1:length(t)]
#  ret = reshape(ret, 1, 1, nt, nt)
#
#  tv = [_dos2gf_tv(dos, ti, τj, β) for ti in t, τj in τ]
#  tv = reshape(tv, 1, 1, nt, ntau)
#
#  mat = [_dos2gf_mat(dos, τi, β) for τi in τ]
#  mat = reshape(mat, 1, 1, ntau)
#
#  return TimeGF(les, ret, tv, mat, grid, true)
#end
#
#function _dos2gf_keldysh_contour(dos, β, grid::TimeGrid)
#  @assert grid.contour.domain == keldysh_contour
#
#  t = realtimes(grid)
#  nt = length(t)
#
#  les_ = [_dos2gf_les(dos, ti, β) for ti in t]
#  les = [i >= j ? les_[i - j + 1] : -conj(les_[j - i + 1]) for i in 1:length(t), j in 1:length(t)]
#  les = reshape(les, 1, 1, nt, nt)
#
#  ret_ = [_dos2gf_ret(dos, ti) for ti in t]
#  ret = [i >= j ? ret_[i - j + 1] : -conj(ret_[j - i + 1]) for i in 1:length(t), j in 1:length(t)]
#  ret = reshape(ret, 1, 1, nt, nt)
#
#  return TimeGF(les, ret, grid, true)
#end
#
#function dos2gf(dos, β, grid::TimeGrid)
#  if grid.contour.domain == full_contour
#    _dos2gf_full_contour(dos, β, grid)
#  elseif grid.contour.domain == keldysh_contour
#    _dos2gf_keldysh_contour(dos, β, grid)
#  else
#    error("unsupported")
#  end
#end
