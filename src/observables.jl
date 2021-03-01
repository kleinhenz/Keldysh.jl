"""
    density(G::TimeGF)

n(j,t) = <c†(j,t)c(j,t)> = -i G<(j,j,t,t)
"""
function density(G::TimeGF)
  les = G[:lesser]

  norb = size(les,1)
  nt = size(les,3) - 1

  n = zeros(ComplexF64, norb, nt+1)

  for j in 1:(nt+1)
    for i in 1:norb
      n[i,j] = -im * les[i,i,j,j]
    end
  end

  return n
end

#"""
#    equilibrium_spectrum(gf::TimeGF, ω)
#
#A(ω) = Im -1/π ∫dt Gʳ(t, 0) exp(iωt)
#"""
#function equilibrium_spectrum(gf::TimeGF, ω)
#  grid = gf.grid
#  t = realtimes(grid)
#  trapz = ft -> sum((ft[i] + ft[i+1]) * (t[i+1] - t[i]) / 2 for i in 1:(length(t) - 1))
#  gfʳ = gf[:retarded][:,1]
#  Aω = [imag((-1.0 / π) * trapz(gfʳ .* exp.(1.0im .* ωi .* t))) for ωi in ω]
#end
#
#function aux_spectrum(gf::TimeGF, ω::Number)
#  grid = gf.grid
#
#  delta_f = (t1, t2) -> -1.0im *(θ(t1.val, t2.val) - 1) * exp(-1.0im * (t1.val.val - t2.val.val) * ω)
#  delta_e = (t1, t2) -> -1.0im *(θ(t1.val, t2.val) - 0) * exp(-1.0im * (t1.val.val - t2.val.val) * ω)
#
#  At = map(zip(grid[forward_branch], reverse(grid[backward_branch]))) do (t⁺, t⁻)
#    Iω_f = -2.0 * integrate(t -> @inbounds(gf[t⁺, t] * delta_f(t, t⁻)), grid, t⁺, t⁻)
#    Iω_e = 2.0 * integrate(t -> @inbounds(delta_e(t⁺, t) * gf[t, t⁻]), grid, t⁺, t⁻)
#    (1 / (2π)) * (Iω_e - Iω_f)
#  end
#end
#
#"""
#    aux_spectrum(gf::TimeGF, ω)
#
#Compute auxiliary current spectral function of gf for frequencies ω
#See Eq. 24 from [1]
#
#[1] Cohen, Guy, David R. Reichman, Andrew J. Millis, and Emanuel Gull. “Green’s
#Functions from Real-Time Bold-Line Monte Carlo.” Physical Review B 89, no. 11
#(March 31, 2014): 115139. https://doi.org/10.1103/PhysRevB.89.115139.
#"""
#function aux_spectrum(gf::TimeGF, ω)
#  reduce(vcat, (aux_spectrum(gf, ωi)' for ωi in ω))
#end
