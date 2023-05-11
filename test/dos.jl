using Keldysh, Test

@testset "dos" begin
  dos_integrator = Keldysh.GaussKronrodDOSIntegrator()
  # Flat DOS
  let ν=5.0, D=2.0, μ=0.5, dos = Keldysh.flat_dos(ν=ν, D=D, μ=μ)
    α = ν*D
    moments = [dos_integrator(ω -> ω^n, dos) for n = 0:3]
    σ2 = D^2*(π^2 + α^2)/(3α^2)
    moments_ref = (D/π) * (1 + coth(α)) * [1.0, μ, μ^2 + σ2, μ^3 + 3μ*σ2]
    @test isapprox(moments, moments_ref, atol=1e-10, rtol=1e-10)
  end

  # Gaussian DOS
  let ν = 2.0, ϵ = 1.0, dos = Keldysh.gaussian_dos(ν=ν, ϵ=ϵ)
    moments = [dos_integrator(ω -> ω^n, dos) for n = 0:6]
    moments_ref = [1.0,
                   ϵ,
                   ϵ^2 + 2ν,
                   ϵ^3 + 6ϵ*ν,
                   ϵ^4 + 12ϵ^2*ν + 12*ν^2,
                   ϵ^5 + 20ϵ^3*ν + 60ϵ*ν^2,
                   ϵ^6 + 30ϵ^4*ν + 180ϵ^2*ν^2 + 120ν^3
                  ]
    @test isapprox(moments, moments_ref, atol=1e-10, rtol=1e-10)
  end

  # Bethe DOS
  let t=2.0, ϵ=0.5, dos = Keldysh.bethe_dos(t=t, ϵ=ϵ)
    moments = [dos_integrator(ω -> ω^n, dos) for n = 0:6]
    moments_ref = [1.0,
                   ϵ,
                   ϵ^2 + t^2,
                   ϵ^3 + 3ϵ*t^2,
                   ϵ^4 + 6ϵ^2 * t^2 + 2t^4,
                   ϵ^5 + 10ϵ^3*t^2 + 10ϵ*t^4,
                   ϵ^6 + 15ϵ^4*t^2 + 30ϵ^2*t^4+ 5t^6]
    @test isapprox(moments, moments_ref, atol=1e-10, rtol=1e-10)
  end

  # Linear chain DOS
  let t=2.0, ϵ=0.5, dos = Keldysh.chain_dos(t=t, ϵ=ϵ)
    moments = [dos_integrator(ω -> ω^n, dos) for n = 0:6]
    moments_ref = [1.0,
                   ϵ,
                   ϵ^2 + 2t^2,
                   ϵ^3 + 6ϵ*t^2,
                   ϵ^4 + 12ϵ^2*t^2 + 6t^4,
                   ϵ^5 + 20ϵ^3*t^2 + 30ϵ*t^4,
                   ϵ^6 + 30ϵ^4*t^2 + 90ϵ^2*t^4 + 20t^6]
    @test isapprox(moments, moments_ref, atol=1e-10, rtol=1e-10)
  end

  # Square lattice DOS
  let t=2.0, ϵ=0.5, dos = Keldysh.square_dos(t=t, ϵ=ϵ)
    moments = [dos_integrator(ω -> ω^n, dos) for n = 0:6]
    moments_ref = [1.0,
                   ϵ,
                   ϵ^2 + 4t^2,
                   ϵ^3 + 12ϵ*t^2,
                   ϵ^4 + 24ϵ^2*t^2 + 36t^4,
                   ϵ^5 + 40ϵ^3*t^2 + 180ϵ*t^4,
                   ϵ^6 + 60ϵ^4*t^2 + 540ϵ^2*t^4 + 400t^6]
    @test isapprox(moments, moments_ref, atol=1e-10, rtol=1e-10)
  end
end
