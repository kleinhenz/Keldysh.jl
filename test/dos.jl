@testset "dos" begin
  # Flat DOS
  let ν=5.0, D=2.0, dos = Keldysh.flat_dos(ν=ν, D=D)
    α = ν*D
    moments = [dos_integrator(ω -> ω^n, dos) for n = 0:3]
    moments_ref = (1/π)*[D*(1 + coth(α)),
                         0.,
                         D^3*(π^2 + α^2)*(1 + coth(α))/(3α^2),
                         0.]
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
  let t=2.0, dos = Keldysh.bethe_dos(t=t)
    moments = [dos_integrator(ω -> ω^n, dos) for n = 0:6]
    moments_ref = [1.0, 0., t^2, 0., 2t^4, 0., 5t^6]
    @test isapprox(moments, moments_ref, atol=1e-10, rtol=1e-10)
  end

  # Linear chain DOS
  let t=2.0, dos = Keldysh.chain_dos(t=t)
    moments = [dos_integrator(ω -> ω^n, dos) for n = 0:6]
    moments_ref = [1.0, 0., 2t^2, 0., 6t^4, 0., 20t^6]
    @test isapprox(moments, moments_ref, atol=1e-10, rtol=1e-10)
  end

  # Square lattice DOS
  let t=2.0, dos = Keldysh.square_dos(t=t)
    moments = [dos_integrator(ω -> ω^n, dos) for n = 0:6]
    moments_ref = [1.0, 0., 4t^2, 0., 36t^4, 0., 400t^6]
    @test isapprox(moments, moments_ref, atol=1e-10, rtol=1e-10)
  end
end
