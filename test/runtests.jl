#!/usr/bin/env julia

using Keldysh, Test, LinearAlgebra

@testset "branch" begin
  tmax = 2.0
  beta = 3.0

  fwd = Branch(forward_branch, tmax)
  back = Branch(backward_branch, tmax)
  imag = Branch(imaginary_branch, beta)

  @test fwd.domain == forward_branch
  @test back.domain == backward_branch
  @test imag.domain == imaginary_branch

  @test length(fwd) == tmax
  @test length(back) == tmax
  @test length(imag) == beta

end

@testset "contour" begin
  c = Contour(full_contour, tmax=2.0, β=5.0)
  for i in 1:3
    @test c.branches[i].domain == BranchEnum(i)
  end

  twist!(c)
  for i in 1:3
    @test c.branches[i].domain == BranchEnum(mod1(i+1, 3))
  end
end

@testset "time_grid" begin
  c = Contour(full_contour, tmax=2.0, β=5.0)
  grid = TimeGrid(c, npts_real=21, npts_imag=51)

  @test grid.step[1] ≈ 0.1
  @test grid.step[2] ≈ -0.1
  @test grid.step[3] ≈ -0.1im
end

@testset "generate_gf" begin
  let tmax = 1.0, β = 1.0, D=10.0, ν = 10.0

    c = twist!(Contour(full_contour, tmax=tmax, β=β))
    grid = TimeGrid(c, npts_real=51, npts_imag=51)
    dos = ω -> (1.0/π) / ((1 + exp(ν * (ω - D))) * (1 + exp(-ν * (ω + D))))

    @time hyb1 = make_gf(grid, time_invariant=false) do t1, t2
      dos2gf(t1, t2, dos=dos, β=β)
    end

    @time hyb2 = make_gf(grid, time_invariant=true) do t1, t2
      dos2gf(t1, t2, dos=dos, β=β)
    end

    @test hyb1 ≈ hyb2
  end

  let tmax = 1.0, β = 1.0, ν = 1/1000, ϵ = 2.0
    c = twist!(Contour(full_contour, tmax=tmax, β=β))
    grid = TimeGrid(c, npts_real=51, npts_imag=51)

    dos = ω -> (1.0 / (2 * sqrt(π * ν))) * exp(-((ω - ϵ)^2)/(4ν)) # ~δ(ω - ϵ)
    hyb1 = dos2gf(grid, dos=dos)
    hyb2 = gf_1level(grid, ϵ=ϵ)

    # gf_1level is gf for a delta function spectrum
    @test isapprox(hyb1, hyb2, atol=ν, norm=x -> norm(x, Inf))
  end
end
