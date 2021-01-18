#!/usr/bin/env julia

using Keldysh, Test, LinearAlgebra

@testset "branch" begin
  let tmax = 2.0, β = 3.0
    fwd = Branch(forward_branch, tmax)
    back = Branch(backward_branch, tmax)
    imag = Branch(imaginary_branch, β)

    @test fwd.domain == forward_branch
    @test back.domain == backward_branch
    @test imag.domain == imaginary_branch

    @test length(fwd) == tmax
    @test length(back) == tmax
    @test length(imag) == β

    @test fwd(0.0).val == 0.0
    @test fwd(1.0).val == tmax

    @test back(0.0).val == tmax
    @test back(1.0).val == 0.0

    @test imag(0.0).val == 0.0
    @test imag(1.0).val == -1.0im * β
  end
end

@testset "contour" begin
  let c = Contour(full_contour, tmax=2.0, β=5.0)
    for i in 1:3
      @test c.branches[i].domain == BranchEnum(i)
    end

    c = twist(c)
    for i in 1:3
      @test c.branches[i].domain == BranchEnum(mod1(i+1, 3))
    end

    @test nbranches(full_contour) == 3
    @test nbranches(keldysh_contour) == 2
    @test nbranches(imaginary_contour) == 1

    for b in instances(BranchEnum)
      @test get_branch(c, b).domain == b
    end
  end

  let tmax = 2.0, β = 5.0, npts_real = 21
    c = Contour(keldysh_contour, tmax=tmax)
    @test forward_branch ∈ c
    @test backward_branch ∈ c
    @test imaginary_branch ∉ c
  end
end

@testset "time_grid" begin
  let tmax = 2.0, β = 5.0, npts_real=21, npts_imag=51
    c = Contour(full_contour, tmax=tmax, β=β)
    grid = TimeGrid(c, npts_real=npts_real, npts_imag=npts_imag)
    @test grid.step[1] ≈ 0.1
    @test grid.step[2] ≈ -0.1
    @test grid.step[3] ≈ -0.1im

    @test length(grid, forward_branch) == npts_real
    @test length(grid, backward_branch) == npts_real
    @test length(grid, imaginary_branch) == npts_imag

    @test map(p -> p.idx, grid) == 1:(2npts_real + npts_imag)

    for i in 1:3
      @test grid.branch_bounds[i][1].val == grid.contour.branches[i](0.0)
      @test grid.branch_bounds[i][2].val == grid.contour.branches[i](1.0)
    end

    @test integrate(t -> 1, grid) ≈ -1.0im * β
  end

  let tmax = 2.0, β = 5.0, npts_real = 21
    c = Contour(keldysh_contour, tmax=tmax)
    grid = TimeGrid(c, npts_real=npts_real)
    @test grid.step[1] ≈ 0.1
    @test grid.step[2] ≈ -0.1
    @test length(grid, forward_branch) == npts_real
    @test length(grid, backward_branch) == npts_real

    @test map(p -> p.idx, grid) == 1:(2npts_real)

    for i in 1:2
      @test grid.branch_bounds[i][1].val == grid.contour.branches[i](0.0)
      @test grid.branch_bounds[i][2].val == grid.contour.branches[i](1.0)
    end

  end

  let tmax = 2.0, c = Contour(keldysh_contour, tmax=tmax)
    grid = TimeGrid(c, npts_real=51)

    Δt1 = TimeGF(grid) do t1, t2
      t1.val.val - t2.val.val
    end

    Δt2 = TimeGF(grid) do t1, t2
      integrate(t -> 1.0, grid, t1, t2)
    end
    @test Δt1 ≈ Δt2
  end

  let tmax = 2.0, c = twist(Contour(keldysh_contour, tmax=tmax))
    grid = TimeGrid(c, npts_real=51)

    Δt1 = TimeGF(grid) do t1, t2
      t1.val.val - t2.val.val
    end

    Δt2 = TimeGF(grid) do t1, t2
      integrate(t -> 1.0, grid, t1, t2)
    end
    @test Δt1 ≈ Δt2

    A = TimeGF(grid, Norb=2) do t1, t2
      I * (t1.val.val - t2.val.val)        
    end
    B = TimeGF(grid, Norb=2) do t1, t2
      I * integrate(t -> 1.0, grid, t1, t2)
    end
    @test A ≈ B
      
  end

end

@testset "generate_gf" begin
  let tmax = 1.0, β = 1.0, ν = 1/1000, ϵ = 2.0
    c = twist(Contour(full_contour, tmax=tmax, β=β))
    grid = TimeGrid(c, npts_real=51, npts_imag=51)
    dos = Keldysh.gaussian_dos(ν=ν, ϵ=ϵ)

    hyb1 = dos2gf(dos, β, grid)
    hyb2 = gf_1level(grid; ϵ, β)

    # gf_1level is gf for a delta function spectrum
    @test isapprox(hyb1, hyb2, atol=ν, norm=x -> norm(x, Inf))
  end
end

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
