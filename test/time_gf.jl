using Keldysh, Test

@testset "time_gf" begin

  function test_arithmetics(G::AbstractTimeGF, f::Function, grid::AbstractTimeGrid)
    tt_range = Iterators.product(grid, grid)

    A = G + G
    @test all(map(((t1, t2),) -> isapprox(A[t1,t2], 2*f(t1, t2), atol=1e-10), tt_range))
    B = A - G
    @test all(map(((t1, t2),) -> isapprox(B[t1,t2], f(t1, t2), atol=1e-10), tt_range))
    C = 3 * G
    @test all(map(((t1, t2),) -> isapprox(C[t1,t2], 3*f(t1, t2), atol=1e-10), tt_range))
    D = G * 3
    @test all(map(((t1, t2),) -> isapprox(D[t1,t2], 3*f(t1, t2), atol=1e-10), tt_range))
    E = -G
    @test all(map(((t1, t2),) -> isapprox(E[t1,t2], -f(t1, t2), atol=1e-10), tt_range))
  end

  @testset "FullContour" begin
    let c = FullContour(tmax=2.0, β=5.0), nt = 21, ntau=51

      ϵ = -1.0
      f = (t1, t2) -> -1.0im * (heaviside(t1.bpoint, t2.bpoint) - fermi(ϵ, c.β)) *
        exp(-1.0im * (t1.bpoint.val - t2.bpoint.val) * ϵ)

      grid = FullTimeGrid(c, nt, ntau)

      G1 = GenericTimeGF(f, grid, 1, true)
      G1_ = similar(G1)
      G1_ = zero(G1)
      test_arithmetics(G1, f, grid)

      @test Keldysh.jump(G1) == 1.0im

      G2 = FullTimeGF(f, grid, 1, fermionic, true)
      G2_ = similar(G2)
      G2_ = zero(G2)
      test_arithmetics(G2, f, grid)

      @test all(map(((t1, t2),) -> isapprox(G1[t1,t2], G2[t1,t2], atol=1e-10),
                    Iterators.product(grid, grid)))

      G3 = TimeInvariantFullTimeGF(f, grid, 1, fermionic, true)
      G3_ = similar(G3)
      G3_ = zero(G3)
      test_arithmetics(G3, f, grid)

      @test all(map(((t1, t2),) -> isapprox(G1[t1,t2], G3[t1,t2], atol=1e-10),
                    Iterators.product(grid, grid)))

      # test interpolation on linear function
      lin_f = (t1, t2) -> (2.0 * t1.val - im * t2.val) + 5.0 * heaviside(t1, t2)
      for t1 in grid
        for t2 in grid
          G1[t1, t2] = lin_f(t1.bpoint, t2.bpoint)
        end
      end

      grid_fine = FullTimeGrid(c, 41, 101)

      @test all(map(Iterators.product(grid_fine, grid_fine)) do (t1, t2)
                  isapprox(G1(t1.bpoint, t2.bpoint), lin_f(t1.bpoint, t2.bpoint), atol=1e-10)
                end
               )
    end
  end

  @testset "KeldyshContour" begin
    let c = KeldyshContour(tmax=2.0), nt = 21

      ϵ = -1.0
      β = 5.0
      f = (t1, t2) -> -1.0im * (heaviside(t1.bpoint, t2.bpoint) - fermi(ϵ, β)) *
        exp(-1.0im * (t1.bpoint.val - t2.bpoint.val) * ϵ)

      grid = KeldyshTimeGrid(c, nt)

      G1 = GenericTimeGF(f, grid, 1, true)
      G1_ = similar(G1)
      G1_ = zero(G1)
      test_arithmetics(G1, f, grid)

      @test Keldysh.jump(G1) == 1.0im

      G2 = KeldyshTimeGF(f, grid, 1, fermionic, true)
      G2_ = similar(G2)
      G2_ = zero(G2)
      test_arithmetics(G2, f, grid)

      @test all(map(((t1, t2),) -> isapprox(G1[t1,t2], G2[t1,t2], atol=1e-10),
                    Iterators.product(grid, grid)))

      G3 = TimeInvariantKeldyshTimeGF(f, grid, 1, fermionic, true)
      G3_ = similar(G3)
      G3_ = zero(G3)
      test_arithmetics(G3, f, grid)

      @test all(map(((t1, t2),) -> isapprox(G1[t1,t2], G3[t1,t2], atol=1e-10),
                    Iterators.product(grid, grid)))

    end
  end

  @testset "ImaginaryContour" begin
    let c = ImaginaryContour(β=5.0), ntau=51

      ϵ = -1.0
      f = (t1, t2) -> -1.0im * (heaviside(t1.bpoint, t2.bpoint) - fermi(ϵ, c.β)) *
        exp(-1.0im * (t1.bpoint.val - t2.bpoint.val) * ϵ)

      grid = ImaginaryTimeGrid(c, ntau)

      G1 = GenericTimeGF(f, grid, 1, true)
      G1_ = similar(G1)
      G1_ = zero(G1)
      test_arithmetics(G1, f, grid)

      G2 = ImaginaryTimeGF(f, grid, 1, fermionic, true)
      G2_ = similar(G2)
      G2_ = zero(G2)
      test_arithmetics(G2, f, grid)

      @test all(map(((t1, t2),) -> isapprox(G1[t1,t2], G2[t1,t2], atol=1e-10),
                    Iterators.product(grid, grid)))

    end
  end

end
