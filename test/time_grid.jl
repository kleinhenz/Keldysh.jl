using Keldysh, Test

@testset "time_grid" begin
  let c = FullContour(tmax=2.0, β=5.0), nt = 21, ntau=51
    grid = FullTimeGrid(c, nt, ntau)
    @test_throws MethodError KeldyshTimeGrid(c, nt)
    @test_throws MethodError ImaginaryTimeGrid(c, ntau)

    @test length(grid) == 2*nt + ntau

    x = [grid[i] for i in 1:length(grid)]
    @test length(x) == 2*nt + ntau

    @test length(grid[forward_branch]) == nt
    @test length(grid[backward_branch]) == nt
    @test length(grid[imaginary_branch]) == ntau

    @test step(grid, forward_branch) == 0.1
    @test step(grid, backward_branch) == -0.1
    @test step(grid, imaginary_branch) == -0.1im

    @test integrate(t -> 1, grid) ≈ -1.0im * c.β

    for b in c.branches
      for t in b.(range(0.0, 1.0, length=3))
        tl = Keldysh.find_lower(grid, t)
        @test tl.bpoint.domain == t.domain
        @test tl.cidx != length(grid)
        tu = grid[tl.cidx + 1]

        @test heaviside(c, t, tl.bpoint)
        @test heaviside(c, tu.bpoint, t)
        @test heaviside(c, t, grid[1].bpoint)
        @test heaviside(c, grid[end].bpoint, t)
      end
    end
  end

  let c = KeldyshContour(tmax=2.0), nt = 21, ntau=51
    grid = KeldyshTimeGrid(c, nt)
    @test_throws MethodError FullTimeGrid(c, nt, ntau)
    @test_throws MethodError ImaginaryTimeGrid(c, ntau)

    @test length(grid) == 2*nt

    x = [grid[i] for i in 1:length(grid)]
    @test length(x) == 2*nt

    @test length(grid[forward_branch]) == nt
    @test length(grid[backward_branch]) == nt
    @test_throws AssertionError grid[imaginary_branch]

    @test step(grid, forward_branch) == 0.1
    @test step(grid, backward_branch) == -0.1
    @test_throws AssertionError step(grid, imaginary_branch)

    N = length(grid)
    A = zeros(ComplexF64, N, N)
    B = zeros(ComplexF64, N, N)
    for t1 in grid
      for t2 in grid
        A[t1.cidx, t2.cidx] = integrate(t -> 1.0, grid, t1, t2)
        B[t1.cidx, t2.cidx] = t1.bpoint.val - t2.bpoint.val
      end
    end
    @test maximum(abs.(A .- B)) < 1e-12

  end

  let c = ImaginaryContour(β=5.0), nt = 21, ntau=51
    grid = ImaginaryTimeGrid(c, ntau)
    @test_throws MethodError FullTimeGrid(c, nt, ntau)
    @test_throws MethodError KeldyshTimeGrid(c, nt)

    @test length(grid) == ntau

    x = [grid[i] for i in 1:length(grid)]
    @test length(x) == ntau

    @test length(grid[imaginary_branch]) == ntau
    @test_throws AssertionError grid[forward_branch]
    @test_throws AssertionError grid[backward_branch]

    @test_throws AssertionError step(grid, forward_branch)
    @test_throws AssertionError step(grid, backward_branch)
    @test step(grid, imaginary_branch) == -0.1im
  end

end
