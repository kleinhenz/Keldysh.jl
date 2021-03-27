@testset "time_gf" begin
  let c = FullContour(tmax=2.0, β=5.0), nt = 21, ntau=51
    grid = FullTimeGrid(c, nt, ntau)

    ϵ = -1.0

    G = GenericTimeGF(grid, 1, true) do t1, t2
      -1.0im * (heaviside(t1.val, t2.val) - fermi(ϵ, c.β)) * exp(-1.0im * (t1.val.val - t2.val.val) * ϵ)
    end

    @test Keldysh.jump(G) == 1.0im

    G_ = FullTimeGF(grid, 1, -1, true) do t1, t2
      -1.0im * (heaviside(t1.val, t2.val) - fermi(ϵ, c.β)) * exp(-1.0im * (t1.val.val - t2.val.val) * ϵ)
    end

    @test all(map(((t1, t2),) -> G[t1,t2] ≈ G_[t1,t2], Iterators.product(grid, grid)))

  end
end
