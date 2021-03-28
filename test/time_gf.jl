@testset "time_gf" begin

  let c = FullContour(tmax=2.0, β=5.0), nt = 21, ntau=51
    ϵ = -1.0
    f = (t1, t2) -> -1.0im * (heaviside(t1.val, t2.val) - fermi(ϵ, c.β)) * exp(-1.0im * (t1.val.val - t2.val.val) * ϵ)

    grid = FullTimeGrid(c, nt, ntau)

    G1 = GenericTimeGF(f, grid, 1, true)

    @test Keldysh.jump(G1) == 1.0im

    G2 = FullTimeGF(f, grid, 1, fermionic, true)

    @test all(map(((t1, t2),) -> isapprox(G1[t1,t2], G2[t1,t2], atol=1e-10), Iterators.product(grid, grid)))

    G3 = TimeInvariantFullTimeGF(f, grid, 1, fermionic, true)

    @test all(map(((t1, t2),) -> isapprox(G1[t1,t2], G3[t1,t2], atol=1e-10), Iterators.product(grid, grid)))
  end

  let c = KeldyshContour(tmax=2.0), nt = 21

    ϵ = -1.0
    β = 5.0
    f = (t1, t2) -> -1.0im * (heaviside(t1.val, t2.val) - fermi(ϵ, β)) * exp(-1.0im * (t1.val.val - t2.val.val) * ϵ)

    grid = KeldyshTimeGrid(c, nt)

    G1 = GenericTimeGF(f, grid, 1, true)

    @test Keldysh.jump(G1) == 1.0im

    G2 = KeldyshTimeGF(f, grid, 1, fermionic, true)

    @test all(map(((t1, t2),) -> isapprox(G1[t1,t2], G2[t1,t2], atol=1e-10), Iterators.product(grid, grid)))

    G3 = TimeInvariantKeldyshTimeGF(f, grid, 1, fermionic, true)

    @test all(map(((t1, t2),) -> isapprox(G1[t1,t2], G3[t1,t2], atol=1e-10), Iterators.product(grid, grid)))

  end

end