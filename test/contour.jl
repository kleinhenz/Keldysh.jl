using Keldysh, Test

@testset "contour" begin
  let c = FullContour(tmax=2.0, β=5.0)
    @test nbranches(c) == 3

    b = map(x -> x.domain, c.branches)
    @test b == (forward_branch, backward_branch, imaginary_branch)

    c = twist(c)
    b = map(x -> x.domain, c.branches)
    @test b == (backward_branch, imaginary_branch, forward_branch)

    for b in (forward_branch, backward_branch, imaginary_branch)
      @test c[b].domain == b
    end

    @test forward_branch ∈ c
    @test backward_branch ∈ c
    @test imaginary_branch ∈ c
  end

  let c = KeldyshContour(tmax=2.0)
    @test nbranches(c) == 2

    b = map(x -> x.domain, c.branches)
    @test b == (forward_branch, backward_branch)

    c = twist(c)
    b = map(x -> x.domain, c.branches)
    @test b == (backward_branch, forward_branch)

    for b in (forward_branch, backward_branch)
      @test c[b].domain == b
    end

    @test forward_branch ∈ c
    @test backward_branch ∈ c
    @test imaginary_branch ∉ c
  end

  let c = ImaginaryContour(β=5.0)
    @test nbranches(c) == 1

    b = map(x -> x.domain, c.branches)
    @test b == (imaginary_branch,)

    for b in (imaginary_branch,)
      @test c[b].domain == b
    end

    @test forward_branch ∉ c
    @test backward_branch ∉ c
    @test imaginary_branch ∈ c
  end

  let c = twist(FullContour(tmax=1.0, β=5.0))
    for ref in [0.0, 0.5, 2.0, 5.0, 5.5, 6.5]
        @test get_ref(c, get_point(c, ref)) == ref
        @test get_ref(c, c(ref)) == ref
    end
  end
end