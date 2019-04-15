#!/usr/bin/env julia

using Keldysh, Test

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
