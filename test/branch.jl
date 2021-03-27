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
