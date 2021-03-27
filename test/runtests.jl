#!/usr/bin/env julia

using Keldysh, Test, LinearAlgebra

include("branch.jl")
include("contour.jl")
include("time_grid.jl")
include("dos.jl")
include("time_gf.jl")

#@testset "generate_gf" begin
#  let tmax = 1.0, β = 1.0, ν = 1/1000, ϵ = 2.0
#    c = twist(Contour(full_contour, tmax=tmax, β=β))
#    grid = TimeGrid(c, npts_real=51, npts_imag=51)
#    dos = Keldysh.gaussian_dos(ν=ν, ϵ=ϵ)
#
#    hyb1 = dos2gf(dos, β, grid)
#    hyb2 = gf_1level(grid; ϵ, β)
#
#    # gf_1level is gf for a delta function spectrum
#    @test isapprox(hyb1.data, hyb2.data, atol=ν, norm=x -> norm(x, Inf))
#  end
#end
#
