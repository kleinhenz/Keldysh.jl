"""
check/get β from contour or explicitly given value
"""
function get_beta(grid::TimeGrid, β)
  im_b = get_branch(grid.contour, imaginary_branch)
  @assert (im_b === nothing && β !== nothing) || (im_b !== nothing && β === nothing)
  β === nothing ? length(im_b) : β
end

function gf_1level(t1::BranchPoint, t2::BranchPoint; ϵ, β)
    -1.0im * (θ(t1, t2) - fermi(ϵ, β)) * exp(-1.0im * (t1.val - t2.val) * ϵ)
end

function gf_1level(grid::TimeGrid; ϵ, β=nothing)
  β = get_beta(grid, β)
  TimeGF(grid) do t1, t2
    gf_1level(t1.val, t2.val, ϵ=ϵ, β=β)
  end
end

function dos2gf(dos, t1::BranchPoint, t2::BranchPoint; β, integrator = dos_integrator)
    theta = θ(t1, t2)
    Δt = t1.val - t2.val
    f = ω -> (ω > 0.0 ? exp(-1.0im * ω * (Δt - 1.0im * (1.0 - theta) * β)) / (exp(-β * ω) + 1) :
                        exp(-1.0im * ω * (Δt + 1.0im * theta * β)) / (exp(β * ω) + 1))
    return -1.0im * (2 * theta - 1) * integrator(f, dos)
end

function dos2gf(dos, grid::TimeGrid;
                β=nothing, integrator=dos_integrator,
                ωmin = -Inf, ωmax = Inf, singularities = Real[])
  β = get_beta(grid, β)
  TimeGF(grid) do t1, t2
    dos2gf(dos, t1.val, t2.val, β=β, integrator=integrator)
  end
end
