struct TimeInvariantFullTimeGF{T, scalar} <: AbstractTimeGF{T, scalar}
  grid::FullTimeGrid
  gtr::AntiHermitianToeplitzStorage{T,scalar}
  les::AntiHermitianToeplitzStorage{T,scalar}
  rm::GenericStorage{T,scalar}
  mat::PeriodicStorage{T,scalar}
  ξ::GFSignEnum
end

norbitals(G::TimeInvariantFullTimeGF) = G.gtr.norb

function TimeInvariantFullTimeGF(::Type{T}, grid::FullTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) where T <: Number
  nt = grid.nt
  ntau = grid.ntau

  gtr = AntiHermitianToeplitzStorage(T, nt, norb, scalar)
  les = AntiHermitianToeplitzStorage(T, nt, norb, scalar)
  rm = GenericStorage(T, ntau, nt, norb, scalar)
  mat = PeriodicStorage(T, ntau, norb, scalar)

  TimeInvariantFullTimeGF(grid, gtr, les, rm, mat, ξ)
end
TimeInvariantFullTimeGF(grid::FullTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) = TimeInvariantFullTimeGF(ComplexF64, grid, norb, ξ, scalar)

@inline function Base.getindex(G::TimeInvariantFullTimeGF, k, l, t1::TimeGridPoint, t2::TimeGridPoint, greater=true)
  greater = t1 == t2 ? greater : heaviside(t1.bpoint, t2.bpoint)

  i = t1.ridx
  j = t2.ridx
  ξ = Int(G.ξ)

  real_branches = (forward_branch, backward_branch)
  if (t1.bpoint.domain ∈ real_branches && t2.bpoint.domain ∈ real_branches)
    return greater ? G.gtr[k,l,i,j] : G.les[k,l,i,j]
  elseif (t1.bpoint.domain == imaginary_branch && (t2.bpoint.domain == forward_branch || t2.bpoint.domain == backward_branch))
    return G.rm[k,l,i,j]
  elseif ((t1.bpoint.domain == forward_branch || t1.bpoint.domain == backward_branch) && t2.bpoint.domain == imaginary_branch)
    ntau = G.grid.ntau
    return copy(-ξ * adjoint(G.rm[l,k,ntau+1-j,i])) # Aoki 19c
  else
    greater ? G.mat[k,l,i,j,greater] : ξ * G.mat[k,l,i,j,greater]
  end
end

function Base.setindex!(G::TimeInvariantFullTimeGF, v, k, l, t1::TimeGridPoint, t2::TimeGridPoint)
  greater = heaviside(t1.bpoint, t2.bpoint)

  i = t1.ridx
  j = t2.ridx
  ξ = Int(G.ξ)

  real_branches = (forward_branch, backward_branch)
  if (t1.bpoint.domain ∈ real_branches && t2.bpoint.domain ∈ real_branches)
    return greater ? G.gtr[k,l,i,j] = v : G.les[k,l,i,j] = v
  elseif (t1.bpoint.domain == imaginary_branch && (t2.bpoint.domain == forward_branch || t2.bpoint.domain == backward_branch))
    return G.rm[k,l,i,j] = v
  elseif ((t1.bpoint.domain == forward_branch || t1.bpoint.domain == backward_branch) && t2.bpoint.domain == imaginary_branch)
    ntau = G.grid.ntau
    return G.rm[l,k,ntau+1-j,i] = -ξ * adjoint(v) # Aoki 19c
  else
    if greater
      G.mat[k,l,i,j] = v
    else
      G.mat[k,l,i,j] = ξ * v
    end
  end
end

function Base.similar(G::T, ξ::GFSignEnum) where T <: TimeInvariantFullTimeGF
  T(G.grid, similar(G.gtr), similar(G.les), similar(G.rm), similar(G.mat), ξ)
end

function Base.zero(G::T, ξ::GFSignEnum) where T <: TimeInvariantFullTimeGF
  T(G.grid, zero(G.gtr), zero(G.les), zero(G.rm), zero(G.mat), ξ)
end

function Base.similar(G::T) where T <: TimeInvariantFullTimeGF
  similar(G, G.ξ)
end

function Base.zero(G::T) where T <: TimeInvariantFullTimeGF
  zero(G, G.ξ)
end

function TimeDomain(G::TimeInvariantFullTimeGF)
  grid = G.grid
  tplus = grid[forward_branch]
  tminus = reverse(grid[backward_branch])
  tau = grid[imaginary_branch]

  nt = length(tplus)
  ntau = length(tau)

  points = Tuple{TimeGridPoint, TimeGridPoint}[]

  for i in 1:nt
    push!(points, (tminus[1], tminus[i])) # greater
    push!(points, (tplus[1], tminus[i])) # lesser
  end

  for i in 1:nt
    for j in 1:ntau
      push!(points, (tau[j], tminus[i])) # right mixing
    end
  end

  for i in 1:ntau
    push!(points, (tau[i], tau[1])) # matsubara
  end

  return TimeDomain(points)
end

function TimeInvariantFullTimeGF(f::Function, ::Type{T}, grid::FullTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) where T <: Number
  G = TimeInvariantFullTimeGF(T, grid, norb, ξ, scalar)
  D = TimeDomain(G)

  for (t1, t2) in D.points
    G[t1, t2] = f(t1, t2)
  end

  return G
end

TimeInvariantFullTimeGF(f::Function, grid::FullTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) = TimeInvariantFullTimeGF(f, ComplexF64, grid, norb, ξ, scalar)

function TimeInvariantFullTimeGF(dos::AbstractDOS, grid::FullTimeGrid)
  β = grid.contour.β
  TimeInvariantFullTimeGF(grid, 1, fermionic, true) do t1, t2
    Keldysh.dos2gf(dos, β, t1.bpoint, t2.bpoint)
  end
end

function interpolate(G::TimeInvariantFullTimeGF{T, true}, t1::BranchPoint, t2::BranchPoint) where T
  if t1.domain == imaginary_branch && t2.domain == imaginary_branch
    return interpolate_imag(G, t1, t2)
  else
    return interpolate_generic(G, t1, t2)
  end
end

function interpolate!(x, G::TimeInvariantFullTimeGF{T, false}, t1::BranchPoint, t2::BranchPoint) where T
  if t1.domain == imaginary_branch && t2.domain == imaginary_branch
    return interpolate_imag!(x, G, t1, t2)
  else
    return interpolate_generic!(x, G, t1, t2)
  end
end
