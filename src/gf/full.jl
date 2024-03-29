struct FullTimeGF{T, scalar} <: AbstractTimeGF{T, scalar}
  grid::FullTimeGrid
  gtr::AntiHermitianStorage{T,scalar}
  les::AntiHermitianStorage{T,scalar}
  rm::GenericStorage{T,scalar}
  mat::PeriodicStorage{T,scalar}
  ξ::GFSignEnum
end

norbitals(G::FullTimeGF) = G.gtr.norb

function FullTimeGF(::Type{T}, grid::FullTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) where T <: Number
  nt = grid.nt
  ntau = grid.ntau

  gtr = AntiHermitianStorage(T, nt, norb, scalar)
  les = AntiHermitianStorage(T, nt, norb, scalar)
  rm = GenericStorage(T, ntau, nt, norb, scalar)
  mat = PeriodicStorage(T, ntau, norb, scalar)

  FullTimeGF(grid, gtr, les, rm, mat, ξ)
end
FullTimeGF(grid::FullTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) = FullTimeGF(ComplexF64, grid, norb, ξ, scalar)

@inline function Base.getindex(G::FullTimeGF, k, l, t1::TimeGridPoint, t2::TimeGridPoint, greater=true)
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
    return -ξ * adjoint(G.rm[l,k,ntau+1-j,i]) # Aoki 19c
  else
    return greater ? G.mat[k,l,i,j,greater] : ξ * G.mat[k,l,i,j,greater]
  end
end

function Base.setindex!(G::FullTimeGF, v, k, l, t1::TimeGridPoint, t2::TimeGridPoint)
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
      return G.mat[k,l,i,j] = v
    else
      return G.mat[k,l,i,j] = ξ * v
    end
  end
end

function Base.similar(G::T, ξ::GFSignEnum) where T <: FullTimeGF
  T(G.grid, similar(G.gtr), similar(G.les), similar(G.rm), similar(G.mat), ξ)
end

function Base.zero(G::T, ξ::GFSignEnum) where T <: FullTimeGF
  T(G.grid, zero(G.gtr), zero(G.les), zero(G.rm), zero(G.mat), ξ)
end

function Base.similar(G::T) where T <: FullTimeGF
  similar(G, G.ξ)
end

function Base.zero(G::T) where T <: FullTimeGF
  zero(G, G.ξ)
end

function Base.:+(G1::T, G2::T) where T <: FullTimeGF
  @assert G1.grid == G2.grid
  @assert G1.ξ == G2.ξ
  return T(G1.grid, G1.gtr + G2.gtr, G1.les + G2.les, G1.rm + G2.rm, G1.mat + G2.mat, G1.ξ)
end

function Base.:-(G1::T, G2::T) where T <: FullTimeGF
  @assert G1.grid == G2.grid
  @assert G1.ξ == G2.ξ
  return T(G1.grid, G1.gtr - G2.gtr, G1.les - G2.les, G1.rm - G2.rm, G1.mat - G2.mat, G1.ξ)
end

function Base.:*(G::T, α::Number) where T <: FullTimeGF
  return T(G.grid, G.gtr * α, G.les * α, G.rm * α, G.mat * α, G.ξ)
end

function Base.:*(α::Number, G::T) where T <: FullTimeGF
  return G * α
end

function Base.:-(G::T) where T <: FullTimeGF
  return T(G.grid, -G.gtr, -G.les, -G.rm, -G.mat, G.ξ)
end

function TimeDomain(G::FullTimeGF)
  grid = G.grid
  tplus = grid[forward_branch]
  tminus = reverse(grid[backward_branch])
  tau = grid[imaginary_branch]

  nt = G.grid.nt
  ntau = G.grid.ntau

  points = Tuple{TimeGridPoint, TimeGridPoint}[]

  for i in 1:nt
    for j in 1:i
      push!(points, (tminus[j], tminus[i])) # greater
      push!(points, (tplus[j], tminus[i])) # lesser
    end
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

function FullTimeGF(f::Function, ::Type{T}, grid::FullTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) where T <: Number
  G = FullTimeGF(T, grid, norb, ξ, scalar)
  D = TimeDomain(G)

  for (t1, t2) in D.points
    G[t1, t2] = f(t1, t2)
  end

  return G
end

FullTimeGF(f::Function, grid::FullTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) = FullTimeGF(f, ComplexF64, grid, norb, ξ, scalar)

function FullTimeGF(dos::AbstractDOS, grid::FullTimeGrid)
  # result is time invariant so do calculation for time invariant case and then copy over
  G = TimeInvariantFullTimeGF(dos, grid)
  FullTimeGF(grid,1,fermionic,true) do t1, t2
    G[t1,t2]
  end
end
