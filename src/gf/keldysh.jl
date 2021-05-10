struct KeldyshTimeGF{T, scalar} <: AbstractTimeGF{T, scalar}
  grid::KeldyshTimeGrid
  gtr::AntiHermitianStorage{T,scalar}
  les::AntiHermitianStorage{T,scalar}
  ξ::GFSignEnum
end

norbitals(G::KeldyshTimeGF) = G.gtr.norb

function KeldyshTimeGF(::Type{T}, grid::KeldyshTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) where T <: Number
  nt = grid.nt

  gtr = AntiHermitianStorage(T, nt, norb, scalar)
  les = AntiHermitianStorage(T, nt, norb, scalar)

  KeldyshTimeGF(grid, gtr, les, ξ)
end
KeldyshTimeGF(grid::KeldyshTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) = KeldyshTimeGF(ComplexF64, grid, norb, ξ, scalar)

@inline function Base.getindex(G::KeldyshTimeGF, t1::TimeGridPoint, t2::TimeGridPoint, greater=true)
  greater = t1 == t2 ? greater : heaviside(t1.bpoint, t2.bpoint)

  i = t1.ridx
  j = t2.ridx

  return greater ? G.gtr[i,j] : G.les[i,j]
end

function Base.setindex!(G::KeldyshTimeGF, v, t1::TimeGridPoint, t2::TimeGridPoint)
  greater = heaviside(t1.bpoint, t2.bpoint)

  i = t1.ridx
  j = t2.ridx

  return greater ? G.gtr[i,j] = v : G.les[i,j] = v
end

function Base.similar(G::T, ξ::GFSignEnum) where T <: KeldyshTimeGF
  T(G.grid, similar(G.gtr), similar(G.les), ξ)
end

function Base.zero(G::T, ξ::GFSignEnum) where T <: KeldyshTimeGF
  T(G.grid, zero(G.gtr), zero(G.les), ξ)
end

function Base.similar(G::T) where T <: KeldyshTimeGF
  similar(G, G.ξ)
end

function Base.zero(G::T) where T <: KeldyshTimeGF
  zero(G, G.ξ)
end

function TimeDomain(G::KeldyshTimeGF)
  grid = G.grid
  tplus = grid[forward_branch]
  tminus = reverse(grid[backward_branch])

  nt = grid.nt

  points = Tuple{TimeGridPoint, TimeGridPoint}[]

  for i in 1:nt
    for j in 1:i
      push!(points, (tminus[j], tminus[i])) # greater
      push!(points, (tplus[j], tminus[i])) # lesser
    end
  end

  return TimeDomain(points)
end

function KeldyshTimeGF(f::Function, ::Type{T}, grid::KeldyshTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) where T <: Number
  G = KeldyshTimeGF(T, grid, norb, ξ, scalar)
  D = TimeDomain(G)

  for (t1, t2) in D.points
    G[t1, t2] = f(t1, t2)
  end

  return G
end

KeldyshTimeGF(f::Function, grid::KeldyshTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) = KeldyshTimeGF(f, ComplexF64, grid, norb, ξ, scalar)

function KeldyshTimeGF(dos::AbstractDOS, β, grid::KeldyshTimeGrid)
  # result is time invariant so do calculation for time invariant case and then copy over
  G = TimeInvariantKeldyshTimeGF(dos, β, grid)
  KeldyshTimeGF(grid,1,fermionic,true) do t1, t2
    G[t1,t2]
  end
end
