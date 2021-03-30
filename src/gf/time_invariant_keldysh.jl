struct TimeInvariantKeldyshTimeGF{T, scalar} <: AbstractTimeGF{T, scalar}
  grid::KeldyshTimeGrid
  gtr::AntiHermitianToeplitzStorage{T,scalar}
  les::AntiHermitianToeplitzStorage{T,scalar}
  ξ::GFSignEnum

  function TimeInvariantKeldyshTimeGF(grid::KeldyshTimeGrid,
                      gtr::AntiHermitianToeplitzStorage{T,scalar},
                      les::AntiHermitianToeplitzStorage{T,scalar},
                      ξ::GFSignEnum=fermionic) where {T, scalar}
    new{T, scalar}(grid, gtr, les, ξ)
  end
end

norbitals(G::TimeInvariantKeldyshTimeGF) = G.gtr.norb

function TimeInvariantKeldyshTimeGF(::Type{T}, grid::KeldyshTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) where T <: Number
  nt = grid.nt

  gtr = AntiHermitianToeplitzStorage(T, nt, norb, scalar)
  les = AntiHermitianToeplitzStorage(T, nt, norb, scalar)

  TimeInvariantKeldyshTimeGF(grid, gtr, les, ξ)
end
TimeInvariantKeldyshTimeGF(grid::KeldyshTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) = TimeInvariantKeldyshTimeGF(ComplexF64, grid, norb, ξ, scalar)

@inline function Base.getindex(G::TimeInvariantKeldyshTimeGF, t1::TimeGridPoint, t2::TimeGridPoint, greater=true)
  greater = t1 == t2 ? greater : heaviside(t1.val, t2.val)

  i = t1.ridx
  j = t2.ridx

  return greater ? G.gtr[i,j] : G.les[i,j]
end

function Base.setindex!(G::TimeInvariantKeldyshTimeGF, v, t1::TimeGridPoint, t2::TimeGridPoint)
  greater = heaviside(t1.val, t2.val)

  i = t1.ridx
  j = t2.ridx

  return greater ? G.gtr[i,j] = v : G.les[i,j] = v
end

function TimeDomain(G::TimeInvariantKeldyshTimeGF)
  grid = G.grid
  tplus = grid[forward_branch]
  tminus = reverse(grid[backward_branch])

  nt = length(tplus)

  points = Tuple{TimeGridPoint, TimeGridPoint}[]

  for i in 1:nt
    push!(points, (tminus[1], tminus[i])) # greater
    push!(points, (tplus[1], tminus[i])) # lesser
  end

  return TimeDomain(points)
end

function TimeInvariantKeldyshTimeGF(f::Function, ::Type{T}, grid::KeldyshTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) where T <: Number
  G = TimeInvariantKeldyshTimeGF(T, grid, norb, ξ, scalar)
  D = TimeDomain(G)

  for (t1, t2) in D.points
    G[t1, t2] = f(t1, t2)
  end

  return G
end

TimeInvariantKeldyshTimeGF(f::Function, grid::KeldyshTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) = TimeInvariantKeldyshTimeGF(f, ComplexF64, grid, norb, ξ, scalar)

function TimeInvariantKeldyshTimeGF(dos::AbstractDOS, β, grid::KeldyshTimeGrid)
  TimeInvariantKeldyshTimeGF(grid, 1, fermionic, true) do t1, t2
    Keldysh.dos2gf(dos, β, t1.val, t2.val)
  end
end
