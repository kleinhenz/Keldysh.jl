struct KeldyshTimeGF{T, scalar} <: AbstractTimeGF{T}
  grid::TimeGrid
  gtr::AntiHermitianStorage{T,scalar}
  les::AntiHermitianStorage{T,scalar}
  nt::Int # TODO move to grid
  ξ::Int

  function KeldyshTimeGF(grid::TimeGrid,
                      gtr::AntiHermitianStorage{T,scalar},
                      les::AntiHermitianStorage{T,scalar},
                      ξ=-1) where {T, scalar}

    @assert grid.contour.domain == keldysh_contour
    @assert ξ == -1 || ξ == 1

    nt = length(grid, forward_branch)
    new{T, scalar}(grid, gtr, les, nt)
  end
end

norbitals(G::KeldyshTimeGF) = G.gtr.norb

function KeldyshTimeGF(::Type{T}, grid::TimeGrid, norb=1, ξ=-1, scalar=false) where T <: Number
  nt = length(grid, forward_branch)

  gtr = AntiHermitianStorage(T, nt, norb, scalar)
  les = AntiHermitianStorage(T, nt, norb, scalar)

  KeldyshTimeGF(grid, gtr, les, ξ)
end
KeldyshTimeGF(grid::TimeGrid, norb=1, ξ=-1, scalar=false) = KeldyshTimeGF(ComplexF64, grid, norb, ξ, scalar)

function Base.getindex(G::KeldyshTimeGF, t1::TimeGridPoint, t2::TimeGridPoint, greater=true)
  greater = t1 == t2 ? greater : heaviside(t1.val, t2.val)

  i = t1.ridx
  j = t2.ridx

  return greater ? G.gtr[i,j] : G.les[i,j]
end

function Base.setindex!(G::KeldyshTimeGF, v, t1::TimeGridPoint, t2::TimeGridPoint)
  greater = heaviside(t1.val, t2.val)

  i = t1.ridx
  j = t2.ridx

  return greater ? G.gtr[i,j] = v : G.les[i,j] = v
end

function TimeDomain(G::KeldyshTimeGF)
  grid = G.grid
  tplus = grid[forward_branch]
  tminus = reverse(grid[backward_branch])

  nt = length(tplus)

  points = Tuple{TimeGridPoint, TimeGridPoint}[]

  for i in 1:nt
    for j in 1:i
      push!(points, (tminus[j], tminus[i])) # greater
      push!(points, (tplus[j], tminus[i])) # lesser
    end
  end

  return TimeDomain(points)
end

function KeldyshTimeGF(f::Function, ::Type{T}, grid::TimeGrid, norb=1, ξ=-1, scalar=false) where T <: Number
  G = KeldyshTimeGF(T, grid, norb, ξ, scalar)
  D = TimeDomain(G)

  for (t1, t2) in D.points
    G[t1, t2] = f(t1, t2)
  end

  return G
end

KeldyshTimeGF(f::Function, grid::TimeGrid, norb=1, ξ=-1, scalar=false) = KeldyshTimeGF(f, ComplexF64, grid, norb, ξ, scalar)

function KeldyshTimeGF(dos::AbstractDOS, β, grid::TimeGrid)
  # result is time invariant so do calculation for time invariant case and then copy over
  G = TimeInvariantKeldyshTimeGF(dos, β, grid)
  KeldyshTimeGF(grid,1,-1,true) do t1, t2
    G[t1,t2]
  end
end
