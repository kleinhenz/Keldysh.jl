struct TimeInvariantKeldyshTimeGF{T, scalar} <: AbstractTimeGF{T}
  grid::TimeGrid
  gtr::TimeInvariantAntiHermitianStorage{T,scalar}
  les::TimeInvariantAntiHermitianStorage{T,scalar}
  nt::Int # TODO move to grid

  function TimeInvariantKeldyshTimeGF(grid::TimeGrid,
                      gtr::TimeInvariantAntiHermitianStorage{T,scalar},
                      les::TimeInvariantAntiHermitianStorage{T,scalar}) where {T, scalar}
    @assert grid.contour.domain == keldysh_contour
    nt = length(grid, forward_branch)
    new{T, scalar}(grid, gtr, les, nt)
  end
end

norbitals(G::TimeInvariantKeldyshTimeGF) = G.gtr.norb

function TimeInvariantKeldyshTimeGF(::Type{T}, grid::TimeGrid, norb=1, scalar=false) where T <: Number
  nt = length(grid, forward_branch)

  gtr = TimeInvariantAntiHermitianStorage(T, nt, norb, scalar)
  les = TimeInvariantAntiHermitianStorage(T, nt, norb, scalar)

  TimeInvariantKeldyshTimeGF(grid, gtr, les)
end
TimeInvariantKeldyshTimeGF(grid::TimeGrid, norb=1, scalar=false) = TimeInvariantKeldyshTimeGF(ComplexF64, grid, norb, scalar)

function Base.getindex(G::TimeInvariantKeldyshTimeGF, t1::TimeGridPoint, t2::TimeGridPoint, greater=true)
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

function TimeInvariantKeldyshTimeGF(f::Function, ::Type{T}, grid::TimeGrid, norb=1, scalar=false) where T <: Number
  G = TimeInvariantKeldyshTimeGF(T, grid, norb, scalar)
  D = TimeDomain(G)

  for (t1, t2) in D.points
    G[t1, t2] = f(t1, t2)
  end

  return G
end

TimeInvariantKeldyshTimeGF(f::Function, grid::TimeGrid, norb=1, scalar=false) = TimeInvariantKeldyshTimeGF(f, ComplexF64, grid, norb, scalar)

function TimeInvariantKeldyshTimeGF(dos::AbstractDOS, β, grid::TimeGrid)
  TimeInvariantKeldyshTimeGF(grid, 1, true) do t1, t2
    dos2gf(dos, β, t1.val, t2.val)
  end
end
