struct TimeInvariantFullTimeGF{T, scalar} <: AbstractTimeGF{T}
  grid::TimeGrid
  gtr::TimeInvariantAntiHermitianStorage{T,scalar}
  les::TimeInvariantAntiHermitianStorage{T,scalar}
  rm::GenericStorage{T,scalar}
  mat::ImaginaryTimeStorage{T,scalar}
  ntau::Int # TODO move to grid
  nt::Int # TODO move to grid

  function TimeInvariantFullTimeGF(grid::TimeGrid,
                      gtr::TimeInvariantAntiHermitianStorage{T,scalar},
                      les::TimeInvariantAntiHermitianStorage{T,scalar},
                      rm::GenericStorage{T,scalar},
                      mat::ImaginaryTimeStorage{T,scalar}) where {T, scalar}

    nt = length(grid, forward_branch)
    ntau = length(grid, imaginary_branch)
    new{T, scalar}(grid, gtr, les, rm, mat, ntau, nt)
  end
end

norbitals(G::TimeInvariantFullTimeGF) = G.gtr.norb

function TimeInvariantFullTimeGF(::Type{T}, grid::TimeGrid, norb=1, scalar=false) where T <: Number
  nt = length(grid, forward_branch)
  ntau = length(grid, imaginary_branch)

  gtr = TimeInvariantAntiHermitianStorage(T, nt, norb, scalar)
  les = TimeInvariantAntiHermitianStorage(T, nt, norb, scalar)
  rm = GenericStorage(T, ntau, nt, norb, scalar)
  mat = ImaginaryTimeStorage(T, ntau, norb, scalar)

  TimeInvariantFullTimeGF(grid, gtr, les, rm, mat)
end
TimeInvariantFullTimeGF(grid::TimeGrid, norb=1, scalar=false) = TimeInvariantFullTimeGF(ComplexF64, grid, norb, scalar)

function Base.getindex(G::TimeInvariantFullTimeGF, t1::TimeGridPoint, t2::TimeGridPoint, greater=true)
  greater = t1 == t2 ? greater : heaviside(t1.val, t2.val)

  i = t1.ridx
  j = t2.ridx

  if ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) &&
      (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
    return greater ? G.gtr[i,j] : G.les[i,j]
  elseif (t1.val.domain == imaginary_branch && (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
    return G.rm[i,j]
  elseif ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) && t2.val.domain == imaginary_branch)
    return conj(G.rm[G.ntau+1-j,i]) # akoi 19c
  else
    i == j && !greater ? -G.mat[i,j] : G.mat[i,j]
  end
end

function Base.setindex!(G::TimeInvariantFullTimeGF, v, t1::TimeGridPoint, t2::TimeGridPoint)
  greater = heaviside(t1.val, t2.val)

  i = t1.ridx
  j = t2.ridx

  if ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) &&
      (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
    return greater ? G.gtr[i,j] = v : G.les[i,j] = v
  elseif (t1.val.domain == imaginary_branch && (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
    return G.rm[i,j] = v
  elseif ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) && t2.val.domain == imaginary_branch)
    return G.rm[ntau+1-j,i] = conj(v) #akoi 19c
  else
    G.mat[i,j] = v
  end
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

function TimeInvariantFullTimeGF(f::Function, ::Type{T}, grid::TimeGrid, norb=1, scalar=false) where T <: Number
  G = TimeInvariantFullTimeGF(T, grid, norb, scalar)
  D = TimeDomain(G)

  for (t1, t2) in D.points
    G[t1, t2] = f(t1, t2)
  end

  return G
end

TimeInvariantFullTimeGF(f::Function, grid::TimeGrid, norb=1, scalar=false) = TimeInvariantFullTimeGF(f, ComplexF64, grid, norb, scalar)

function TimeInvariantFullTimeGF(dos::AbstractDOS, grid::TimeGrid)
  t = realtimes(grid)
  τ = imagtimes(grid)
  β = τ[end]

  nt = length(t)
  ntau = length(τ)

  TimeInvariantFullTimeGF(grid, 1, true) do t1, t2
    dos2gf(dos, β, t1.val, t2.val)
  end
end
