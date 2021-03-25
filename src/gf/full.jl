struct FullTimeGF{T, scalar} <: AbstractTimeGF{T}
  grid::TimeGrid
  gtr::AntiHermitianStorage{T,scalar}
  les::AntiHermitianStorage{T,scalar}
  rm::GenericStorage{T,scalar}
  mat::ImaginaryTimeStorage{T,scalar}
  ntau::Int # TODO move to grid
  nt::Int # TODO move to grid
  ξ::Int

  function FullTimeGF(grid::TimeGrid,
                      gtr::AntiHermitianStorage{T,scalar},
                      les::AntiHermitianStorage{T,scalar},
                      rm::GenericStorage{T,scalar},
                      mat::ImaginaryTimeStorage{T,scalar},
                      ξ=-1) where {T, scalar}

    @assert grid.contour.domain == full_contour
    @assert ξ == -1 || ξ == 1

    nt = length(grid, forward_branch)
    ntau = length(grid, imaginary_branch)
    new{T, scalar}(grid, gtr, les, rm, mat, ntau, nt, ξ)
  end
end

norbitals(G::FullTimeGF) = G.gtr.norb

function FullTimeGF(::Type{T}, grid::TimeGrid, norb=1, ξ=-1, scalar=false) where T <: Number
  nt = length(grid, forward_branch)
  ntau = length(grid, imaginary_branch)

  gtr = AntiHermitianStorage(T, nt, norb, scalar)
  les = AntiHermitianStorage(T, nt, norb, scalar)
  rm = GenericStorage(T, ntau, nt, norb, scalar)
  mat = ImaginaryTimeStorage(T, ntau, norb, scalar)

  FullTimeGF(grid, gtr, les, rm, mat, ξ)
end
FullTimeGF(grid::TimeGrid, norb=1, ξ=-1, scalar=false) = FullTimeGF(ComplexF64, grid, norb, ξ, scalar)

function Base.getindex(G::FullTimeGF, t1::TimeGridPoint, t2::TimeGridPoint, greater=true)
  greater = t1 == t2 ? greater : heaviside(t1.val, t2.val)

  i = t1.ridx
  j = t2.ridx
  ξ = G.ξ

  if ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) &&
      (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
    return greater ? G.gtr[i,j] : G.les[i,j]
  elseif (t1.val.domain == imaginary_branch && (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
    return G.rm[i,j]
  elseif ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) && t2.val.domain == imaginary_branch)
    return -ξ * conj(G.rm[G.ntau+1-j,i]) # akoi 19c
  else
    greater ? G.mat[i,j] : ξ * G.mat[i,j]
  end
end

function Base.setindex!(G::FullTimeGF, v, t1::TimeGridPoint, t2::TimeGridPoint)
  greater = heaviside(t1.val, t2.val)

  i = t1.ridx
  j = t2.ridx
  ξ = G.ξ

  if ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) &&
      (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
    return greater ? G.gtr[i,j] = v : G.les[i,j] = v
  elseif (t1.val.domain == imaginary_branch && (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
    return G.rm[i,j] = v
  elseif ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) && t2.val.domain == imaginary_branch)
    return G.rm[ntau+1-j,i] = -ξ * conj(v) #akoi 19c
  else
    if greater
      G.mat[i,j] = v
    else
      G.mat[i,j] = ξ * v
    end
  end
end

function TimeDomain(G::FullTimeGF)
  grid = G.grid
  tplus = grid[forward_branch]
  tminus = reverse(grid[backward_branch])
  tau = grid[imaginary_branch]

  nt = length(tplus)
  ntau = length(tau)

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

function FullTimeGF(f::Function, ::Type{T}, grid::TimeGrid, norb=1, ξ=-1, scalar=false) where T <: Number
  G = FullTimeGF(T, grid, norb, ξ, scalar)
  D = TimeDomain(G)

  for (t1, t2) in D.points
    G[t1, t2] = f(t1, t2)
  end

  return G
end

FullTimeGF(f::Function, grid::TimeGrid, norb=1, ξ=-1, scalar=false) = FullTimeGF(f, ComplexF64, grid, norb, ξ, scalar)

function FullTimeGF(dos::AbstractDOS, grid::TimeGrid)
  # result is time invariant so do calculation for time invariant case and then copy over
  G = TimeInvariantFullTimeGF(dos, grid)
  FullTimeGF(grid,1,-1,true) do t1, t2
    G[t1,t2]
  end
end
