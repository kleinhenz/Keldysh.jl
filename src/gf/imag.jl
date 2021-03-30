struct ImaginaryTimeGF{T, scalar} <: AbstractTimeGF{T}
  grid::ImaginaryTimeGrid
  mat::CirculantStorage{T,scalar}
  ξ::GFSignEnum

  function ImaginaryTimeGF(grid::ImaginaryTimeGrid,
                           mat::CirculantStorage{T, scalar},
                           ξ::GFSignEnum=fermionic) where {T, scalar}
    new{T, scalar}(grid, mat, ξ)
  end
end

norbitals(G::ImaginaryTimeGF) = G.mat.norb

function ImaginaryTimeGF(::Type{T}, grid::ImaginaryTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) where T <: Number
  ntau = grid.ntau
  mat = CirculantStorage(T, ntau, norb, scalar)
  ImaginaryTimeGF(grid, mat, ξ)
end
ImaginaryTimeGF(grid::ImaginaryTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) = ImaginaryTimeGF(ComplexF64, grid, norb, ξ, scalar)

@inline function Base.getindex(G::ImaginaryTimeGF, t1::TimeGridPoint, t2::TimeGridPoint, greater=true)
  greater = t1 == t2 ? greater : heaviside(t1.val, t2.val)

  i = t1.ridx
  j = t2.ridx
  ξ = Int(G.ξ)
  ntau = G.grid.ntau

  greater ? G.mat[i,j, greater] : ξ * G.mat[i,j, greater]
end

function Base.setindex!(G::ImaginaryTimeGF, v, t1::TimeGridPoint, t2::TimeGridPoint)
  greater = heaviside(t1.val, t2.val)

  i = t1.ridx
  j = t2.ridx
  ξ = Int(G.ξ)

  if greater
    G.mat[i,j] = v
  else
    G.mat[i,j] = ξ * v
  end
end

function TimeDomain(G::ImaginaryTimeGF)
  grid = G.grid
  tau = grid[imaginary_branch]

  ntau = G.grid.ntau

  points = Tuple{TimeGridPoint, TimeGridPoint}[]

  for i in 1:ntau
    push!(points, (tau[i], tau[1])) # matsubara
  end

  return TimeDomain(points)
end

function ImaginaryTimeGF(f::Function, ::Type{T}, grid::ImaginaryTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) where T <: Number
  G = ImaginaryTimeGF(T, grid, norb, ξ, scalar)
  D = TimeDomain(G)

  for (t1, t2) in D.points
    G[t1, t2] = f(t1, t2)
  end

  return G
end

ImaginaryTimeGF(f::Function, grid::ImaginaryTimeGrid, norb=1, ξ::GFSignEnum=fermionic, scalar=false) = ImaginaryTimeGF(f, ComplexF64, grid, norb, ξ, scalar)

function ImaginaryTimeGF(dos::AbstractDOS, grid::ImaginaryTimeGrid)
  β = grid.contour.β
  ImaginaryTimeGF(grid, 1, fermionic, true) do t1, t2
    Keldysh.dos2gf(dos, β, t1.val, t2.val)
  end
end
