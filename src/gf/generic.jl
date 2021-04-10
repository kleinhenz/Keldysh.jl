struct GenericTimeGF{T, scalar, U <: AbstractTimeGrid} <: AbstractTimeGF{T, scalar}
  grid::U
  data::GenericStorage{T, scalar}
end

norbitals(G::GenericTimeGF) = G.data.norb

function GenericTimeGF(::Type{T}, grid::AbstractTimeGrid, norb=1, scalar=false) where T <: Number
  N = length(grid)
  data = GenericStorage(T, N, N, norb, scalar)
  GenericTimeGF(grid, data)
end
GenericTimeGF(grid::AbstractTimeGrid, norb=1, scalar=false) = GenericTimeGF(ComplexF64, grid, norb, scalar)

function GenericTimeGF(f::Function, ::Type{T}, grid::AbstractTimeGrid, norb=1, scalar=false) where T <: Number
  N = length(grid)
  G = GenericTimeGF(T, grid, norb, scalar)

  for t1 in grid
    for t2 in grid
      G[t1, t2] = f(t1, t2)
    end
  end
  return G
end
GenericTimeGF(f::Function, grid::AbstractTimeGrid, norb=1, scalar=false) = GenericTimeGF(f, ComplexF64, grid, norb, scalar)

@inline function Base.getindex(G::GenericTimeGF, t1::TimeGridPoint, t2::TimeGridPoint, greater=true)
  val = G.data[t1.idx, t2.idx]
  (!greater && t1.idx == t2.idx) && (val += jump(G))
  return val
end

function Base.setindex!(G::GenericTimeGF, v, t1::TimeGridPoint, t2::TimeGridPoint)
  G.data[t1.idx, t2.idx] = v
end

function jump(G::GenericTimeGF)
  t0_plus = branch_bounds(G.grid, forward_branch)[1]
  t0_minus = branch_bounds(G.grid, backward_branch)[2]
  return G[t0_plus, t0_minus] - G[t0_plus, t0_plus]
end

function Base.similar(G::T) where T <: GenericTimeGF
  T(G.grid, similar(G.data))
end

function Base.zero(G::T) where T <: GenericTimeGF
  T(G.grid, zero(G.data))
end

function GenericTimeGF(dos::AbstractDOS, grid::FullTimeGrid)
  G = TimeInvariantFullTimeGF(dos, grid)
  GenericTimeGF(grid,1,true) do t1, t2
    G[t1,t2]
  end
end

function GenericTimeGF(dos::AbstractDOS, grid::ImaginaryTimeGrid)
  G = ImaginaryTimeGF(dos, grid)
  GenericTimeGF(grid,1,true) do t1, t2
    G[t1,t2]
  end
end

function GenericTimeGF(dos::AbstractDOS, β, grid::KeldyshTimeGrid)
  G = TimeInvariantKeldyshTimeGF(dos, β, grid)
  GenericTimeGF(grid,1,true) do t1, t2
    G[t1,t2]
  end
end
