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

@inline function Base.getindex(G::GenericTimeGF, k, l, t1::TimeGridPoint, t2::TimeGridPoint, greater=true)
  val = G.data[k, l, t1.cidx, t2.cidx]
  (!greater && t1.cidx == t2.cidx) && (val += jump(G))
  return val
end

function Base.setindex!(G::GenericTimeGF, v, k, l, t1::TimeGridPoint, t2::TimeGridPoint)
  G.data[k, l, t1.cidx, t2.cidx] = v
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

function Base.:+(G1::T, G2::T) where T <: GenericTimeGF
  @assert G1.grid == G2.grid
  return T(G1.grid, G1.data + G2.data)
end

function Base.:-(G1::T, G2::T) where T <: GenericTimeGF
  @assert G1.grid == G2.grid
  return T(G1.grid, G1.data - G2.data)
end

function Base.:*(G::T, α::Number) where T <: GenericTimeGF
  return T(G.grid, G.data * α)
end

function Base.:*(α::Number, G::T) where T <: GenericTimeGF
  return G * α
end

function Base.:-(G::T) where T <: GenericTimeGF
  return T(G.grid, -G.data)
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
