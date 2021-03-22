struct TimeGF{T, scalar} <: AbstractTimeGF{T}
  grid::TimeGrid
  data::Array{T, 4}

  function TimeGF(grid::TimeGrid, data::Array{T,4}, scalar=false) where T
    norb = size(data,1)
    N = size(data,3)
    @assert size(data) == (norb, norb, N, N)
    @assert !(scalar && norb > 1)
    norb = size(data,1)
    return new{T,scalar}(grid, data)
  end
end

norbitals(G::TimeGF) = size(G.data,1)

function TimeGF(::Type{T}, grid::TimeGrid, norb=1, scalar=false) where T <: Number
  N = length(grid)
  data = zeros(T, norb, norb, N, N)
  TimeGF(grid, data, scalar)
end
TimeGF(grid::TimeGrid, norb=1, scalar=false) = TimeGF(ComplexF64, grid, norb, scalar)

function TimeGF(f::Function, ::Type{T}, grid::TimeGrid, norb=1, scalar=false) where T <: Number
  N = length(grid)
  G = TimeGF(T, grid, norb, scalar)

  for t1 in grid
    for t2 in grid
      G[t1, t2] = f(t1, t2)
    end
  end
  return G
end
TimeGF(f::Function, grid::TimeGrid, norb=1, scalar=false) = TimeGF(f, ComplexF64, grid, norb, scalar)

AbstractArray4 = AbstractArray{T,4} where T
AbstractArray3 = AbstractArray{T,3} where T

function TimeGF(les::AbstractArray4,
                ret::AbstractArray4,
                tv::AbstractArray4,
                mat::AbstractArray3,
                grid::TimeGrid,
                scalar=false)

  @assert grid.contour.domain == full_contour

  t = realtimes(grid)
  tau = imagtimes(grid)

  nt = length(t)
  ntau = length(tau)
  norb = size(les, 1)

  @assert size(les) == (norb, norb, nt, nt)
  @assert size(ret) == (norb, norb, nt, nt)
  @assert size(tv) == (norb, norb, nt, ntau)
  @assert size(mat) == (norb, norb, ntau)

  G = TimeGF(grid, norb, scalar) do t1, t2
    greater = heaviside(t1.val, t2.val)
    i = t1.ridx
    j = t2.ridx

    x =
    if ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) &&
        (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
      greater ? ret[:,:,i,j] .+ les[:,:,i,j] : les[:,:,i,j]
    elseif (t1.val.domain == imaginary_branch && (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
      conj(tv[:,:,j, ntau+1-i]) # akoi 19c
    elseif ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) && t2.val.domain == imaginary_branch)
      tv[:,:,i,j]
    else
      greater ? 1.0im * mat[:,:,i - j + 1] : -1.0im * mat[:,:,i - j + ntau]
    end

    return norb == 1 ? x[] : x
  end

  return G
end

function TimeGF(les::AbstractArray4,
                ret::AbstractArray4,
                grid::TimeGrid,
                scalar=false)

  @assert grid.contour.domain == keldysh_contour

  t = realtimes(grid)

  nt = length(t)
  norb = size(les, 1)

  @assert size(les) == (norb,norb,nt,nt)
  @assert size(ret) == (norb,norb,nt,nt)

  G = TimeGF(grid, norb, scalar) do t1, t2
    greater = heaviside(t1.val, t2.val)
    i = t1.ridx
    j = t2.ridx
    x = greater ? ret[:,:,i,j] .+ les[:,:,i,j] : les[:,:,i,j]
    return norb == 1 ? x[] : x
  end

  return G
end

#### Indexing ###
Base.@propagate_inbounds function Base.getindex(G::TimeGF{T,scalar}, i::Int, j::Int) where {T,scalar}
  if scalar
    return G.data[1,1,i,j]
  else
    return G.data[:,:,i,j]
  end
end

Base.@propagate_inbounds function Base.setindex!(G::TimeGF{T,scalar}, v, i::Int, j::Int) where {T,scalar}
  if scalar
    G.data[1, 1, i, j] = v
  else
    G.data[:,:,i,j] = v
  end
end

# indexing with TimeGridPoint
function Base.getindex(G::TimeGF, t1::TimeGridPoint, t2::TimeGridPoint, greater=true)
  val = @inbounds G[t1.idx, t2.idx]
  (!greater && t1.idx == t2.idx) && (val += jump(G))
  return val
end

function Base.setindex!(G::TimeGF, v, t1::TimeGridPoint, t2::TimeGridPoint)
  @inbounds G[t1.idx, t2.idx] = v
end

function jump(G::TimeGF)
  t0_plus = branch_bounds(G.grid, forward_branch)[1]
  t0_minus = branch_bounds(G.grid, backward_branch)[2]
  return G[t0_plus, t0_minus] - G[t0_plus, t0_plus]
end
