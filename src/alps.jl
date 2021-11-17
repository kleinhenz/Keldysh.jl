struct ALPSComplex
  data::Array{ComplexF64}
end

function ALPSComplex(g::HDF5.Dataset)
  data = read(g)

  # flip all dimensions since data is stored as row-major
  data = permutedims(data, reverse(1:ndims(data)))
  sz = size(data)

  # interpret last dimension as real/complex
  if haskey(attributes(g), "__complex__")
    @assert sz[end] == 2
    R = CartesianIndices(sz[1:end-1])
    data = data[R,1] + 1.0im * data[R, 2]
  end

  return ALPSComplex(data)
end

Base.read(g::HDF5.Group, ::Type{ALPSComplex}) = ALPSComplex(g)

function Base.write(parent::Union{HDF5.File, HDF5.Group}, name::String, X::ALPSComplex)
  # flip all dimensions since data is stored as row-major
  data = permutedims(X.data, reverse(1:ndims(X.data)))

  # reinterpret real/complex as new axis
  data = reshape(reinterpret(Float64, data), 2, size(data)...)

  parent[name] = data
  attributes(parent[name])["__complex__"] = Int8(1)
end

struct ALPSTimeGrid
  grid
  function ALPSTimeGrid(grid::T) where T <: AbstractTimeGrid
    new(grid)
  end
end

function ALPSTimeGrid(g::HDF5.Group)
  nb = read(g, "contour/size")
  branches = Branch[]

  branches = ntuple(nb) do i
    b = BranchEnum(read(g, "contour/branch$(i-1)/type") + 1)
    l = read(g, "contour/branch$(i-1)/len")
    return Branch(b, l)
  end

  values = ALPSComplex(g["values"]).data
  branch_enums = read(g["branch_enums"])

  npts_fwd = count(branch_enums .== 0)
  npts_back = count(branch_enums .== 1)
  ntau = count(branch_enums .== 2)

  @assert npts_fwd == npts_back
  nt = npts_fwd

  get_branch = domain -> begin
    idx = findfirst(b -> b.domain == domain, branches)
    isnothing(idx) ? nothing : branches[idx]
  end

  if nb == 3
    tmax = length(get_branch(forward_branch))
    β = length(get_branch(imaginary_branch))
    c = FullContour(branches, tmax, β)
    grid = FullTimeGrid(c, nt, ntau)
  elseif nb == 2
    tmax = length(get_branch(forward_branch))
    c = KeldyshContour(branches, tmax)
    grid = KeldyshTimeGrid(c, nt)
  else
    β = length(get_branch(imaginary_branch))
    grid = ImaginaryTimeGrid(c, ntau)
  end

  @assert map(t -> t.bpoint.val, grid) ≈ values

  return ALPSTimeGrid(grid)
end

Base.read(g::HDF5.Group, ::Type{ALPSTimeGrid}) = ALPSTimeGrid(g)

function Base.write(parent::Union{HDF5.File, HDF5.Group}, name::String, X::ALPSTimeGrid)
  g = HDF5.create_group(parent, name)

  grid = X.grid
  c = grid.contour

  g["contour/size"] = UInt(nbranches(c))
  for (i, branch) in enumerate(c.branches)
    g["contour/branch$(i-1)/type"] = Int32(Int(branch.domain) - 1)
    g["contour/branch$(i-1)/len"] = length(branch)
  end

  values = ALPSComplex(map(t -> t.bpoint.val, grid))
  refs = map(t -> t.bpoint.ref, grid)
  branch_enums = map(t -> Int32(Int(t.bpoint.domain) - 1), grid)

  write(g, "values", values)
  write(g, "refs", refs)
  write(g, "branch_enums", branch_enums)
end


struct ALPSTimeGF
  G
  function ALPSTimeGF(G::T) where T <: AbstractTimeGF{ComplexF64, true}
    new(G)
  end
end

function ALPSTimeGF(g::HDF5.Group)
  grid = ALPSTimeGrid(g["mesh/1"]).grid
  data = ALPSComplex(g["data"]).data

  G = GenericTimeGF(grid, 1, true) do t1, t2
    data[t1.cidx, t2.cidx]
  end

  return ALPSTimeGF(G)
end

Base.read(g::HDF5.Group, ::Type{ALPSTimeGF}) = ALPSTimeGF(g)

function Base.write(parent::Union{HDF5.File, HDF5.Group}, name::String, X::ALPSTimeGF)
  G = X.G
  grid = G.grid
  N = length(grid)

  data = Matrix{ComplexF64}(undef, N, N)

  for t1 in grid
    for t2 in grid
      data[t1.cidx, t2.cidx] = G[t1, t2]
    end
  end

  g = HDF5.create_group(parent, name)

  write(g, "data", ALPSComplex(data))
  write(g, "mesh/N", 2)

  for i in 1:2
    write(g, "mesh/$i", ALPSTimeGrid(grid))
  end
end
