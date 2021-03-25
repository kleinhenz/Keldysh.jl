struct ALPSComplex end

function Base.read(g::HDF5.Dataset, ::Type{ALPSComplex})
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

  return data
end

function write_alpscomplex(parent::Union{HDF5.File, HDF5.Group}, name::String, data::AbstractArray{Complex{T}}) where T <: HDF5.ScalarType
  # flip all dimensions since data is stored as row-major
  data = permutedims(data, reverse(1:ndims(data)))

  # reinterpret real/complex as new axis
  data = reshape(reinterpret(T, data), 2, size(data)...)

  parent[name] = data
  attributes(parent[name])["__complex__"] = Int8(1)
end

function Base.read(g::HDF5.Group, ::Type{Contour})
  nb = read(g, "size")
  branches = Branch[]
  for i in 0:nb-1
    b = BranchEnum(read(g, "branch$i/type") + 1)
    l = read(g, "branch$i/len")
    push!(branches, Branch(b,l))
  end
  return Contour(branches)
end

function Base.write(parent::Union{HDF5.File, HDF5.Group}, name::String, c::Contour)
  g = HDF5.create_group(parent, name)
  g["size"] = UInt(nbranches(c))
  for (i, branch) in enumerate(c.branches)
    g["branch$(i-1)/type"] = Int32(Int(branch.domain) - 1)
    g["branch$(i-1)/len"] = length(branch)
  end
end

function Base.read(g::HDF5.Group, ::Type{TimeGrid})
  c = read(g["contour"], Contour)
  branch_enums = read(g["branch_enums"])
  values = read(g["values"], ALPSComplex)

  npts_fwd = count(branch_enums .== 0)
  npts_back = count(branch_enums .== 1)
  npts_imag = count(branch_enums .== 2)

  @assert npts_fwd == npts_back

  npts_real = npts_fwd
  grid = TimeGrid(c, npts_real=npts_real, npts_imag=npts_imag)

  @assert map(t -> t.val.val, grid) ≈ values

  return grid
end

function Base.write(parent::Union{HDF5.File, HDF5.Group}, name::String, grid::TimeGrid)
  g = HDF5.create_group(parent, name)

  write(g, "contour", grid.contour)

  values = map(t -> t.val.val, grid)
  refs = map(t -> t.val.ref, grid)
  branch_enums = map(t -> Int32(Int(t.val.domain) - 1), grid)

  write_alpscomplex(g, "values", values)
  write(g, "refs", refs)
  write(g, "branch_enums", branch_enums)
end

function Base.read(g::HDF5.Group, ::Type{TimeGF})
  data = read(g["data"], ALPSComplex)
  grid = read(g["mesh/1"], TimeGrid)
  return TimeGF(data, grid)
end

function Base.write(parent::Union{HDF5.File, HDF5.Group}, name::String, G::TimeGF)
  g = HDF5.create_group(parent, name)
  write_alpscomplex(g, "data", G.data)
  write(g, "mesh/1", G.grid)
  write(g, "mesh/2", G.grid)
  write(g, "mesh/N", 2)
end
