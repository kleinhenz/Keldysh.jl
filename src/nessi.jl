struct NessiGFData
  les::Array{ComplexF64, 4}
  ret::Array{ComplexF64, 4}
  tv::Array{ComplexF64, 4}
  mat::Array{ComplexF64, 3}
  nt::Int
  ntau::Int
  norb::Int
end

function _nessi_read_les(h5g::HDF5.Group)
  nt = read(h5g, "nt")[]
  norb = read(h5g, "size1")[]

  les_ = read(h5g, "les")
  les = zeros(ComplexF64, norb, norb, nt+1, nt+1)

  indices = findall(triu(trues(nt+1,nt+1)))
  les[:,:,indices] .= les_

  for i in 1:norb
    for j in 1:norb
      les[i,j,:,:] .-= triu(les[i,j,:,:],1)'
    end
  end

  return les
end

function _nessi_read_ret(h5g::HDF5.Group)
  nt = read(h5g, "nt")[]
  norb = read(h5g, "size1")[]

  ret_ = read(h5g, "ret")
  ret = zeros(ComplexF64, norb, norb, nt+1, nt+1)

  indices = findall(triu(trues(nt+1,nt+1)))
  ret[:,:,indices] .= ret_

  for i in 1:norb
    for j in 1:norb
      ret[i,j,:,:] .-= triu(ret[i,j,:,:],1)'
    end
  end

  ret = permutedims(ret, (1,2,4,3))
  return ret
end

function _nessi_read_tv(h5g::HDF5.Group)
  nt = read(h5g, "nt")[]
  ntau = read(h5g, "ntau")[]
  norb = read(h5g, "size1")[]

  tv_ = read(h5g, "tv")
  tv = permutedims(reshape(tv_, (norb, norb, ntau+1, nt+1)),(1,2,4,3))
  return tv
end

function _nessi_read_mat(h5g::HDF5.Group)
  mat = read(h5g, "mat")
  return mat
end

function NessiGFData(h5g::HDF5.Group)
  les = _nessi_read_les(h5g)
  ret = _nessi_read_ret(h5g)
  tv  = _nessi_read_tv(h5g)
  mat = _nessi_read_mat(h5g)

  nt = read(h5g, "nt")[]
  ntau = read(h5g, "ntau")[]
  norb = read(h5g, "size1")[]

  return NessiGFData(les, ret, tv, mat, nt, ntau, norb)
end

function TimeGF(data::NessiGFData, grid::TimeGrid)
  if grid.contour.domain == full_contour
    return TimeGF(data.les, data.ret, data.tv, data.mat, grid)
  elseif grid.contour.domain == keldysh_contour
    return TimeGF(data.les, data.ret, grid)
  else
    error("unsupported")
  end
end
