struct NessiGFData
  les::Array{ComplexF64, 4}
  ret::Array{ComplexF64, 4}
  tv::Array{ComplexF64, 4}
  mat::Array{ComplexF64, 3}
  nt::Int
  ntau::Int
  norb::Int
  sig::Int
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
  sig = read(h5g, "sig")[]

  return NessiGFData(les, ret, tv, mat, nt, ntau, norb, sig)
end

function FullTimeGF(data::NessiGFData, grid::FullTimeGrid)
  nt = grid.nt
  ntau = grid.ntau

  @assert nt == data.nt+1
  @assert ntau == data.ntau+1
  ξ = data.sig == 1 ? bosonic : fermionic

  G = FullTimeGF(grid, data.norb, ξ, false)

  tplus = grid[forward_branch]
  tminus = reverse(grid[backward_branch])
  tau = grid[imaginary_branch]

  # lesser
  for i in 1:nt
    for j in 1:i
      t1 = tplus[j]
      t2 = tminus[i]
      G[t1,t2] = data.les[:,:,j,i]
    end
  end

  # greater
  for i in 1:nt
    for j in 1:i
      t1 = tminus[j]
      t2 = tminus[i]
      G[t1,t2] = (data.les[:,:,j,i] .+ data.ret[:,:,j,i])
    end
  end

  # matsubara
  for i in 1:ntau
    t1 = tau[i]
    t2 = tau[1]
    G[t1,t2] = 1.0im * data.mat[:,:,i]
  end


  # left-mixing
  # NOTE: NessiGFData stores left-mixing while FullTimeGF stores right-mixing
  # We rely on setindex! method to apply symmetry Hermitian symmetry
  for i in 1:nt
    for j in 1:ntau
      t1 = tminus[i]
      t2 = tau[j]
      G[t1, t2] = data.tv[:,:,i,j]
    end
  end

  return G
end
