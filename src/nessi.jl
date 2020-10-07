function nessi_read_les(h5g::HDF5Group)
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

function nessi_read_ret(h5g::HDF5Group)
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

function nessi_read_tv(h5g::HDF5Group)
  nt = read(h5g, "nt")[]
  ntau = read(h5g, "ntau")[]
  norb = read(h5g, "size1")[]

  tv_ = read(h5g, "tv")
  tv = permutedims(reshape(tv_, (norb, norb, ntau+1, nt+1)),(1,2,4,3))
  return tv
end

function nessi_read_mat(h5g::HDF5Group)
  mat = read(h5g, "mat")
  return mat
end

function nessi_read_gf(h5g::HDF5Group)
  les = nessi_read_les(h5g)
  ret = nessi_read_ret(h5g)
  tv  = nessi_read_tv(h5g)
  mat = nessi_read_mat(h5g)
  return (;les, ret, tv, mat)
end

function TimeGF(les::AbstractMatrix, ret::AbstractMatrix, tv::AbstractMatrix, mat::AbstractVector, grid::TimeGrid)
  @assert grid.contour.domain == full_contour

  t = realtimes(grid)
  tau = imagtimes(grid)

  nt = length(t)
  ntau = length(tau)

  @assert size(les) == (nt,nt)
  @assert size(ret) == (nt,nt)
  @assert size(tv) == (nt,ntau)
  @assert size(mat) == (ntau,)

  G = TimeGF(grid) do t1, t2
    greater = heaviside(t1.val, t2.val)
    if ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) &&
        (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
      i = findfirst(isapprox(real(t1.val.val)), t)
      j = findfirst(isapprox(real(t2.val.val)), t)
      greater ? ret[i,j] + les[i,j] : les[i,j]
    elseif (t1.val.domain == imaginary_branch && (t2.val.domain == forward_branch || t2.val.domain == backward_branch))
      i = findfirst(isapprox(imag(t1.val.val)), -tau)
      j = findfirst(isapprox(real(t2.val.val)), t)
      conj(tv[j, ntau+1-i]) # akoi 19c
    elseif ((t1.val.domain == forward_branch || t1.val.domain == backward_branch) && t2.val.domain == imaginary_branch)
      i = findfirst(isapprox(real(t1.val.val)), t)
      j = findfirst(isapprox(imag(t2.val.val)), -tau)
      tv[i,j]
    else
      i = findfirst(isapprox(imag(t1.val.val)), -tau)
      j = findfirst(isapprox(imag(t2.val.val)), -tau)
      greater ? 1.0im * mat[i - j + 1] : -1.0im * mat[i - j + ntau]
    end
  end

  return G
end
