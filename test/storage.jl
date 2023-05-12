using Keldysh, LinearAlgebra, Test

@testset "storage" begin
  let norb=2, N=10
    X = rand(norb, norb, N, N)
    Y = rand(norb, norb, N, N)
    Z = -(X .+ 2 .* Y .- 3 .* X)

    A = Keldysh.GenericStorage(X)
    B = Keldysh.GenericStorage(Y)
    C = -(A + 2*B - 3*A)

    D = similar(A)
    for (i,j) in Iterators.product(1:N,1:N)
      D[i,j] = Z[:,:,i,j]
    end

    E = zero(A)
    for (k,l,i,j) in Iterators.product(1:norb, 1:norb, 1:N, 1:N)
      E[k,l,i,j] = Z[k,l,i,j]
    end

    @test size(A) == (norb, norb, N, N)

    # test getindex
    for (i,j) in Iterators.product(1:N,1:N)
      @test C[i,j] == Z[:,:,i,j]
      @test D[i,j] == Z[:,:,i,j]
      @test E[i,j] == Z[:,:,i,j]
    end
    for (k,l,i,j) in Iterators.product(1:norb, 1:norb, 1:N, 1:N)
      @test C[k,l,i,j] == Z[k,l,i,j]
      @test D[k,l,i,j] == Z[k,l,i,j]
      @test E[k,l,i,j] == Z[k,l,i,j]
    end
  end

  let norb=2, N=10
    get_rand_herm = () -> begin
      X = rand(ComplexF64, norb, norb, N, N)
      X = im * (X + conj(permutedims(X, (2,1,4,3))))
      return X
    end

    X = get_rand_herm()
    Y = get_rand_herm()

    Z = -(X .+ 2 .* Y .- 3 .* X)

    A = Keldysh.AntiHermitianStorage(X)
    B = Keldysh.AntiHermitianStorage(Y)
    C = -(A + 2*B - 3*A)

    Dup = similar(A)
    Dlo = similar(A)
    for (i,j) in Iterators.product(1:N, 1:N)
      i >= j && (Dup[i,j] = Z[:,:,i,j])
      i <= j && (Dlo[i,j] = Z[:,:,i,j])
    end

    Eup = zero(A)
    Elo = zero(A)
    for (k,l,i,j) in Iterators.product(1:norb, 1:norb, 1:N, 1:N)
      i >= j && (Eup[k,l,i,j] = Z[k,l,i,j])
      i <= j && (Elo[k,l,i,j] = Z[k,l,i,j])
    end

    for (i,j) in Iterators.product(1:N, 1:N)
      @test C[i,j] == Z[:,:,i,j]
      @test Dup[i,j] == Z[:,:,i,j]
      @test Dlo[i,j] == Z[:,:,i,j]
      @test Eup[i,j] == Z[:,:,i,j]
      @test Elo[i,j] == Z[:,:,i,j]
    end

    for (k,l,i,j) in Iterators.product(1:norb, 1:norb, 1:N, 1:N)
      @test C[k,l,i,j] == Z[k,l,i,j]
      @test Dup[k,l,i,j] == Z[k,l,i,j]
      @test Dlo[k,l,i,j] == Z[k,l,i,j]
      @test Eup[k,l,i,j] == Z[k,l,i,j]
      @test Elo[k,l,i,j] == Z[k,l,i,j]
    end
  end

  let norb=2, N=10
    get_rand_herm_toeplitz = () -> begin
      v = rand(ComplexF64, norb, norb, N)
      v[:,:,1] = im*(v[:,:,1] + v[:,:,1]')
      X = zeros(ComplexF64, norb, norb, N, N)
      for (i,j) in Iterators.product(1:N, 1:N)
        X[:,:,i,j] = i >= j ? v[:,:,i-j+1] : -v[:,:,j-i+1]'
      end
      return X
    end

    X = get_rand_herm_toeplitz()
    Y = get_rand_herm_toeplitz()

    Z = -(X .+ 2 .* Y .- 3 .* X)

    A = Keldysh.AntiHermitianToeplitzStorage(X)
    B = Keldysh.AntiHermitianToeplitzStorage(Y)
    C = -(A + 2*B - 3*A)

    Dup = similar(A)
    Dlo = similar(A)
    for i in 1:N
      Dup[i,1] = Z[:,:,i,1]
      Dlo[1,i] = Z[:,:,1,i]
    end

    Eup = zero(A)
    Elo = zero(A)
    for (k,l,i) in Iterators.product(1:norb, 1:norb, 1:N)
      Eup[k,l,i,1] = Z[k,l,i,1]
      Elo[k,l,1,i] = Z[k,l,1,i]
    end

    for (i,j) in Iterators.product(1:N, 1:N)
      @test C[i,j] == Z[:,:,i,j]
      @test Dup[i,j] == Z[:,:,i,j]
      @test Dlo[i,j] == Z[:,:,i,j]
      @test Eup[i,j] == Z[:,:,i,j]
      @test Elo[i,j] == Z[:,:,i,j]
    end

    for (k,l,i,j) in Iterators.product(1:norb, 1:norb, 1:N, 1:N)
      @test C[k,l,i,j] == Z[k,l,i,j]
      @test Dup[k,l,i,j] == Z[k,l,i,j]
      @test Dlo[k,l,i,j] == Z[k,l,i,j]
      @test Eup[k,l,i,j] == Z[k,l,i,j]
      @test Elo[k,l,i,j] == Z[k,l,i,j]
    end
  end

  let norb=2, N=10
    get_rand_periodic = () -> begin
      v = rand(norb, norb, N)
      X = zeros(norb, norb, N, N)
      for (i,j) in Iterators.product(1:N, 1:N)
        X[:,:,i,j] = i >= j ? v[:,:,i-j+1] : v[:,:,i-j+N]
      end
      return X
    end

    X = get_rand_periodic()
    Y = get_rand_periodic()
    Z = -(X .+ 2 .* Y .- 3 .* X)

    A = Keldysh.PeriodicStorage(X[:,:,:,1])
    B = Keldysh.PeriodicStorage(Y[:,:,:,1])
    C = -(A + 2*B - 3*A)

    Dup = similar(A)
    for i in 1:N
      Dup[i,1] = Z[:,:,i,1]
    end

    Eup = zero(A)
    for (k,l,i) in Iterators.product(1:norb, 1:norb, 1:N)
      Eup[k,l,i,1] = Z[k,l,i,1]
    end

    for (i,j) in Iterators.product(1:N, 1:N)
      @test C[i,j] == Z[:,:,i,j]
      @test Dup[i,j] == Z[:,:,i,j]
      @test Eup[i,j] == Z[:,:,i,j]
    end

    for (k,l,i,j) in Iterators.product(1:norb, 1:norb, 1:N, 1:N)
      @test C[k,l,i,j] == Z[k,l,i,j]
      @test Dup[k,l,i,j] == Z[k,l,i,j]
      @test Eup[k,l,i,j] == Z[k,l,i,j]
    end

  end
end
