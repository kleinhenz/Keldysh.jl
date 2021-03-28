using LinearAlgebra
using BenchmarkTools
using Keldysh

suite = BenchmarkGroup()

function setup_indexing_benchmark()
  suite["indexing"] = BenchmarkGroup()

  c = FullContour(tmax=2.0, β=5.0)
  grid = FullTimeGrid(c, 41, 101)
  N = length(grid)

  X = zeros(1,1,N,N)
  Y = Keldysh.GenericStorage(N,N,1,true)

  G1 = GenericTimeGF(grid, 1, true)
  G2 = FullTimeGF(grid, 1, fermionic, true)

  b = @benchmarkable $X[1,1,i,j] setup=(i=rand(1:$N); j=rand(1:$N))
  suite["indexing"]["Array"] = b

  b = @benchmarkable $Y[i,j] setup=(i=rand(1:$N); j=rand(1:$N))
  suite["indexing"]["GenericStorage"] = b

  b = @benchmarkable $G1[t1,t2] setup=(t1=$grid[rand(1:$N)]; t2=$grid[rand(1:$N)])
  suite["indexing"]["GenericTimeGF"] = b

  b = @benchmarkable $G2[t1,t2] setup=(t1=$grid[rand(1:$N)]; t2=$grid[rand(1:$N)])
  suite["indexing"]["FullTimeGF"] = b

  tune!(suite["indexing"])
end

function main()
  setup_indexing_benchmark()

  compare = false
  if length(ARGS) > 0

    if length(ARGS) != 2
      error("usage: julia test/benchmark.jl test/params.ref.json test/results.ref.json")
    end

    params_file = ARGS[1]
    p = BenchmarkTools.load(params_file)[1]
    loadparams!(suite, p, :evals, :samples)

    ref_file = ARGS[2]
    ref = BenchmarkTools.load(ref_file)[1]
    compare = true
  end

  results = run(suite, verbose=true)
  display(median(results))
  println()

  BenchmarkTools.save("params.json", params(suite))
  BenchmarkTools.save("results.json", results)

  if compare
    cmp = judge(median(results), median(ref))
    display(cmp)
  end
end

main()
