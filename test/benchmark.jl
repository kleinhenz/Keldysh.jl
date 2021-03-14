using LinearAlgebra
using BenchmarkTools
using Keldysh

suite = BenchmarkGroup()

function setup_indexing_benchmark()
  c = Contour(full_contour, tmax=2.0, β=5.0)
  grid = TimeGrid(c, npts_real=41, npts_imag=101)
  N = length(grid)

  dos = Keldysh.bethe_dos()
  G = dos2gf(dos, 5.0, grid);

  # set evals=1 so we are benchmarking random accesses rather than repeated accesses to the same element
  b = @benchmarkable $G[t1,t2] setup=(t1=$grid[rand(1:$N)]; t2=$grid[rand(1:$N)]) evals=1 samples=1000000
  tune!(b)
  suite["indexing"] = b
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
