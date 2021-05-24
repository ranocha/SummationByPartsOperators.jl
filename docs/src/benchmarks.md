# Benchmarks

Here are some simple benchmarks. Take them with a grain of salt since they run
on virtual machines in the cloud to generate the documentation automatically.


## First-derivative operators

```@example
using BenchmarkTools
using SummationByPartsOperators
using BandedMatrices, LinearAlgebra, Random, SparseArrays

BLAS.set_num_threads(1) # make sure that BLAS is serial to be fair

T = Float64
xmin, xmax = T(0), T(1)
N = 10^3
der_order = 1 # first-derivative operators
acc_order = 6 # the (interior) order of accuracy is six
source = MattssonSvärdShoeybi2008()

D_periodic_serial  = periodic_derivative_operator(der_order, acc_order, xmin, xmax, N, Val{:serial}())
D_nonperiodic_serial  = derivative_operator(source, der_order, acc_order, xmin, xmax, N, Val{:serial}())
D_nonperiodic_sparse  = sparse(D_nonperiodic_serial)
D_nonperiodic_banded  = BandedMatrix(D_nonperiodic_serial)

Random.seed!(12345)
u = randn(T, N)
dest = similar(u)

function doit(D, text, dest, u)
  println(text)
  sleep(0.1)
  show(stdout, MIME"text/plain"(), @benchmark mul!($dest, $D, $u))
  println()
end

doit(D_periodic_serial, "D_periodic_serial:", dest, u)
doit(D_nonperiodic_serial, "D_nonperiodic_serial:", dest, u)
doit(D_nonperiodic_sparse, "D_nonperiodic_sparse:", dest, u)
doit(D_nonperiodic_banded, "D_nonperiodic_banded:", dest, u)
```


## Dissipation operators

```@example
using BenchmarkTools
using SummationByPartsOperators
using BandedMatrices, LinearAlgebra, Random, SparseArrays

BLAS.set_num_threads(1) # make sure that BLAS is serial to be fair

T = Float64
xmin, xmax = T(0), T(1)
N = 10^3
acc_order = 8
source_D = MattssonSvärdShoeybi2008()
source_Di = MattssonSvärdNordström2004()

D_serial  = derivative_operator(source_D, 1, acc_order, xmin, xmax, N, Val{:serial}())

Di_serial  = dissipation_operator(source_Di, D_serial)
Di_sparse  = sparse(Di_serial)
Di_full    = Matrix(Di_serial)

Random.seed!(12345)
u = randn(T, N)
dest = similar(u)

function doit(D, text, dest, u)
  println(text)
  sleep(0.1)
  show(stdout, MIME"text/plain"(), @benchmark mul!($dest, $D, $u))
  println()
end

doit(D_serial, "D_serial:", dest, u)
doit(Di_serial, "Di_serial:", dest, u)
doit(Di_sparse, "Di_sparse:", dest, u)
doit(Di_full, "Di_full:", dest, u)
```
