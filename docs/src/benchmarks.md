# Benchmarks

Here are some simple benchmarks. Take them with a grain of salt since they run
on virtual machines in the cloud to generate the documentation automatically.


## First-derivative operators

Periodic domains:
```@example
using BenchmarkTools
using LinearAlgebra, SparseArrays
using SummationByPartsOperators, DiffEqOperators

BLAS.set_num_threads(1) # make sure that BLAS is serial to be fair

T = Float64
xmin, xmax = T(0), T(1)

D_SBP = periodic_derivative_operator(derivative_order=1, accuracy_order=2,
                                     xmin=xmin, xmax=xmax, N=100)
x = grid(D_SBP)
D_DEO = CenteredDifference(derivative_order(D_SBP), accuracy_order(D_SBP),
                           step(x), length(x)) * PeriodicBC(eltype(D_SBP))

D_sparse = sparse(D_SBP)

u = randn(eltype(D_SBP), length(x)); du = similar(u);
@show D_SBP * u ≈ D_DEO * u ≈ D_sparse * u

function doit(D, text, du, u)
  println(text)
  sleep(0.1)
  show(stdout, MIME"text/plain"(), @benchmark mul!($du, $D, $u))
  println()
end

doit(D_SBP, "D_SBP:", du, u)
doit(D_DEO, "D_DEO:", du, u)
doit(D_sparse, "D_sparse:", du, u)
```


Bounded domains:
```@example
using BenchmarkTools
using LinearAlgebra, SparseArrays
using SummationByPartsOperators, BandedMatrices

BLAS.set_num_threads(1) # make sure that BLAS is serial to be fair

T = Float64
xmin, xmax = T(0), T(1)

D_SBP = derivative_operator(MattssonNordström2004(), derivative_order=1,
                            accuracy_order=6, xmin=xmin, xmax=xmax, N=10^3)
D_sparse = sparse(D_SBP)
D_banded = BandedMatrix(D_SBP)

u = randn(eltype(D_SBP), size(D_SBP, 1)); du = similar(u);
@show D_SBP * u ≈ D_sparse * u ≈ D_banded * u

function doit(D, text, du, u)
  println(text)
  sleep(0.1)
  show(stdout, MIME"text/plain"(), @benchmark mul!($du, $D, $u))
  println()
end

doit(D_SBP, "D_SBP:", du, u)
doit(D_sparse, "D_sparse:", du, u)
doit(D_banded, "D_banded:", du, u)
```


## Dissipation operators

```@example
using BenchmarkTools
using LinearAlgebra, SparseArrays
using SummationByPartsOperators, BandedMatrices

BLAS.set_num_threads(1) # make sure that BLAS is serial to be fair

T = Float64
xmin, xmax = T(0), T(1)

D_SBP = derivative_operator(MattssonNordström2004(), derivative_order=1,
                            accuracy_order=6, xmin=xmin, xmax=xmax, N=10^3)
Di_SBP  = dissipation_operator(MattssonSvärdNordström2004(), D_SBP)
Di_sparse = sparse(Di_SBP)
Di_banded = BandedMatrix(Di_SBP)
Di_full   = Matrix(Di_SBP)

u = randn(eltype(D_SBP), size(D_SBP, 1)); du = similar(u);
@show Di_SBP * u ≈ Di_sparse * u ≈ Di_banded * u ≈ Di_full * u

function doit(D, text, du, u)
  println(text)
  sleep(0.1)
  show(stdout, MIME"text/plain"(), @benchmark mul!($du, $D, $u))
  println()
end

doit(D_SBP, "D_SBP:", du, u)
doit(Di_SBP, "Di_SBP:", du, u)
doit(Di_sparse, "Di_sparse:", du, u)
doit(Di_banded, "Di_banded:", du, u)
doit(Di_full, "Di_full:", du, u)
```
