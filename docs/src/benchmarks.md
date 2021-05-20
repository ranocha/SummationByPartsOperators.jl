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

D_periodic_serial  = periodic_derivative_operator(der_order, acc_order, xmin, xmax, N+1, Val{:serial}())
D_periodic_threads = periodic_derivative_operator(der_order, acc_order, xmin, xmax, N+1, Val{:threads}())
D_nonperiodic_serial  = derivative_operator(source, der_order, acc_order, xmin, xmax, N, Val{:serial}())
D_nonperiodic_sparse  = sparse(D_nonperiodic_serial)
D_nonperiodic_banded  = BandedMatrix(D_nonperiodic_serial)
D_nonperiodic_threads = derivative_operator(source, der_order, acc_order, xmin, xmax, N, Val{:threads}())

Random.seed!(12345)
u = randn(T, N)
dest = similar(u)

println("D_periodic_serial:");     sleep(0.1); display(@benchmark mul!($dest, $D_periodic_serial, $u))
println("D_periodic_threads:");    sleep(0.1); display(@benchmark mul!($dest, $D_periodic_threads, $u))

println("D_nonperiodic_serial:");  sleep(0.1); display(@benchmark mul!($dest, $D_nonperiodic_serial, $u))
println("D_nonperiodic_threads:"); sleep(0.1); display(@benchmark mul!($dest, $D_nonperiodic_threads, $u))
println("D_nonperiodic_sparse:");  sleep(0.1); display(@benchmark mul!($dest, $D_nonperiodic_sparse, $u))
println("D_nonperiodic_banded:");  sleep(0.1); display(@benchmark mul!($dest, $D_nonperiodic_banded, $u))
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
D_threads = derivative_operator(source_D, 1, acc_order, xmin, xmax, N, Val{:threads}())

Di_serial  = dissipation_operator(source_Di, D_serial)
Di_sparse  = sparse(Di_serial)
Di_full    = Matrix(Di_serial)
Di_threads = dissipation_operator(source_Di, D_threads)

Random.seed!(12345)
u = randn(T, N)
dest = similar(u)

println("D_serial:");  sleep(0.1); display(@benchmark mul!($dest, $D_serial, $u))
println("Di_serial:"); sleep(0.1); display(@benchmark mul!($dest, $Di_serial, $u))
println("Di_sparse:"); sleep(0.1); display(@benchmark mul!($dest, $Di_sparse, $u))
println("Di_full:");   sleep(0.1); display(@benchmark mul!($dest, $Di_full, $u))

println("D_threads:");  sleep(0.1); display(@benchmark mul!($dest, $D_threads, $u))
println("Di_threads:"); sleep(0.1); display(@benchmark mul!($dest, $Di_threads, $u))
```
