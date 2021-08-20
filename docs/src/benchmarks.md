# Benchmarks

Here are some simple benchmarks. Take them with a grain of salt since they run
on virtual machines in the cloud to generate the documentation automatically.


## First-derivative operators

#### Periodic domains

Let's set up some benchmark code.

```@example first-derivative-periodic
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
```

First, we benchmark the implementation from SummationByPartsOperators.jl.
```@example first-derivative-periodic
doit(D_SBP, "D_SBP:", du, u)
```

Next, we compare this to the runtime obtained using a sparse matrix representation
of the derivative operator. Depending on the hardware etc., this can be an order
of magnitude slower than the optimized implementation from SummationByPartsOperators.jl.
```@example first-derivative-periodic
doit(D_sparse, "D_sparse:", du, u)
```

Finally, we benchmark the implementation of the same derivative operator in
DiffEqOperators.jl.
```@example first-derivative-periodic
doit(D_DEO, "D_DEO:", du, u)
```


#### Bounded domains

We start again by setting up some benchmark code.

```@example first-derivative-bounded
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
```

First, we benchmark the implementation from SummationByPartsOperators.jl.
```@example first-derivative-bounded
doit(D_SBP, "D_SBP:", du, u)
```

Again, we compare this to a representation of the derivative operator as a
sparse matrix. No surprise - it is again much slower, as in periodic domains.
```@example first-derivative-bounded
doit(D_sparse, "D_sparse:", du, u)
```

FInally, we compare it to a representation as banded matrix. Disappointingly,
this is still much slower than the optimized implementation from
SummationByPartsOperators.jl.
```@example first-derivative-bounded
doit(D_banded, "D_banded:", du, u)
```


## Dissipation operators

We follow the same structure as before. At first, we set up some benchmark code.

```@example dissipation
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
```

At first, let us benchmark the derivative and dissipation operators implemented
in SummationByPartsOperators.jl.
```@example dissipation
doit(D_SBP, "D_SBP:", du, u)
doit(Di_SBP, "Di_SBP:", du, u)
```

Next, we compare the results to sparse matrix representations. It will not
come as a surprise that these are again much (around an order of magnitude)
slower.
```@example dissipation
doit(Di_sparse, "Di_sparse:", du, u)
doit(Di_banded, "Di_banded:", du, u)
```

Finally, let's benchmark the same computation if a full (dense) matrix is used
to represent the derivative operator. This is obviously a bad idea but 🤷
```@example dissipation
doit(Di_full, "Di_full:", du, u)
```


## Structure-of-Arrays (SoA) and Array-of-Structures (AoS)

[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
tries to provide efficient support of

- `StaticVector`s from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl)
- [StructArrays.jl](https://github.com/JuliaArrays/StructArrays.jl)

To demonstrate this, let us set up some benchmark code.

```@example soa-aos
using BenchmarkTools
using StaticArrays, StructArrays
using LinearAlgebra, SparseArrays
using SummationByPartsOperators, BandedMatrices

BLAS.set_num_threads(1) # make sure that BLAS is serial to be fair

struct Vec5{T} <: FieldVector{5,T}
  x1::T
  x2::T
  x3::T
  x4::T
  x5::T
end

# Apply `mul!` to each component of a plain array of structures one after another
function mul_aos!(du, D, u, args...)
  for i in 1:size(du, 1)
    mul!(view(du, i, :), D, view(u, i, :), args...)
  end
end

T = Float64
xmin, xmax = T(0), T(1)

D_SBP = derivative_operator(MattssonNordström2004(), derivative_order=1,
                            accuracy_order=4, xmin=xmin, xmax=xmax, N=101)
D_sparse = sparse(D_SBP)
D_full   = Matrix(D_SBP)
```

At first, we benchmark the application of the operators implemented in
SummationByPartsOperators.jl and their representations as sparse and dense
matrices in the scalar case. As before, the sparse matrix representation
is around an order of magnitude slower and the dense matrix representation
is far off.
```@example soa-aos
println("Scalar case")
u = randn(T, size(D_SBP, 1)); du = similar(u)
println("D_SBP")
show(stdout, MIME"text/plain"(), @benchmark mul!($du, $D_SBP, $u))
println("\nD_sparse")
show(stdout, MIME"text/plain"(), @benchmark mul!($du, $D_sparse, $u))
println("\nD_full")
show(stdout, MIME"text/plain"(), @benchmark mul!($du, $D_full, $u))
```

Next, we use a plain array of structures (AoS) in the form of a two-dimensional
array and our custom `mul_aos!` implementation that loops over each component,
using `mul!` on `view`s.
Here, the differences between the timings are less pronounced.
```@example soa-aos
println("Plain Array of Structures")
u_aos_plain = randn(T, 5, size(D_SBP, 1)); du_aos_plain = similar(u_aos_plain)
println("D_SBP")
show(stdout, MIME"text/plain"(), @benchmark mul_aos!($du_aos_plain, $D_SBP, $u_aos_plain))
println("\nD_sparse")
show(stdout, MIME"text/plain"(), @benchmark mul_aos!($du_aos_plain, $D_sparse, $u_aos_plain))
println("\nD_full")
show(stdout, MIME"text/plain"(), @benchmark mul_aos!($du_aos_plain, $D_full, $u_aos_plain))
```

Now, we use an array of structures (AoS) based on `reinterpret` and standard
`mul!`. This is much more efficient for the implementation in SummationByPartsOperators.jl.
In Julia v1.6, this is also more efficient for sparse matrices but less efficient
for dense matrices (compared to the plain AoS approach with `mul_aos!` above).
```@example soa-aos
println("Array of Structures (reinterpreted array)")
u_aos_r = reinterpret(reshape, Vec5{T}, u_aos_plain); du_aos_r = similar(u_aos_r)
@show D_SBP * u_aos_r ≈ D_sparse * u_aos_r ≈ D_full * u_aos_r
mul!(du_aos_r, D_SBP, u_aos_r)
@show reinterpret(reshape, T, du_aos_r) ≈ du_aos_plain
println("D_SBP")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_aos_r, $D_SBP, $u_aos_r))
println("\nD_sparse")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_aos_r, $D_sparse, $u_aos_r))
println("\nD_full")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_aos_r, $D_full, $u_aos_r))
```

Next, we still use an array of structures (AoS), but copy the data into a plain
`Array` instead of using the `reinterpret`ed versions. There is no significant
difference to the previous version in this case.
```@example soa-aos
println("Array of Structures")
u_aos = Array(u_aos_r); du_aos = similar(u_aos)
@show D_SBP * u_aos ≈ D_sparse * u_aos ≈ D_full * u_aos
mul!(du_aos, D_SBP, u_aos)
@show du_aos ≈ du_aos_r
println("D_SBP")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_aos, $D_SBP, $u_aos))
println("\nD_sparse")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_aos, $D_sparse, $u_aos))
println("\nD_full")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_aos, $D_full, $u_aos))
```

Finally, let's look at a structure of arrays (SoA). Interestingly, this is
slower than the array of structures we used above. On Julia v1.6, the sparse
matrix representation performs particularly bad in this case.
```@example soa-aos
println("Structure of Arrays")
u_soa = StructArray(u_aos); du_soa = similar(u_soa)
@show D_SBP * u_soa ≈ D_sparse * u_soa ≈ D_full * u_soa
mul!(du_soa, D_SBP, u_soa)
@show du_soa ≈ du_aos
println("D_SBP")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_soa, $D_SBP, $u_soa))
println("\nD_sparse")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_soa, $D_sparse, $u_soa))
println("\nD_full")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_soa, $D_full, $u_soa))
```
