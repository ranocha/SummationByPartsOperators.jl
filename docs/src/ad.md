# Automatic/algorithmic differentiation (AD)

[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
is written using generic Julia code whenever possible. This means that
standard AD tools just work. For example, computing the Jacobian using
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) is possible
using the following code.

```jldoctest
julia> using SummationByPartsOperators, ForwardDiff

julia> D = periodic_derivative_operator(derivative_order = 1, accuracy_order = 2,
                                        xmin = 0.0, xmax = 1.0, N = 8)
Periodic first-derivative operator of order 2 on a grid in [0.0, 1.0] using 8 nodes,
stencils with 1 nodes to the left, 1 nodes to the right, and coefficients of Fornberg (1998)
  Calculation of Weights in Finite Difference Formulas.
  SIAM Rev. 40.3, pp. 685-691.

julia> u = rand(size(D, 2));

julia> J = ForwardDiff.jacobian(u -> D * u, u)
8×8 Matrix{Float64}:
  0.0   4.0   0.0   0.0   0.0   0.0   0.0  -4.0
 -4.0   0.0   4.0   0.0   0.0   0.0   0.0   0.0
  0.0  -4.0   0.0   4.0   0.0   0.0   0.0   0.0
  0.0   0.0  -4.0   0.0   4.0   0.0   0.0   0.0
  0.0   0.0   0.0  -4.0   0.0   4.0   0.0   0.0
  0.0   0.0   0.0   0.0  -4.0   0.0   4.0   0.0
  0.0   0.0   0.0   0.0   0.0  -4.0   0.0   4.0
  4.0   0.0   0.0   0.0   0.0   0.0  -4.0   0.0

julia> J ≈ Matrix(D)
true
```

This works of course also for non-periodic SBP operators, e.g.,

```jldoctest
julia> using SummationByPartsOperators, ForwardDiff

julia> D = derivative_operator(MattssonNordström2004(),
                               derivative_order = 1, accuracy_order = 2,
                               xmin = 0.0, xmax = 1.0, N = 8)
SBP first-derivative operator of order 2 on a grid in [0.0, 1.0] using 8 nodes
and coefficients of Mattsson, Nordström (2004)
  Summation by parts operators for finite difference approximations of second
    derivatives.
  Journal of Computational Physics 199, pp. 503-540.

julia> u = rand(size(D, 2));

julia> J = ForwardDiff.jacobian(u -> D * u, u)
8×8 Matrix{Float64}:
 -7.0   7.0   0.0   0.0   0.0   0.0   0.0  0.0
 -3.5   0.0   3.5   0.0   0.0   0.0   0.0  0.0
  0.0  -3.5   0.0   3.5   0.0   0.0   0.0  0.0
  0.0   0.0  -3.5   0.0   3.5   0.0   0.0  0.0
  0.0   0.0   0.0  -3.5   0.0   3.5   0.0  0.0
  0.0   0.0   0.0   0.0  -3.5   0.0   3.5  0.0
  0.0   0.0   0.0   0.0   0.0  -3.5   0.0  3.5
  0.0   0.0   0.0   0.0   0.0   0.0  -7.0  7.0

julia> J ≈ Matrix(D)
true
```

However, this does not work for Fourier derivative operators - and all other
operators involving an FFT - since [FFTW.jl](https://github.com/JuliaMath/FFTW.jl)
cannot handle dual numbers and a simple `reinterpret` trick does also not help
since [FFTW.jl](https://github.com/JuliaMath/FFTW.jl) requires unit strides.


## Jacobian-vector products

There is a nat trick that you can use if you are only interested in Jacobian-vector
products. [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) does not offer
such a functionality at the time of writing, but Simon Byrne suggested the following
implementation.

```@example ForwardDiff_StructArrays
using SummationByPartsOperators, ForwardDiff, StructArrays

function StructDual(x::AbstractVector{T}, w::AbstractVector{T}) where {T}
  @assert length(x) == length(w)
  # This was the original suggestion. However, it is currently not stable
  # under broadcasting. Thus, we use a slightly different version.
  # partials = StructArray{ForwardDiff.Partials{1, T}}(
  #     (StructArray{Tuple{T}}(
  #         (w,)
  #     ),)
  # )
  partials = reinterpret(reshape, ForwardDiff.Partials{1, T}, w)
  duals = StructArray{ForwardDiff.Dual{Nothing, T, 1}}((x, partials))
  return duals
end

function ForwardDiff.value(dx::StructArray{D}) where {D <: ForwardDiff.Dual}
  return dx.value
end

function ForwardDiff.partials(dx::StructArray{<: ForwardDiff.Dual{Tag, T, 1}}, i) where {Tag, T}
  # This was the original suggestion. We need to update it (see above).
  # return getproperty(dx.partials.values, i)
  @assert i == 1
  return reinterpret(reshape, T, dx.partials)
end
```

You can use it as follows to compute the Jacobian-vector product

```math
J_f(u) \cdot v
```

for the function ``f`` given by

```math
f(u) = D u.
```

```@example ForwardDiff_StructArrays
D = fourier_derivative_operator(xmin = 0.0, xmax = 1.0, N = 8)

u = randn(size(D, 2)); # the point `u` where we want to compute the derivative
v = randn(size(D, 2)); # the direction `v` in which we want to compute the derivative
u_v = StructDual(u, v); # combined StructArray containing the value and direction

f_df = D * u_v # compute the function value and its derivative

@assert ForwardDiff.value(f_df) ≈ D * u

@assert ForwardDiff.partials(f_df, 1) ≈ D * v # the Jacobian of `f(u) = D * u` is `D`
```

You can of course also use this with nonlinear functions, e.g.,

```@example ForwardDiff_StructArrays
f(u, D) = u .* (D * (u.^2))

f_df = f(u_v, D)
```

The Jacobian of this function is

```@example ForwardDiff_StructArrays
using LinearAlgebra

J = Diagonal(D * u.^2) + 2 .* u .* Matrix(D) * Diagonal(u)

@assert ForwardDiff.value(f_df) ≈ f(u, D)

@assert ForwardDiff.partials(f_df, 1) ≈ J * v
```


## Reproducibility

These results were obtained using the following versions.
```@example
using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(["SummationByPartsOperators", "ForwardDiff", "StructArrays"],
           mode=PKGMODE_MANIFEST)
nothing # hide
```
