# Basic interface

Here, we discuss the basic interface of
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl).
We assume you are already familiar with the concept of SBP operators
in general and the [introduction](@ref intro-introduction) describing
how to construct specific operators.


## Applying SBP operators

All SBP operators implement the general interface of matrix vector
multiplication in Julia. The most simple version is to just use `*`,
e.g.,

```@repl
using SummationByPartsOperators

D = derivative_operator(MattssonNordström2004(),
                        derivative_order = 1, accuracy_order = 2,
                        xmin = 0.0, xmax = 1.0, N = 9)

x = grid(D)

u = @. sin(pi * x)

D * u

@allocated D * u
```

As you can see above, calling `D * u` allocates a new vector for the
result. If you want to apply an SBP operator multiple times and need
good performance, you should consider using an in-place update instead.
Julia provides the function `mul!` for this purpose.

```@repl
using LinearAlgebra, InteractiveUtils

@doc mul!
```

To improve the performance, you can pre-allocate an output vector
and call the non-allocating function `mul!`.

```@repl
using SummationByPartsOperators

D = derivative_operator(MattssonNordström2004(),
                        derivative_order = 1, accuracy_order = 2,
                        xmin = 0.0, xmax = 1.0, N = 9)

x = grid(D)

u = @. sin(pi * x)

du = similar(u); mul!(du, D, u)

du ≈ D * u

@allocated mul!(du, D, u)
```

All operators provided by
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
implement this 3-argument version of `mul!`.
Most operators also implement the 5-argument version of `mul!` that
can be used to scale the output and add it to some multiple of the
result vector.

```@repl
using SummationByPartsOperators

D = derivative_operator(MattssonNordström2004(),
                        derivative_order = 1, accuracy_order = 2,
                        xmin = 0.0, xmax = 1.0, N = 9)

x = grid(D); u = @. sin(pi * x); du = similar(u); mul!(du, D, u);

mul!(du, D, u, 2) # equivalent to du .= 2 * D * u

du ≈ 2 * D * u

@allocated mul!(du, D, u, 2)

du_background = rand(length(du)); du .= du_background

mul!(du, D, u, 2, 3) # equivalent to du .= 2 * D * u + 3 * du

du ≈ 2 * D * u + 3 * du_background

@allocated mul!(du, D, u, 2, 3)
```


## Integration and the mass/norm matrix

SBP operators come with a mass matrix yielding a quadrature rule. In
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl),
all operators typically have diagonal mass/norm matrices.
You can access them via [`mass_matrix`](@ref), e.g.,

```@repl
using SummationByPartsOperators

D = derivative_operator(MattssonNordström2004(),
                        derivative_order = 1, accuracy_order = 2,
                        xmin = 0.0, xmax = 1.0, N = 9)

mass_matrix(D)

D = periodic_derivative_operator(derivative_order = 1,
                                 accuracy_order = 2,
                                 xmin = 0.0, xmax = 1.0,
                                 N = 8)

mass_matrix(D)
```

If you want to use the quadrature associated with a mass matrix,
you do not need to form it explicitly. Instead, it is recommended
to use the function [`integrate`](@ref), e.g.,

```@repl
using SummationByPartsOperators

D = derivative_operator(MattssonNordström2004(),
                        derivative_order = 1, accuracy_order = 2,
                        xmin = 0.0, xmax = 1.0, N = 9)

M = mass_matrix(D)

x = grid(D)

u = x.^2

integrate(u, D)

integrate(u, D) ≈ sum(M * u)

integrate(u, D) ≈ integrate(x -> x^2, x, D)
```

For example, you can proceed as follows to compute the error of the
SBP operator when computing a derivative as follows.



```@repl
using SummationByPartsOperators

D = derivative_operator(MattssonNordström2004(),
                        derivative_order = 1, accuracy_order = 2,
                        xmin = 0.0, xmax = 1.0, N = 9)

M = mass_matrix(D)

x = grid(D)

difference = D * x.^3 - 3 * x.^2

error_l2 = sqrt(integrate(abs2, difference, D))
```
