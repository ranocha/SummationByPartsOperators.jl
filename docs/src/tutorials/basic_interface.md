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

