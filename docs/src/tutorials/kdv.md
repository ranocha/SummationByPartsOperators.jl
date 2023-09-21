# Korteweg-de Vries equation

Let's consider the Korteweg-de Vries (KdV) equation

```math
\begin{aligned}
    \partial_t u(t,x) + \partial_x \frac{u(t,x)^2}{2} + \partial_x^3 u(t,x) &= 0, && t \in (0,T), x \in (x_{min}, x_{max}), \\
    u(0,x) &= u_0(x), && x \in (x_{min}, x_{max}), \\
\end{aligned}
```

with periodic boundary conditions. The KdV equation has the quadratic invariant

```math
    J = \frac{1}{2} \int u(t,x)^2 \mathrm{d}x.
```

A classical trick to conserve this invariant is to use following split form

```math
    u_t + \frac{1}{3} (u^2)_x + \frac{1}{3} u u_x + \partial_x^3 u = 0.
```

Indeed, integration by parts with periodic boundary conditions yields

```math
\begin{aligned}
    \partial_t J
    &=
    \int u u_t
    =
    -\frac{1}{3} \int u (u^2)_x - \frac{1}{3} \int u^2 u_x - \int u \partial_x^3 u
    \\
    &=
    0 + \frac{1}{3} \int u_x u^2 - \frac{1}{3} \int u^2 u_x + 0
    =
    0.
\end{aligned}
```

## Basic example using finite difference SBP operators

Let's create an appropriate discretization of this equation step by step. At first,
we load packages that we will use in this example.

```@example kdv
using SummationByPartsOperators, OrdinaryDiffEq
using LaTeXStrings; using Plots: Plots, plot, plot!, savefig
```

Next, we specify the initial data as Julia function as well as the
spatial domain and the time span. Here, we use an analytic soliton solution
of the KdV equation for the initial data.

```@example kdv
# traveling wave solution (soliton)
get_xmin() = 0.0  # left boundary of the domain
get_xmax() = 80.0 # right boundary of the domain
get_c() = 2 / 3   # wave speed
function usol(t, x)
    xmin = get_xmin()
    xmax = get_xmax()
    μ = (xmax - xmin) / 2
    c = get_c()
    A = 3 * c
    K = sqrt(1/c - 1)
    x_t = mod(x - c*t - xmin, xmax - xmin) + xmin - μ

    A / cosh(sqrt(3*A) / 6 * x_t)^2
end

tspan = (0.0, (get_xmax() - get_xmin()) / (3 * get_c()) + 10 * (get_xmax() - get_xmin()) / get_c())
```

Next, we implement the semidiscretization using the interface of
[OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl)
which is part of [DifferentialEquations.jl](https://diffeq.sciml.ai/latest/).
For simplicity, we just use the out-of-place version here since we do not have
to worry about appropriate temporary buffers when using automatic differentiation
in implicit time integration methods.

```@example kdv
function kdv(u, parameters, t)
    D1 = parameters.D1
    D3 = parameters.D3

    # conservative semidiscretization using a split form
    return (-1 / 3) * (u .* (D1 * u) + D1 * (u.^2)) - D3 * u
end
```

Next, we choose first- and third-derivative SBP operators `D1, D3`, evaluate
the initial data on the grid, and set up the semidiscretization as an ODE problem.

```@example kdv
N = 128 # number of grid points
D1 = periodic_derivative_operator(derivative_order=1, accuracy_order=8,
                                  xmin=get_xmin(), xmax=get_xmax(), N=N)
D3 = periodic_derivative_operator(derivative_order=3, accuracy_order=8,
                                  xmin=get_xmin(), xmax=get_xmax(), N=N)
u0 = usol.(first(tspan), grid(D1))
parameters = (; D1, D3)

ode = ODEProblem(kdv, u0, tspan, parameters);
```

Finally, we can solve the ODE using a Rosenbrock method with adaptive time stepping.
We use such a linearly implicit time integration method since the third-order
derivative makes the system stiff.

```@example kdv
sol = solve(ode, Rodas5(), saveat=range(first(tspan), stop=last(tspan), length=200))

plot(xguide=L"x", yguide=L"u")
plot!(evaluate_coefficients(sol[1], D1), label=L"u_0")
plot!(evaluate_coefficients(sol[end], D1), label=L"u_\mathrm{numerical}")
savefig("example_kdv.png");
```

![](example_kdv.png)


## Advanced visualization

Let's create an animation of the numerical solution.

```julia
using Printf; using Plots: Animation, frame, gif

let anim = Animation()
    idx = 1
    x, u = evaluate_coefficients(sol[idx], D1)
    fig = plot(x, u, xguide=L"x", yguide=L"u", xlim=extrema(x), ylim=(-0.05, 2.05),
              label="", title=@sprintf("\$t = %6.2f \$", sol.t[idx]))
    for idx in 1:length(sol.t)
        fig[1] = x, sol.u[idx]
        plot!(title=@sprintf("\$t = %6.2f \$", sol.t[idx]))
        frame(anim)
    end
    gif(anim, "example_kdv.gif")
end
```

![example_kdv_animation](https://user-images.githubusercontent.com/12693098/186075685-a4f12cf0-df6d-486a-ba5c-04d7a7466743.gif)


## Package versions

These results were obtained using the following versions.

```@example kdv
using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(["SummationByPartsOperators", "OrdinaryDiffEq"],
           mode=PKGMODE_MANIFEST)
```
