# Linear advection equation

```math
\begin{aligned}
    \partial_t u(t,x) + \partial_x (a(x) u(t,x)) &= 0, && t \in (0,T), x \in (x_{min}, x_{max}), \\
    u(0,x) &= u_0(x), && x \in (x_{min}, x_{max}), \\
    \text{boundary conditions}, &&& x \in \partial (x_{min}, x_{max}).
\end{aligned}
```

The boundary conditions depend on the sign of the transport velocity ``a``
at the boundary. In particular, specifying a Dirichlet type boundary condition
is only allowed for inflow boundaries, e.g. ``a(x_{min}) > 0`` at ``x = x_{min}``.

```@example
using SummationByPartsOperators, OrdinaryDiffEq
using Plots, LaTeXStrings, Printf

# general parameters
xmin = -1.
xmax = +1.
tspan = (0., 8.0)
afunc(x) = one(x)
u0func(x) = sinpi(x)
# Dirichlet type boundary conditions; they are used only at inflow boundaries
left_bc(t) = t >= 3 ? sinpi(t) : zero(t)
right_bc(t) = zero(t)

# discretization parameters
interior_order = 4
N = 101
# whether a split form should be applied or not
split_form = Val(false)

# setup spatial semidiscretization
D = derivative_operator(MattssonSvärdShoeybi2008(), 1, interior_order, xmin, xmax, N)
# whether or not artificial dissipation should be applied: nothing, dissipation_operator(D)
Di = nothing
semidisc = VariableLinearAdvectionNonperiodicSemidiscretisation(D, Di, afunc, split_form, left_bc, right_bc)
ode = semidiscretise(u0func, semidisc, tspan)

# solve ode
sol = solve(ode, SSPRK104(), dt=D.Δx, adaptive=false,
            saveat=range(first(tspan), stop=last(tspan), length=200))

# visualise the result
plot(xguide=L"x", yguide=L"u")
plot!(evaluate_coefficients(sol[end], semidisc), label="")
savefig("example_linear_advection.png")
```

![](example_linear_advection.png)
