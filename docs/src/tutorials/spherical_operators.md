# Spherical operators

Spherical symmetry reduces a problem on a ball or cylinder to the radial
interval ``[0, R]``. The reduced equations retain the metric weight ``r^p``:
``p = 1`` describes rotational/cylindrical symmetry and ``p = 2`` describes
spherical symmetry in three spatial dimensions. The two first-order spatial
actions are

```math
G\phi \approx \partial_r \phi,
\qquad
Dv \approx \frac{1}{r^p}\partial_r\left(r^p v\right).
```

A smooth scalar field ``\phi`` is even under reflection through the origin,
while the radial component ``v`` of a smooth vector field is odd. In
particular, regular odd fields vanish at ``r = 0``. The spherical operators are
constructed by extending data with the appropriate parity to a mirrored
Cartesian grid, applying a Cartesian SBP derivative, and restricting the
result back to the non-negative radial half-grid.

Both scalar and radial-vector fields are stored as ordinary vectors with one
entry per radial node. A collocated grid contains ``r = 0``; the first entry of
a regular radial-vector field is therefore zero. A staggered grid has all
nodes strictly away from the origin, so the parity condition is incorporated
by the operator construction and no value at ``r = 0`` is stored.

The gradient and divergence have separate scalar and vector norm matrices,
``S`` and ``V``. Together they satisfy

```math
S D + G^T V = B,
```

where ``B`` represents the weighted outer-boundary contribution. It is zero
except for ``B_{N,N} = R^p``; radial regularity removes the origin
contribution. For ``p = 0``, use [`derivative_operator`](@ref) to construct an
ordinary Cartesian SBP operator instead.

## Constructing and applying the operators

```@example spherical_operators
using LinearAlgebra, SummationByPartsOperators

operators = spherical_operators(VretinarisSchnetter2026();
                                accuracy_order = 2,
                                N = 40,
                                R = 1.0,
                                p = 2,
                                grid = :collocated)

r = grid(operators)
G = gradient_operator(operators)
D = divergence_operator(operators)
S = scalar_mass(operators)
V = vector_mass(operators)
B = boundary_matrix(operators)

S * D + G' * V ≈ B
```

`gradient_operator` and `divergence_operator` return matrix-like views. They
support matrix-vector multiplication, `mul!`, `Matrix`, and `sparse`.

```@example spherical_operators
scalar = exp.(-20 .* r .^ 2)
radial_vector = r .* exp.(-20 .* r .^ 2)
gradient = G * scalar
divergence = D * radial_vector;
```

## Evolving a spherical Gaussian pulse

We reproduce the scalar-wave setup of
[Vretinaris and Schnetter (2026)](https://doi.org/10.1103/hg1v-qxsj).
Sections III--V of the paper describe the operator construction, convergence
properties, and numerical evolutions in more detail; a preprint is also
available as [arXiv:2606.05155](https://arxiv.org/abs/2606.05155).

Writing ``\Pi = \partial_t\Phi`` and ``\Psi = \partial_r\Phi``, the spherical
wave equation becomes

```math
\partial_t \Pi = \frac{1}{r^2}\partial_r(r^2\Psi),
\qquad
\partial_t \Psi = \partial_r\Pi.
```

The exact regular solution used to initialize the evolution is

```math
\Phi(t,r) = \frac{g(t+r)-g(t-r)}{r},
\qquad
g(u) = A\exp\left[-\frac{(u-r_c)^2}{d^2}\right].
```

As in the paper, we take ``A=1``, ``d=2``, ``r_c=10``, ``p=2``, a nominal
domain ``[0,40]``, and evolve to ``t=25``. The pulse moves inward, passes
smoothly through the coordinate origin around ``t=10``, and then moves
outward. We use ``h=1/4`` and fourth-order operators here. The staggered grid
starts at ``h/2`` and ends at ``40-h/2``; the collocated grid includes both
``r=0`` and ``r=40``.

```@example spherical_operators
using OrdinaryDiffEqHighOrderRK: DP8, ODEProblem, solve
import Plots
using Plots: plot, plot!, savefig, @animate, gif

h = 1 / 4
nominal_radius = 40.0
number_of_intervals = Int(nominal_radius / h)

collocated = spherical_operators(VretinarisSchnetter2026();
                                   accuracy_order = 4,
                                   N = number_of_intervals,
                                   R = nominal_radius,
                                   p = 2,
                                   grid = :collocated)
staggered = spherical_operators(VretinarisSchnetter2026();
                                 accuracy_order = 4,
                                 N = number_of_intervals,
                                 R = nominal_radius - h / 2,
                                 p = 2,
                                 grid = :staggered)

amplitude = 1.0
width = 2.0
center = 10.0
gaussian(u) = amplitude * exp(-((u - center) / width)^2)
gaussian_prime(u) = -2 * (u - center) / width^2 * gaussian(u)
gaussian_second(u) =
    (4 * (u - center)^2 / width^4 - 2 / width^2) * gaussian(u)

function initial_fields(r)
    if iszero(r)
        return (Pi = 2 * gaussian_second(0.0), Psi = 0.0)
    end
    numerator = gaussian(r) - gaussian(-r)
    Pi = (gaussian_prime(r) - gaussian_prime(-r)) / r
    Psi = (gaussian_prime(r) + gaussian_prime(-r)) / r - numerator / r^2
    return (; Pi, Psi)
end

function spherical_wave!(dw, w, parameters, t)
    (; G, D, S, B, r, n) = parameters
    Pi = @view w[1:n]
    Psi = @view w[(n + 1):(2n)]
    dPi = @view dw[1:n]
    dPsi = @view dw[(n + 1):(2n)]

    mul!(dPi, D, Psi)
    mul!(dPsi, G, Pi)
    iszero(first(r)) && (dPsi[1] = 0) # odd regularity on a collocated grid

    # Reflecting Psi(t,R)=0 boundary condition in energy-stable SAT form.
    dPi[end] -= B[end, end] / S[end, end] * Psi[end]
    return nothing
end

function solve_pulse(operators)
    r = grid(operators)
    initial = initial_fields.(r)
    Pi0 = getproperty.(initial, :Pi)
    Psi0 = getproperty.(initial, :Psi)
    n = length(r)
    parameters = (; G = gradient_operator(operators),
                  D = divergence_operator(operators),
                  S = scalar_mass(operators),
                  B = boundary_matrix(operators), r, n)
    problem = ODEProblem(spherical_wave!, vcat(Pi0, Psi0),
                         (0.0, 25.0), parameters)
    solution = solve(problem, DP8(); dt = 0.5h, adaptive = false,
                     saveat = range(0.0, 25.0; length = 101))
    return (; operators, r, n, solution)
end

collocated_run = solve_pulse(collocated)
staggered_run = solve_pulse(staggered);
```

The fixed time step ``\Delta t=0.5h`` and eighth-order DP8 integrator also
follow the paper. The reflective SAT cancels the outer-boundary term in the
SBP identity. Thus the discrete energy

```math
E_h = \frac{1}{2}\left(\Pi^T S\Pi + \Psi^T V\Psi\right)
```

is conserved up to time-integration and floating-point errors.

The following snapshots show the radial energy density. They make the
inward-moving packet, its interaction with the origin, and the outgoing
packet visible for both grid arrangements.

```@example spherical_operators
function radial_energy_density(run, t)
    w = run.solution(t)
    Pi = @view w[1:run.n]
    Psi = @view w[(run.n + 1):(2run.n)]
    return @. 0.5 * run.r^2 * (Pi^2 + Psi^2)
end

snapshot_times = (0.0, 10.0, 25.0)
snapshot_plot = plot(layout = (2, 3), size = (960, 500),
                     link = :both, legend = false)
for (row, (name, run)) in enumerate((("collocated", collocated_run),
                                     ("staggered", staggered_run)))
    for (column, t) in enumerate(snapshot_times)
        panel = (row - 1) * length(snapshot_times) + column
        plot!(snapshot_plot[panel], run.r, radial_energy_density(run, t);
              xlims = (0, nominal_radius), ylims = (0, 0.55),
              xlabel = row == 2 ? "r" : "",
              ylabel = column == 1 ? "radial energy density" : "",
              title = "$name, t = $(Int(t))", linewidth = 2)
    end
end
savefig(snapshot_plot, "spherical_wave_snapshots.svg")
nothing # hide
```

![](spherical_wave_snapshots.svg)

The animation overlays the collocated and staggered results. It plots
``r\Pi``, the same rescaling used in the original implementation to remove
the dominant ``1/r`` amplitude variation.

```@example spherical_operators
animation = @animate for t in collocated_run.solution.t
    collocated_state = collocated_run.solution(t)
    staggered_state = staggered_run.solution(t)
    plot(collocated_run.r,
         collocated_run.r .* collocated_state[1:collocated_run.n];
         label = "collocated", linewidth = 2,
         xlims = (0, nominal_radius), ylims = (-0.8, 0.8),
         xlabel = "r", ylabel = "r Pi", title = "t = $(round(t; digits = 2))")
    plot!(staggered_run.r,
          staggered_run.r .* staggered_state[1:staggered_run.n];
          label = "staggered", linewidth = 2, linestyle = :dash)
end
gif(animation, "spherical_wave_pulse.gif"; fps = 20)
nothing # hide
```

![](spherical_wave_pulse.gif)

## Supported Cartesian sources

The exact spherical closure solve has been tested for both collocated and
staggered grids and for ``p = 1, 2`` with the following Cartesian source/order
combinations. Other combinations are rejected by the constructor.

| Cartesian source | Supported accuracy orders |
|:--|:--|
| `MattssonNordström2004()` | 2, 4, 6 |
| `MattssonSvärdNordström2004()` | 2, 4, 6 |
| `MattssonSvärdShoeybi2008()` | 2, 4, 6 |
| `Mattsson2012()` | 2, 4, 6 |
| `Mattsson2014()` | 2 |
| `MattssonAlmquistCarpenter2014Extended()` | 2 |
| `MattssonAlmquistCarpenter2014Optimal()` | 2 |
| `Mattsson2017(:minus)` | 4, 6 |
| `DienerDorbandSchnetterTiglio2007()` | 2, 4 |

The currently available spherical constructions use uniform radial grids and
support accuracy orders two, four, and six.
