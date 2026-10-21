# Spherical operators

Spherical symmetry reduces a radial problem to the interval ``[0, R]`` while
retaining the metric weight ``r^p``. The corresponding SBP discretization is a
pair of operators: the gradient maps an even scalar field to an odd radial
vector field, and the divergence maps that vector field back to an even scalar
field. These fields have separate norm matrices, so a spherical pair is not a
single ordinary derivative operator.

`p = 0` is ordinary Cartesian geometry, not a spherical construction. Use
[`derivative_operator`](@ref) for that case. The spherical constructor supports
the metric exponents `p = 1` and `p = 2`, both covered by the tests.

```@repl spherical_operators
using LinearAlgebra, SparseArrays, SummationByPartsOperators

operators = spherical_operators(VretinarisSchnetter2026();
                                accuracy_order = 4,
                                N = 20,
                                R = 1.0,
                                p = 2,
                                grid = :collocated)

r = grid(operators)
G = gradient_operator(operators)
D = divergence_operator(operators)
S = scalar_mass(operators)
V = vector_mass(operators)

S * D + G' * V ≈ boundary_matrix(operators)
```

`gradient_operator`, `odd_derivative_operator`, and `divergence_operator`
return matrix-like views. They support matrix-vector multiplication, `mul!`,
`Matrix`, and `sparse`, while keeping their field roles explicit.

```@repl spherical_operators
u = exp.(-20 .* r.^2)
du = similar(u)
mul!(du, G, u)

scalar_mass(operators) == S
vector_mass(operators) == V
```

Use `grid = :staggered` to place the radial nodes away from the origin. The
collocated grid includes the origin and applies the odd-field regularity
condition there. The currently available constructions support accuracy orders
four and six on uniform radial grids.
