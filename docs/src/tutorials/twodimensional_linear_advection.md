# Linear advection equation in two dimensions

In this tutorial, we consider the linear advection equation in two dimensions

```math
\begin{aligned}
    \partial_t u(t,x,y) + a_1\partial_x u(t,x,y) + a_2\partial_y u(t,x,y) &= 0, && t \in (0,T), (x, y)\in\Omega \\
    u(0,x,y) &= u_0(x,y), && (x,y) \in\Omega, \\
    u(t,x,y) &= g(t), && t \in (0,T), (x,y) \in\Gamma_-,
\end{aligned}
```

Here, ``\Omega`` is a two-dimensional domain, and ``\Gamma_- = \{(x, y)\in\partial\Omega|\boldsymbol{a}\cdot\boldsymbol{n}(\boldsymbol{x})\}`` is the inflow boundary,
where ``\boldsymbol{n}(\boldsymbol{x})`` is the outwards pointing normal vector at ``\boldsymbol{x}$ and $\boldsymbol{a} = [a_1, a_2]^T``.
We will use a rectangular domain ``\Omega = [x_{\min}, x_{\max}]\times[y_{\min}, y_{\max}]`` and a tensor-product summation-by-parts operator, which is based on
one-dimensional SBP operators in each direction. Based on one-dimensional SBP operators ``D_1`` on ``N_x`` nodes in ``[x_{\min}, x_{\max}]`` and ``D_2`` on ``N_y`` nodes
in ``[y_{\min}, y_{\max}]``, we can construct a two-dimensional SBP operator ``D`` on ``N_x\cdot N_y`` nodes utilizing Kronecker products.

```@example twodimensional_advection
using SummationByPartsOperators, OrdinaryDiffEq
using LaTeXStrings; using Plots: Plots, scatter, savefig

# Define the domain
xmin, xmax = -1.0, 1.0
ymin, ymax = -2.0, 2.0
N_x, N_y = 30, 60

# Construct one-dimensional SBP operators
D_1 = derivative_operator(MattssonNordström2004(), derivative_order = 1, accuracy_order = 4,
                          xmin = xmin, xmax = xmax, N = N_x)
D_2 = derivative_operator(MattssonNordström2004(), derivative_order = 1, accuracy_order = 4,
                          xmin = ymin, xmax = ymax, N = N_y)

# Construct the two-dimensional SBP operator
D = tensor_product_operator_2D(D_1, D_2)
```

The nodes are stored as a `Vector` of `SVector`s. We can visualize the grid as follows:

```@example twodimensional_advection
nodes = grid(D)
scatter(first.(nodes), last.(nodes), label = "Grid nodes", xlabel = L"x", ylabel = L"y")
savefig("two_dimensional_grid.png");
```

![](two_dimensional_grid.png)

A multi-dimensional SBP operators stores differentiation matrices in each direction. In the case of
a tensor-product operator, these matrices are stored as sparse matrices. They can be accessed by indexing the operator:

```@example twodimensional_advection
D_x = D[1]
```

```@example twodimensional_advection
D_y = D[2]
```

As in the one-dimensional case, a multi-dimensional SBP operator, stores the weights of quadrature rule, which performs numerical
integration in the interior of the domain. In contrast to the one-dimensional case, we also need to compute non-trivial integrals
along the boundary, which
