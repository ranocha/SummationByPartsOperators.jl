# SummationByPartsOperators.jl

The Julia library
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
provides a unified interface of different discretization approaches including
finite difference, Fourier pseudospectral, continuous Galerkin, and discontinuous
Galerkin methods.
This unified interface is based on the notion of summation-by-parts (SBP)
operators. Originally developed for finite difference methods, SBP operators
are discrete derivative operators designed specifically to get provably stable
(semi-) discretizations, mimicking energy/entropy estimates from the continuous
level discretely and paying special attention to boundary conditions.

SummationByPartsOperators.jl is mainly written to be useful for both students
learning the basic concepts and researchers developing new numerical algorithms
based on SBP operators. Thus, this package uses Julia's multiple dispatch and
strong type system to provide a unified framework of all of these seemingly
different discretizations while being reasonably optimized at the same time,
achieving good performance without sacrificing flexibility.


## Installation

[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
is a registered Julia package. Thus, you can install it from the Julia REPL via
```julia
julia> using Pkg; Pkg.add("SummationByPartsOperators")
```

If you want to update SummationByPartsOperators.jl, you can use
```julia
julia> using Pkg; Pkg.update("SummationByPartsOperators")
```
As usual, if you want to update SummationByPartsOperators.jl and all other
packages in your current project, you can execute
```julia
julia> using Pkg; Pkg.update()
```


## Basic examples

Compute the derivative on a periodic domain using a central finite difference operator.
```julia
julia> using SummationByPartsOperators

julia> using Plots: plot, plot!

julia> D = periodic_derivative_operator(derivative_order=1, accuracy_order=2,
                                        xmin=0.0, xmax=2.0, N=20)
Periodic first-derivative operator of order 2 on a grid in [0.0, 2.0] using 20 nodes,
stencils with 1 nodes to the left, 1 nodes to the right, and coefficients of Fornberg (1998)
  Calculation of Weights in Finite Difference Formulas.
  SIAM Rev. 40.3, pp. 685-691.

julia> x = grid(D); u = sinpi.(x);

julia> plot(x, D * u, label="numerical")

julia> plot!(x, π .* cospi.(x), label="analytical")
```
You should see a plot like the following.

![](https://user-images.githubusercontent.com/12693098/118977199-2ef4b280-b976-11eb-8e02-aec722d75bfa.png)


Compute the derivative on a bounded domain using an SBP finite difference operator.
```julia
julia> using SummationByPartsOperators

julia> using Plots: plot, plot!

julia> D = derivative_operator(MattssonNordström2004(), derivative_order=1, accuracy_order=2,
                               xmin=0.0, xmax=1.0, N=21)
SBP first-derivative operator of order 2 on a grid in [0.0, 1.0] using 21 nodes
and coefficients of Mattsson, Nordström (2004)
  Summation by parts operators for finite difference approximations of second
    derivatives.
  Journal of Computational Physics 199, pp. 503-540.

julia> x = grid(D); u = exp.(x);

julia> plot(x, D * u, label="numerical")

julia> plot!(x, exp.(x), label="analytical")
```
You should see a plot like the following.

![](https://user-images.githubusercontent.com/12693098/118978404-93fcd800-b977-11eb-80b3-3dbfce5ecfd6.png)


## Referencing

If you use
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
for your research, please cite it using the bibtex entry
```bibtex
@article{ranocha2021sbp,
  title={{SummationByPartsOperators.jl}: {A} {J}ulia library of provably stable
         semidiscretization techniques with mimetic properties},
  author={Ranocha, Hendrik},
  journal={Journal of Open Source Software},
  year={2021},
  month={08},
  doi={10.21105/joss.03454},
  volume={6},
  number={64},
  pages={3454},
  publisher={The Open Journal},
  url={https://github.com/ranocha/SummationByPartsOperators.jl}
}
```
Please also cite the appropriate references for specific SBP operators
you use, which can be obtained via [`source_of_coefficients`](@ref).


## License and contributing

This project is licensed under the MIT license (see [License](@ref)).
Since it is an open-source project, we are very happy to accept contributions
from the community. Please refer to the section [Contributing](@ref) for more
details.
