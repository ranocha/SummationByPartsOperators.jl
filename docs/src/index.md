# SummationByPartsOperators.jl

[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
is a Julia library of summation-by-parts (SBP) operators, which are discrete
derivative operators developed to get provably stable semidiscretizations,
paying special attention to boundary conditions. Discretizations included in this
framework are finite difference, Fourier pseudospectral, continuous Galerkin,
and discontinuous Galerkin methods. The main aim of
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
is to be useful for researchers and students to learn the basic concepts by
providing a unified framework of all of these seemingly different discretizations.
At the same time, the implementation is optimized to achieve good performance
without sacrificing flexibility.


## Installation

[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
is a registered Julia package. Thus, you can install it from the Julia REPL via
```julia
julia> using Pkg; Pkg.add("SummationByPartsOperators")
```


## Basic examples

Compute the derivative on a periodic domain using a central finite difference operator.
```julia
julia> using SummationByPartsOperators

julia> using Plots: plot, plot!

julia> D = periodic_derivative_operator(derivative_order=1, accuracy_order=2,
                                        xmin=0.0, xmax=2.0, N=21)
Periodic 1st derivative operator of order 2 {T=Float64, Parallel=Val{:serial}}
on a grid in [0.0, 2.0] using 21 nodes,
stencils with 1 nodes to the left, 1 nodes to the right, and coefficients from
  Fornberg (1998)
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
SBP 1st derivative operator of order 2 {T=Float64, Parallel=Val{:serial}}
on a grid in [0.0, 1.0] using 21 nodes
and coefficients given in
  Mattsson, Nordström (2004)
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
@misc{ranocha2021sbp,
  title={{SummationByPartsOperators.jl}: {A} {J}ulia library of provably stable
         semidiscretization techniques with mimetic properties},
  author={Ranocha, Hendrik},
  year={2021},
  howpublished={\url{https://github.com/ranocha/SummationByPartsOperators.jl},
  doi={10.5281/zenodo.4773575}
}
```


## License and contributing

This project is licensed under the MIT license (see [License](@ref)).
Since it is an open-source project, we are very happy to accept contributions
from the community. Please refer to the section [Contributing](@ref) for more
details.
