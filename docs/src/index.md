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
