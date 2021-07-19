---
title: 'SummationByPartsOperators.jl: A Julia library of provably stable discretization techniques with mimetic properties'
tags:
  - Julia
  - numerical analysis
  - differential equations
  - summation-by-parts
  - energy stability
  - finite differences
  - discontinuous Galerkin methods
  - Fourier collocation methods
authors:
  - name: Hendrik Ranocha
    orcid: 0000-0002-3456-2277
    affiliation: 1
affiliations:
 - name: Applied Mathematics Münster, University of Münster, Germany
   index: 1
date: 19 July 2021
bibliography: paper.bib
---


# Summary

[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
is a Julia library of summation-by-parts (SBP) operators, which are discrete
derivative operators developed to get provably stable (semi-) discretizations,
paying special attention to boundary conditions. Discretizations included in this
framework are finite difference, Fourier pseudospectral, continuous Galerkin,
and discontinuous Galerkin methods.
The main aim of SummationByPartsOperators.jl is to be useful for both students
learning the basic concepts and researchers developing new numerical algorithms based
on SBP operators. Therefore, SummationByPartsOperators.jl provides a unified
framework of all of these seemingly different discretizations. At the same time,
the implementation is reasonably optimized to achieve good performance without
sacrificing flexibility.


# Statement of need

Partial differential equations (PDEs) are widely used in science and engineering
to create mathematical models of real-world processes. Since PDEs often need to
be solved numerically, a vast amount of numerical methods has been developed.
Since it is impossible to keep up with all recent research, sub-communities
focussing on specific applications or numerical methods emerged. Sometimes,
these communities develop different vocabulary and notations, making it hard for
newcomers (or even experienced researchers) to see similarities and connections
between seemingly unrelated approaches. To transfer new ideas and developments
and knowledge from one community to another, common abstractions can be helpful.
The concept of SBP operators is such an abstraction.
In recent years, SBP operators have attracted a lot of attention, in particular
for PDEs modeling advection-dominated problems, where they enabled the construction
of energy- or entropy-stable numerical methods, including finite difference,
discontinuous Galerkin, continuous Galerkin, and (pseudo-) spectral methods.
Their success is based on mimetic properties which enable the transfer of
results obtained for differential equations at the continuous level to the
discrete level. In particular, SBP operators are designed to mimic integration-by-parts
discretely as summation-by-parts, enabling discrete analogues of energy/entropy
methods for PDEs.

SummationByPartsOperators.jl is written entirely in Julia [@bezanson2017julia].
Making use of multiple dispatch and generic types, SummationByPartsOperators.jl
provides a unified interface for different SBP operators. At the same time,
the implementations are reasonably fast (again, due to multiple dispatch and specialized
implementations for each operator class). Together, this facilitates the development
of new algorithms and research in numerical analysis, which is the primary goal
of this package. In addition, SummationByPartsOperators.jl has been used in a
number of graduate-level numerical analysis courses, allowing students to
understand the connections between different SBP methods by presenting them in
a unified framework. In addition, some of the operators were not available in
open source software previously (to the best of the author's knowledge).


# Features

SummationByPartsOperators.jl implements numerical methods based on SBP operators
of the following classes:

- finite difference methods
- Fourier collocation methods
- continuous Galerkin methods
- discontinuous Galerkin methods

Since a discrete derivative operator is a linear operator, all of
these SBP operators implement the basic interface of such linear operators
(`AbstractMatrix` in Julia) such as multiplication by vectors and addition of
operators. Finite difference and Fourier operators on periodic domains also
allow the construction of rational functions of operators and their efficient
implementation using the fast Fourier transform [@frigo2005design].

In addition to basic SBP derivative operators, SummationByPartsOperators.jl
contains a number of related operators, such as

- SBP artificial dissipation operators
- spectral viscosity operators for Fourier methods
- modal filter operators for Fourier methods and Legendre pseudospectral methods

Using Julia's rich type system, all of these operators are implemented as their
own types. This enables several optimizations such as a memory requirement
independent of the number of grid points. In contrast, implementations based
on sparse/banded matrices have a memory requirement growing linearly with the
number of grid points. In addition, the structure of the operators can be taken
into account for operator-vector multiplications, usually resulting in speed-ups
of an order of magnitude or more on consumer hardware. For example, the application
an optimized the sixth-order (in the interior) finite difference SBP operator
[@almquist2017optimized] on a grid with 1000 nodes takes roughly 330 ns
on a consumer CPU from 2017 (Intel® Core™ i7-8700K) using version v0.5.5 of
SummationByPartsOperators.jl. In contrast, the same operation takes roughly
3.9 microseconds using a sparse matrix format used in other implementations of
this operator [@almquist2017optimized]. This benchmark is based on the following
code, which also provides a very basic example of SummationByPartsOperators.jl.
```julia
julia> # install the package if necessary
       using Pkg; Pkg.add("SummationByPartsOperators")

julia> using SummationByPartsOperators # load the package

julia> # create a finite difference SBP operator
       D = derivative_operator(MattssonAlmquistVanDerWeide2018Accurate(),
          derivative_order=1, accuracy_order=6, xmin=0.0, xmax=1.0, N=10^3)
SBP first-derivative operator of order 6 on a grid in [0.0, 1.0] using 1000 nodes
and coefficients of Mattsson, Almquist, van der Weide (2018)
  Boundary optimized diagonal-norm SBP operators ('Accurate').
  Journal of Computational Physics 374, pp. 1261-1266.

julia> # evaluate the function `sinpi` on the discrete grid
       x = grid(D); u = sinpi.(x);

julia> du = D * u; # use `D` to approximate the derivative of `u`

julia> # compute the discrete L² error of the approximation
       integrate(u -> u^2, du - pi * cospi.(x), D) |> sqrt
4.238102975456189e-13
```

Following good software development practices, SummationByPartsOperators.jl
makes use of continuous integration and automated tests required before merging
pull requests. Documentation is provided in form of docstrings, a general
introduction, and tutorials. In addition, SummationByPartsOperators.jl is a
registered Julia package and can be installed using the built-in package manager,
handling dependencies and version requirements automatically.


# Related research and software

There are of course many open-source software packages providing discretizations
of differential equations. However, many of them focus on a single class of
numerical methods or a specific application, e.g.,

- finite difference methods ([DiffEqOperators.jl](https://github.com/SciML/DiffEqOperators.jl),
  a part of DifferentialEquations.jl [@rackauckas2017differentialequations])
- finite volume methods ([Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) [@ramadhan2020oceananigans],
  [Kinetic.jl](https://github.com/vavrines/Kinetic.jl) [@xiao2021kinetic])
- spectral methods ([ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl) [@olver2014practical],
  [FourierFlows.jl](https://github.com/FourierFlows/FourierFlows.jl) [@constantinou2021fourierflows])
- finite element methods ([Gridap.jl](https://github.com/gridap/Gridap.jl) [@badia2020gridap])
- discontinuous spectral element methods
  ([Trixi.jl](https://github.com/trixi-framework/Trixi.jl) [@schlottkelakemper2021purely;@schlottkelakemper2020trixi])

We are not aware of any open-source software library implementing all of the
SBP classes using a unified interface or even several finite difference
SBP operators on finite domains, which are usually heavily
optimized [@mattsson2014optimal;@mattsson2018boundary] and not available in
other open source packages. Sometimes, restricted sets of coefficients are
available online [@almquist2017optimized;@oreilly2019sbp], but there is no other
extensive collection of these methods.

Of course, there is a plethora of additional open source software implementing
numerical methods for PDEs and each package has its own design criteria and goals.
SummationByPartsOperators.jl provides a unified interface of different SBP
operators. Thus, there is a partial overlap with some of the packages mentioned
above such as finite difference operators on periodic domains (DiffEqOperators.jl,
but with a different handling of bounded domains) or Fourier methods on periodic
domains (ApproxFun.jl, but with a different interface and extensions). In addition,
many packages focus on a specific application such as some specific fluid models
(Oceananigans.jl, FourierFlows.jl) or hyperbolic PDEs (Trixi.jl). In contrast,
SummationByPartsOperators.jl focuses on the numerical methods and provides them
in a form usable for rather general PDEs. For example, there is ongoing work to
use the basic operators provided by SummationByPartsOperators.jl in Trixi.jl.

Some of the research projects that have made use of SummationByPartsOperators.jl
(most of which have led to its further development) include numerical analysis
of and algorithms for

- nonlinear dispersive wave equations
  [@ranocha2021broad;@ranocha2021rate]
- hyperbolic conservation laws
  [@offner2019error;@lefloch2021kinetic;@ranocha2020preventing]
- ordinary differential equations
  [@ranocha2021class;@ranocha2020energy;@ranocha2021strong]
- Helmholtz Hodge decomposition and analysis of plasma waves
  [@ranocha2020discrete]


# Acknowledgements

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation)
under Grant SO~363/14-1 and Germany's Excellence Strategy EXC 2044-390685587,
Mathematics Münster: Dynamics-Geometry-Structure as well as King Abdullah
University of Science and Technology (KAUST).


# References
