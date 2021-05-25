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
date: 25 May 2021 # TODO: date
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
be solved numerically, a vast amount of numerical methods has been developed,
resulting sometimes in completely disjoint communities. To transfer developments
and knowledge from one community to another, common abstractions can be helpful.
The concept of summation-by-parts (SBP) operators is such an abstraction.
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
the implementations are reasonably fast. Together, this facilitates the development
of new algorithms and research in numerical analysis, which is the primary goal
of this package. In addition, SummationByPartsOperators.jl has been used in a
number of graduate-level numerical analysis courses.


# Features

SummationByPartsOperators.jl implements numerical methods based on SBP operators
of the following classes:

- finite difference methods
- Fourier collocation methods
- continuous Galerkin methods
- discontinuous Galerkin methods

Since a discrete derivative operator is basically a linear operator, all of
these SBP operators implement the basic interface of such linear operators
(`AbstractMatrix` in Julia) such as multiplication by vectors and addition of
operators. Finite difference and Fourier operators on periodic domains also
allow the construction of rational functions of operators and their efficient
implementation using the fast Fourier transform [@frigo2005design].

In addition to basic SBP derivative operators, SummationByPartsOperators.jl
contains a number of artificial dissipation and filtering operators, such as

- SBP artificial dissipation operators
- spectral viscosity operators for Fourier methods
- modal filter operators for Fourier methods and Legendre pseudospectral methods

Following good software development practices, SummationByPartsOperators.jl
makes use of continuous integration and automated tests required before merging
pull requests. Documentation is provided both in form of docstrings, a general
introduction, and tutorials.


# Related research and software

There are of course many open-source software packages providing discretizations
of differential equations. However, many of them focus on a single class of
numerical methods, e.g.,

- finite difference methods ([DiffEqOperators.jl](https://github.com/SciML/DiffEqOperators.jl), a part of DifferentialEquations.jl [@rackauckas2017differentialequations])
- finite volume methods [@ramadhan2020oceananigans;@xiao2021kinetic]
- spectral methods [@olver2014practical;@constantinou2021fourierflows]
- finite element methods [@badia2020gridap]
- discontinuous spectral element methods [@schlottkelakemper2021purely;@schlottkelakemper2020trixi]

We are not aware of any open-source software library implementing all of the
SBP classes using a unified interface or even several finite difference (FD)
SBP operators on finite domains. Such FD SBP operators are usually heavily
optimized [@mattsson2014optimal;@mattsson2018boundary] and not available in
other open source packages. Sometimes, restricted sets of coefficients are
available online [@almquist2017optimized;@oreilly2019sbp], but there is no other
extensive collection of these methods.

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
