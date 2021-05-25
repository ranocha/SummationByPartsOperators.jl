
"""
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

Check out the [documentation](https://ranocha.github.io/SummationByPartsOperators.jl/stable)
for further information. Some noticable functions to start with are
[`derivative_operator`](@ref),
[`legendre_derivative_operator`](@ref),
[`periodic_derivative_operator`](@ref),
[`fourier_derivative_operator`](@ref),
[`dissipation_operator`](@ref),
and [`grid`](@ref).

If you use this package for your research, please cite it using
```bibtex
@misc{ranocha2021sbp,
  title={{SummationByPartsOperators.jl}: {A} {J}ulia library of provably stable
         semidiscretization techniques with mimetic properties},
  author={Ranocha, Hendrik},
  year={2021},
  howpublished={\\url{https://github.com/ranocha/SummationByPartsOperators.jl},
  doi={10.5281/zenodo.4773575}
}
```
"""
module SummationByPartsOperators

using LinearAlgebra
using SparseArrays

using ArgCheck: @argcheck
using FFTW
using LoopVectorization: @turbo, @tturbo
using Reexport: @reexport
using Requires
using StaticArrays
using UnPack: @unpack
using Unrolled

@reexport using DiffEqBase
using DiffEqCallbacks

import LinearAlgebra: mul!
@reexport using PolynomialBases
import PolynomialBases: integrate, evaluate_coefficients, evaluate_coefficients!,
                        compute_coefficients, compute_coefficients!


# types
abstract type AbstractDerivativeOperator{T} end
abstract type AbstractNonperiodicDerivativeOperator{T} <: AbstractDerivativeOperator{T} end
abstract type AbstractPeriodicDerivativeOperator{T} <: AbstractDerivativeOperator{T} end
abstract type AbstractDerivativeCoefficients{T} end
abstract type AbstractMassMatrix{T} end
abstract type AbstractSemidiscretization end
abstract type AbstractFilter{T<:Real} end
abstract type AbstractFilterFunction end
"""
    SourceOfCoefficients

All sources of coefficients (articles) are subtypes of this abstract type.
"""
abstract type SourceOfCoefficients end


# source files
include("general_operators.jl")
include("periodic_operators.jl")
include("SBP_operators.jl")
include("dissipation_operators.jl")
include("var_coef_operators.jl")
include("coupling.jl")

function __init__()
  @require BandedMatrices="aae01518-5342-5314-be14-df237901396f" include("banded_matrices.jl")
end

include("filter.jl")
include("fourier_operators.jl")
include("fourier_operators_2d.jl")
include("legendre_operators.jl")
include("SBP_coefficients/MattssonNordström2004.jl")
include("SBP_coefficients/MattssonSvärdNordström2004.jl")
include("SBP_coefficients/MattssonSvärdShoeybi2008.jl")
include("SBP_coefficients/Mattsson2012.jl")
include("SBP_coefficients/Mattsson2014.jl")
include("SBP_coefficients/MattssonAlmquistCarpenter2014Extended.jl")
include("SBP_coefficients/MattssonAlmquistCarpenter2014Optimal.jl")
include("SBP_coefficients/Mattsson2017.jl")
include("SBP_coefficients/MattssonAlmquistVanDerWeide2018Minimal.jl")
include("SBP_coefficients/MattssonAlmquistVanDerWeide2018Accurate.jl")

include("conservation_laws/general_laws.jl")
include("conservation_laws/burgers.jl")
include("conservation_laws/cubic.jl")
include("conservation_laws/variable_linear_advection.jl")
include("conservation_laws/quartic_nonconvex.jl")
include("second_order_eqs/wave_eq.jl")

# exports
export PeriodicDerivativeOperator, PeriodicDissipationOperator,
       PeriodicRationalDerivativeOperator, PeriodicDerivativeOperatorQuotient,
       DerivativeOperator, DissipationOperator,
       VarCoefDerivativeOperator, SourceOfCoefficients,
       FourierDerivativeOperator, FourierConstantViscosity,
       FourierPolynomialDerivativeOperator, FourierRationalDerivativeOperator,
       FourierDerivativeOperator2D,
       LegendreDerivativeOperator, LegendreSecondDerivativeOperator
export FilterCallback, ConstantFilter, ExponentialFilter
export derivative_order, accuracy_order, source_of_coefficients, grid, semidiscretize
export mass_matrix
export integrate, left_boundary_weight, right_boundary_weight,
       derivative_left, derivative_right,
       mul_transpose_derivative_left!, mul_transpose_derivative_right!,
       evaluate_coefficients, evaluate_coefficients!,
       compute_coefficients, compute_coefficients!
export periodic_central_derivative_operator, periodic_derivative_operator, derivative_operator,
       dissipation_operator, var_coef_derivative_operator,
       fourier_derivative_operator, spectral_viscosity_operator, super_spectral_viscosity_operator,
       legendre_derivative_operator, legendre_second_derivative_operator
export UniformMesh1D, UniformPeriodicMesh1D
export couple_continuously, couple_continuosly, couple_discontinuosly, couple_discontinuously # TODO: deprecated typo
export mul!

export Fornberg1998, Holoborodko2008, BeljaddLeFlochMishraParés2017
export MattssonNordström2004, MattssonSvärdNordström2004, MattssonSvärdShoeybi2008,
       Mattsson2012, Mattsson2014,
       MattssonAlmquistCarpenter2014Extended, MattssonAlmquistCarpenter2014Optimal,
       Mattsson2017,
       MattssonAlmquistVanDerWeide2018Minimal, MattssonAlmquistVanDerWeide2018Accurate
export Tadmor1989, MadayTadmor1989, Tadmor1993,
       TadmorWaagan2012Standard, TadmorWaagan2012Convergent

export BurgersPeriodicSemidiscretization, BurgersNonperiodicSemidiscretization,
       CubicPeriodicSemidiscretization, CubicNonperiodicSemidiscretization,
       VariableLinearAdvectionNonperiodicSemidiscretization,
       WaveEquationNonperiodicSemidiscretization,
       QuarticNonconvexPeriodicSemidiscretization

end # module
