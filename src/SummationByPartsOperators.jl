
"""
    SummationByPartsOperators

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
for further information. Some noticeable functions to start with are
[`derivative_operator`](@ref),
[`legendre_derivative_operator`](@ref),
[`periodic_derivative_operator`](@ref),
[`fourier_derivative_operator`](@ref),
[`dissipation_operator`](@ref),
and [`grid`](@ref).

If you use this package for your research, please cite it using
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
"""
module SummationByPartsOperators

using LinearAlgebra
using SparseArrays

using AutoHashEquals: @auto_hash_equals
using ArgCheck: @argcheck
using StaticArrayInterface: StaticArrayInterface, StaticInt, static_length
using FFTW: FFTW, plan_rfft, plan_brfft, plan_irfft
using LoopVectorization: LoopVectorization, @turbo, @tturbo
using MuladdMacro: @muladd
using RecursiveArrayTools: recursive_bottom_eltype
using Reexport: @reexport
if !isdefined(Base, :get_extension)
    using Requires: @require
end
using StaticArrays: SVector, StaticVector, FieldVector
using SimpleUnPack: @unpack
using Unrolled: @unroll

@reexport using SciMLBase: SciMLBase, DiscreteCallback, ODEProblem,
                           SecondOrderODEProblem

import LinearAlgebra: mul!
@reexport using PolynomialBases
import PolynomialBases: grid, mass_matrix, mass_matrix_boundary, integrate,
                        evaluate_coefficients, evaluate_coefficients!,
                        compute_coefficients, compute_coefficients!

# types
abstract type AbstractDerivativeOperator{T} end
abstract type AbstractNonperiodicDerivativeOperator{T} <: AbstractDerivativeOperator{T} end
abstract type AbstractPeriodicDerivativeOperator{T} <: AbstractDerivativeOperator{T} end
abstract type AbstractMatrixDerivativeOperator{T} <:
              AbstractNonperiodicDerivativeOperator{T} end
abstract type AbstractMultidimensionalMatrixDerivativeOperator{Dim, T} <:
              AbstractMatrixDerivativeOperator{T} end
abstract type AbstractDerivativeCoefficients{T} end
abstract type AbstractMassMatrix{T} end
abstract type AbstractSemidiscretization end
abstract type AbstractFilter{T <: Real} end
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
    @static if !isdefined(Base, :get_extension)
        @require BandedMatrices="aae01518-5342-5314-be14-df237901396f" include("../ext/SummationByPartsOperatorsBandedMatricesExt.jl")
        @require DiffEqCallbacks="459566f4-90b8-5000-8ac3-15dfb0a30def" include("../ext/SummationByPartsOperatorsDiffEqCallbacksExt.jl")
        @require ForwardDiff="f6369f11-7733-5829-9624-2563aa707210" include("../ext/SummationByPartsOperatorsForwardDiffExt.jl")
        @require StructArrays="09ab397b-f2b6-538f-b94a-2f83cf4a842a" include("../ext/SummationByPartsOperatorsStructArraysExt.jl")
    end
end

include("filter.jl")
include("fourier_operators.jl")
include("fourier_operators_2d.jl")
include("legendre_operators.jl")
include("matrix_operators.jl")
include("multidimensional_matrix_operators.jl")
include("tensor_product_operators.jl")
include("function_space_operators.jl")
include("upwind_operators.jl")
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
include("SBP_coefficients/MattssonNiemeläWinters2026.jl")
include("SBP_coefficients/DienerDorbandSchnetterTiglio2007.jl")
include("SBP_coefficients/SharanBradyLivescu2022.jl")
include("SBP_coefficients/WilliamsDuru2024.jl")

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
       VarCoefDerivativeOperator, SourceOfCoefficients, SourceOfCoefficientsCombination,
       FourierDerivativeOperator, FourierConstantViscosity,
       FourierPolynomialDerivativeOperator, FourierRationalDerivativeOperator,
       FourierDerivativeOperator2D,
       LegendreDerivativeOperator, LegendreSecondDerivativeOperator,
       MatrixDerivativeOperator, MultidimensionalMatrixDerivativeOperator,
       TensorProductOperator,
       UpwindOperators, PeriodicUpwindOperators
export FilterCallback, ConstantFilter, ExponentialFilter
export SafeMode, FastMode, ThreadedMode
export derivative_order, accuracy_order, source_of_coefficients, grid, semidiscretize
export mass_matrix, mass_matrix_boundary
export integrate, integrate_boundary,
       restrict_interior, restrict_boundary,
       left_boundary_weight, right_boundary_weight,
       scale_by_mass_matrix!, scale_by_inverse_mass_matrix!,
       derivative_left, derivative_right,
       mul_transpose_derivative_left!, mul_transpose_derivative_right!,
       evaluate_coefficients, evaluate_coefficients!,
       compute_coefficients, compute_coefficients!
export periodic_central_derivative_operator, periodic_derivative_operator,
       derivative_operator,
       dissipation_operator, var_coef_derivative_operator,
       fourier_derivative_operator,
       legendre_derivative_operator, legendre_second_derivative_operator,
       upwind_operators, function_space_operator, tensor_product_operator_2D
export normals, boundary_indices
export UniformMesh1D, UniformPeriodicMesh1D
export couple_continuously, couple_discontinuously
export mul!

export Fornberg1998, Holoborodko2008, LanczosLowNoise, BeljaddLeFlochMishraParés2017
export MattssonNordström2004, MattssonSvärdNordström2004, MattssonSvärdShoeybi2008,
       Mattsson2012, Mattsson2014,
       MattssonAlmquistCarpenter2014Extended, MattssonAlmquistCarpenter2014Optimal,
       Mattsson2017,
       MattssonAlmquistVanDerWeide2018Minimal, MattssonAlmquistVanDerWeide2018Accurate,
       MattssonNiemeläWinters2026,
       DienerDorbandSchnetterTiglio2007,
       SharanBradyLivescu2022,
       WilliamsDuru2024
export Tadmor1989, MadayTadmor1989, Tadmor1993,
       TadmorWaagan2012Standard, TadmorWaagan2012Convergent
export GlaubitzNordströmÖffner2023

export BurgersPeriodicSemidiscretization, BurgersNonperiodicSemidiscretization,
       CubicPeriodicSemidiscretization, CubicNonperiodicSemidiscretization,
       VariableLinearAdvectionNonperiodicSemidiscretization,
       VariableLinearAdvectionPeriodicSemidiscretization,
       WaveEquationNonperiodicSemidiscretization,
       QuarticNonconvexPeriodicSemidiscretization

# explicit precompilation only on Julia v1.9 and newer
@static if VERSION >= v"1.9.0-beta4"
    include("precompile.jl")
end

end # module
