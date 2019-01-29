module SummationByPartsOperators

using LinearAlgebra
using SparseArrays

using FFTW
using ArgCheck
using Reexport
using Requires
using Unrolled
using Parameters
using StaticArrays

@reexport using DiffEqBase
using DiffEqCallbacks

import Base: *, -
import LinearAlgebra: mul!
export mul!
@reexport using PolynomialBases
import PolynomialBases: integrate, evaluate_coefficients, evaluate_coefficients!,
                        compute_coefficients, compute_coefficients!


# types
abstract type AbstractDerivativeOperator{T} end
abstract type AbstractPeriodicDerivativeOperator{T} <: AbstractDerivativeOperator{T} end
abstract type AbstractDerivativeCoefficients{T} end
abstract type AbstractMassMatrix{T} end
abstract type AbstractSemidiscretisation end #TODO: HyperbolicDiffEq.jl; also semidiscretise
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

function __init__()
  @require BandedMatrices="aae01518-5342-5314-be14-df237901396f" include("banded_matrices.jl")
end

include("filter.jl")
include("fourier_operators.jl")
include("legendre_operators.jl")
include("SBP_coefficients/MattssonNordström2004.jl")
include("SBP_coefficients/MattssonSvärdNordström2004.jl")
include("SBP_coefficients/MattssonSvärdShoeybi2008.jl")
include("SBP_coefficients/Mattsson2012.jl")
include("SBP_coefficients/Mattsson2014.jl")
include("SBP_coefficients/MattssonAlmquistCarpenter2014Extended.jl")
include("SBP_coefficients/MattssonAlmquistCarpenter2014Optimal.jl")

include("conservation_laws/general_laws.jl")
include("conservation_laws/burgers.jl")
include("conservation_laws/cubic.jl")
include("conservation_laws/variable_linear_advection.jl")
include("second_order_eqs/wave_eq.jl")


# exports
export PeriodicDerivativeOperator, PeriodicDissipationOperator,
       DerivativeOperator, DissipationOperator,
       VarCoefDerivativeOperator, SourceOfCoefficients,
       FourierDerivativeOperator, FourierConstantViscosity,
       LegendreDerivativeOperator
export FilterCallback, ConstantFilter,
       ExponentialFilter
export derivative_order, accuracy_order, source_of_coeffcients, grid, semidiscretise
export mass_matrix
export integrate, derivative_left, derivative_right,
       evaluate_coefficients, evaluate_coefficients!,
       compute_coefficients, compute_coefficients!
export periodic_central_derivative_operator, periodic_derivative_operator, derivative_operator,
       dissipation_operator, var_coef_derivative_operator,
       fourier_derivative_operator, spectral_viscosity_operator, super_spectral_viscosity_operator,
       legendre_derivative_operator

export Fornberg1998, Holoborodko2008, BeljaddLeFlochMishraParés2017
export MattssonNordström2004, MattssonSvärdNordström2004, MattssonSvärdShoeybi2008,
       Mattsson2012, Mattsson2014,
       MattssonAlmquistCarpenter2014Extended, MattssonAlmquistCarpenter2014Optimal
export Tadmor1989, MadayTadmor1989, Tadmor1993,
       TadmorWaagan2012Standard, TadmorWaagan2012Convergent

export BurgersPeriodicSemidiscretisation, BurgersNonperiodicSemidiscretisation,
       CubicPeriodicSemidiscretisation, CubicNonperiodicSemidiscretisation,
       VariableLinearAdvectionNonperiodicSemidiscretisation,
       WaveEquationNonperiodicSemidiscretisation

end # module
