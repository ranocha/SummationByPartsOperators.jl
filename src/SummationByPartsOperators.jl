__precompile__()

module SummationByPartsOperators

using Unrolled
using ArgCheck
using Parameters
using StaticArrays

import Base: *, -
import PolynomialBases: integrate


# types
abstract type AbstractDerivativeOperator{T} end
abstract type AbstractDerivativeCoefficients{T} end
abstract type AbstractMassMatrix{T} end
"""
    SourceOfCoefficients

All sources of coefficients (articles) are subtypes of this abstract type.
"""
abstract type SourceOfCoefficients end


# source files
include("general_operators.jl")
include("periodic_operators.jl")
include("SBP_operators.jl")
include("SBP_coefficients/MattssonSvärdShoeybi2008.jl")
include("SBP_coefficients/Mattsson2014.jl")
include("SBP_coefficients/MattssonAlmquistCarpenter2014Extended.jl")


# exports
export PeriodicDerivativeOperator, DerivativeOperator, SourceOfCoefficients
export derivative_order, accuracy_order, source_of_coeffcients, grid
export mul!, integrate, derivative_left, derivative_right
export periodic_central_derivative_operator, periodic_derivative_operator, derivative_operator

export MattssonSvärdShoeybi2008, Mattsson2014, MattssonAlmquistCarpenter2014Extended

end # module
