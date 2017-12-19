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


# source files
include("general_operators.jl")
include("periodic_operators.jl")
include("SBP_operators.jl")


# exports
export PeriodicDerivativeOperator, DerivativeOperator
export derivative_order, accuracy_order, source_of_coeffcients
export mul!
export periodic_central_derivative_operator, periodic_derivative_operator, derivative_operator

export MattssonSvärdShoeybi2008

end # module
