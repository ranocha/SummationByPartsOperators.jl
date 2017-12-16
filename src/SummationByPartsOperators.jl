module SummationByPartsOperators

using ArgCheck
using Parameters
using StaticArrays

import PolynomialBases: integrate


# types
abstract type AbstractDerivativeOperator{T} end
abstract type AbstractMassMatrix{T} end


# source files
include("periodic_operators.jl")


# exports


end # module
