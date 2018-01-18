__precompile__()

module SummationByPartsOperators

using Requires
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
include("dissipation_operators.jl")
include("var_coef_operators.jl")
@require BandedMatrices begin
        include("banded_matrices.jl")
end
include("SBP_coefficients/MattssonNordström2004.jl")
include("SBP_coefficients/MattssonSvärdNordström2004.jl")
include("SBP_coefficients/MattssonSvärdShoeybi2008.jl")
include("SBP_coefficients/Mattsson2012.jl")
include("SBP_coefficients/Mattsson2014.jl")
include("SBP_coefficients/MattssonAlmquistCarpenter2014Extended.jl")
include("SBP_coefficients/MattssonAlmquistCarpenter2014Optimal.jl")


# exports
export PeriodicDerivativeOperator, DerivativeOperator, DissipationOperator,
       VarCoefDerivativeOperator, SourceOfCoefficients
export derivative_order, accuracy_order, source_of_coeffcients, grid
export mass_matrix
export mul!, integrate, derivative_left, derivative_right
export periodic_central_derivative_operator, periodic_derivative_operator, derivative_operator,
        dissipation_operator, var_coef_derivative_operator

export MattssonNordström2004, MattssonSvärdNordström2004, MattssonSvärdShoeybi2008,
        Mattsson2012, Mattsson2014,
        MattssonAlmquistCarpenter2014Extended, MattssonAlmquistCarpenter2014Optimal

end # module
