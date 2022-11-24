
"""
    WIP(version::Symbol)

Coefficients of the first-order upwind SBP operators given in
- WIP

Higher-order coefficients are taken from
- Mattsson (2017)
  Diagonal-norm upwind SBP operators.
  Journal of Computational Physics 335, pp. 283-310.

You can choose between the different versions `:central`, `:plus`, and `:minus`.
"""
struct WIP <: SourceOfCoefficients
  kind::Symbol

  function WIP(kind::Symbol)
    if (kind !== :plus) && (kind !== :minus) && (kind !== :central)
      throw(ArgumentError("The only choices are :plus, :minus, and :central, not :$kind."))
    end
    new(kind)
  end
end

# function Base.show(io::IO, source::WIP)
#   if get(io, :compact, false)
#     summary(io, source)
#   else
#     print(io,
#         "Mattsson (2017) \n",
#         "  Diagonal-norm upwind SBP operators. \n",
#         "  Journal of Computational Physics 335, pp. 283-310. \n",
#         "  (upwind coefficients ", source.kind, ")")
#   end
# end


function first_derivative_coefficients(source::WIP, order::Int,
                                       T=Float64, mode=FastMode())
    if order == 1
        left_boundary_plus = (
            DerivativeCoefficientRow{T,1,2}(SVector(T(-1),    # -2 for continuously coupled = same
                                                    T(1), )), #  2 for continuously coupled = same
        )
        right_boundary_plus = (
            DerivativeCoefficientRow{T,1,2}(SVector(T(0),
                                                    T(0), )),
        )
        upper_coef_plus = SVector( T(1), )
        central_coef_plus = T(-1)
        lower_coef_plus = SVector{0,T}()
        left_weights = SVector( T(1), ) # 1//2 for continuously coupled = same
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        left_boundary_minus = (
            DerivativeCoefficientRow{T,1,2}(SVector(T(0),
                                                    T(0), )),
        )
        right_boundary_minus = (
            DerivativeCoefficientRow{T,1,2}(SVector(T(1),      #  2 for continuously coupled = same
                                                    T(-1), )), # -2 for continuously coupled = same
        )
        upper_coef_minus     = .- lower_coef_plus
        central_coef_minus   = .- central_coef_plus
        lower_coef_minus     = .- upper_coef_plus

        left_boundary_central  = (left_boundary_plus  .+ left_boundary_minus)  ./ 2
        right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
        upper_coef_central     = widening_plus(upper_coef_plus, upper_coef_minus) / 2
        central_coef_central   = (central_coef_plus   + central_coef_minus) / 2
        lower_coef_central     = widening_plus(lower_coef_plus, lower_coef_minus) / 2

        if source.kind === :plus
          left_boundary  = left_boundary_plus
          right_boundary = right_boundary_plus
          upper_coef     = upper_coef_plus
          central_coef   = central_coef_plus
          lower_coef     = lower_coef_plus
        elseif source.kind === :minus
          left_boundary  = left_boundary_minus
          right_boundary = right_boundary_minus
          upper_coef     = upper_coef_minus
          central_coef   = central_coef_minus
          lower_coef     = lower_coef_minus
        elseif source.kind === :central
          left_boundary  = left_boundary_central
          right_boundary = right_boundary_central
          upper_coef     = upper_coef_central
          central_coef   = central_coef_central
          lower_coef     = lower_coef_central
        end
        DerivativeCoefficients(left_boundary, right_boundary,
                               left_boundary_derivatives, right_boundary_derivatives,
                               lower_coef, central_coef, upper_coef,
                               left_weights, right_weights, mode, 1, order, source)
  else
        return first_derivative_coefficients(Mattsson2017(source.kind), order, T, mode)
    end
end
