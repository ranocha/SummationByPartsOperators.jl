
"""
    Mattsson2017

Coefficients of the upwind SBP operators given in
  Mattsson (2017)
  Diagonal-norm upwind SBP operators.
  Journal of Computational Physics 335, pp.283-310.
"""
struct Mattsson2017 <: SourceOfCoefficients
  kind::Symbol

  function Mattsson2017(kind::Symbol)
    if (kind !== :plus) && (kind !== :minus) && (kind !== :central)
      throw(ArgumentError("The only choices are :plus, :minus, and :central, not :$kind."))
    end
    new(kind)
  end
end

function Base.show(io::IO, source::Mattsson2017)
print(io,
    "  Upwind coefficients (", source.kind, ") of \n",
    "  Mattsson (2017) \n",
    "  Diagonal-norm upwind SBP operators. \n",
    "  Journal of Computational Physics 335, pp.283-310. \n")
end


function first_derivative_coefficients(source::Mattsson2017, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary_plus = (
            DerivativeCoefficientRow{T,1,3}(SVector(T(-3),
                                                    T(5),
                                                    T(-2), )),
            DerivativeCoefficientRow{T,1,4}(SVector(T(-1//5),
                                                    T(-1),
                                                    T(8//5),
                                                    T(-2//5), )),
        )
        right_boundary_plus = (
            DerivativeCoefficientRow{T,1,2}(SVector(T(1),
                                                    T(-1), )),
            DerivativeCoefficientRow{T,1,2}(SVector(T(1),
                                                    T(-1), )),
        )
        upper_coef_plus = SVector( T(2),
                                   T(-1//2), )
        central_coef_plus = T(-3//2)
        lower_coef_plus = SVector{0,T}()
        left_weights = SVector( T(1//4),
                                T(5//4), )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        left_boundary_minus = (
            DerivativeCoefficientRow{T,1,2}(SVector(T(-1),
                                                    T(1), )),
            DerivativeCoefficientRow{T,1,2}(SVector(T(-1),
                                                    T(1), )),
        )
        right_boundary_minus = (
            DerivativeCoefficientRow{T,1,3}(SVector(T(3),
                                                    T(-5),
                                                    T(2), )),
            DerivativeCoefficientRow{T,1,4}(SVector(T(1//5),
                                                    T(1),
                                                    T(-8//5),
                                                    T(2//5), )),
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
                                left_weights, right_weights, parallel, 1, order, source)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end


# function second_derivative_coefficients(source::Mattsson2017, order::Int, T=Float64, parallel=Val{:serial}())
#     if order == 2
#         left_boundary = (
#             # d1
#             DerivativeCoefficientRow{T,1,3}(SVector(T(1),
#                                                     T(-2),
#                                                     T(1) )),
#         )
#         right_boundary = left_boundary
#         upper_coef = SVector(T(1))
#         central_coef = T(-2)
#         lower_coef = upper_coef
#         left_weights = SVector(T(1//2))
#         right_weights = left_weights
#         left_boundary_derivatives = (
#             DerivativeCoefficientRow{T,1,3}(SVector(T(-3//2),
#                                                     T(2),
#                                                     T(-1//2) )),
#         )
#         right_boundary_derivatives = (-left_boundary_derivatives[1], )

#         DerivativeCoefficients(left_boundary, right_boundary,
#                                 left_boundary_derivatives, right_boundary_derivatives,
#                                 lower_coef, central_coef, upper_coef,
#                                 left_weights, right_weights, parallel, 2, order, source)
#     else
#         throw(ArgumentError("Order $order not implemented/derived."))
#     end
# end
