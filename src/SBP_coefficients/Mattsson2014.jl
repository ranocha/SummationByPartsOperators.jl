
"""
    Mattsson2014

Coefficients of the SBP operators given in
  Mattsson (2014)
  Diagonal-norm summation by parts operators for fiite difference approximations
    of third and fourth derivatives.
  Journal of Computational Physics 274, pp. 432-454.
"""
struct Mattsson2014 <: SourceOfCoefficients end

function Base.show(io::IO, ::Mattsson2014)
    print(io,
        "  Mattsson (2014) \n",
        "  Diagonal-norm summation by parts operators for fiite difference approximations \n",
        "    of third and fourth derivatives. \n",
        "  Journal of Computational Physics 274, pp. 432-454. \n")
end


function first_derivative_coefficients(source::Mattsson2014, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,2}(SVector(T(-1),
                                                    T(1) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(1//2))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector(T(1//2))
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 1, order, source)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end


function second_derivative_coefficients(source::Mattsson2014, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,3}(SVector(T(1),
                                                    T(-2),
                                                    T(1) )),
        )
        right_boundary = left_boundary
        upper_coef = SVector(T(1))
        central_coef = T(-2)
        lower_coef = upper_coef
        left_weights = SVector(T(1//2))
        right_weights = left_weights
        left_boundary_derivatives = (
            # first derivative
            DerivativeCoefficientRow{T,1,3}(SVector(T(-3//2),
                                                    T(2),
                                                    T(-1//2) )),
        )
        right_boundary_derivatives = (
            -left_boundary_derivatives[1],
        )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 2, order, source)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end


function third_derivative_coefficients(source::Mattsson2014, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,4}(SVector(T(1//4),
                                                    T(-5//8),
                                                    T(1//2),
                                                    T(-1//8) )),
            # d2
            DerivativeCoefficientRow{T,1,4}(SVector(T(-11//16),
                                                    T(2),
                                                    T(-31//16),
                                                    T(5//8) )),
            # d3
            DerivativeCoefficientRow{T,1,5}(SVector(T(-1//2),
                                                    T(15//16),
                                                    T(1//8),
                                                    T(-17//16),
                                                    T(1//2) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(-1), T(1//2))
        central_coef = T(0)
        lower_coef = -upper_coef
        left_weights = SVector(T(1//2))
        right_weights = left_weights
        left_boundary_derivatives = (
            # first derivative
            DerivativeCoefficientRow{T,1,3}(SVector(T(-3//2),
                                                    T(2),
                                                    T(-1//2) )),
            # second derivative
            DerivativeCoefficientRow{T,1,3}(SVector(T(1),
                                                    T(-2),
                                                    T(1) )),
        )
        right_boundary_derivatives = (
            -left_boundary_derivatives[1],
            left_boundary_derivatives[2]
        )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 3, order, source)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end


function fourth_derivative_coefficients(source::Mattsson2014, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,4}(SVector(T(8//5),
                                                    T(-24//5),
                                                    T(24//5),
                                                    T(-8//5) )),
            # d2
            DerivativeCoefficientRow{T,1,4}(SVector(T(-2//5),
                                                    T(6//5),
                                                    T(-6//5),
                                                    T(2//5) )),
            # d3
            DerivativeCoefficientRow{T,1,5}(SVector(T(2//5),
                                                    T(-11//5),
                                                    T(21//5),
                                                    T(-17//5),
                                                    T(1) )),
            # d4
            DerivativeCoefficientRow{T,1,6}(SVector(T(1//5),
                                                    T(2//5),
                                                    T(-17//5),
                                                    T(-29//5),
                                                    T(-4),
                                                    T(1) )),
        )
        right_boundary = left_boundary
        upper_coef = SVector(T(-4), T(1))
        central_coef = T(6)
        lower_coef = upper_coef
        left_weights = SVector(T(1//2))
        right_weights = left_weights
        left_boundary_derivatives = (
            # first derivative
            DerivativeCoefficientRow{T,1,3}(SVector(T(-3//2),
                                                    T(2),
                                                    T(-1//2) )),
            # second derivative
            DerivativeCoefficientRow{T,1,3}(SVector(T(1),
                                                    T(-2),
                                                    T(1) )),
            # third derivative
            DerivativeCoefficientRow{T,1,4}(SVector(T(-1),
                                                    T(3),
                                                    T(-3),
                                                    T(1) )),
        )
        right_boundary_derivatives = (
            -left_boundary_derivatives[1],
            left_boundary_derivatives[2],
            -left_boundary_derivatives[3],
        )

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 4, order, source)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end
