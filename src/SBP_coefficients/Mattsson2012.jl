
"""
    Mattsson2012

Coefficients of the SBP operators given in
  Mattsson (2012)
  Summation by Parts Operators for Finite Difference Approximations of
    Second-Derivatives with Variable Coefficients.
  Journal of Scientific Computing 51, pp. 650-682.
"""
struct Mattsson2012 <: SourceOfCoefficients end

function Base.show(io::IO, ::Mattsson2012)
    print(io,
        "  Mattsson (2012) \n",
        "  Summation by Parts Operators for Finite Difference Approximations of\n",
        "    Second-Derivatives with Variable Coefficients. \n",
        "  Journal of Scientific Computing 51, pp. 650-682. \n",
        "See also (first derivatives) \n",
        "  Mattsson, Nordström (2004) \n",
        "  Summation by parts operators for finite difference approximations of second \n",
        "    derivaties. \n",
        "  Journal of Computational Physics 199, pp.503-540. \n")
end


function var_coef_derivative_coefficients(source::Mattsson2012, derivative_order::Int, accuracy_order::Int, grid, parallel=Val{:serial}())
    @argcheck derivative_order == 2
    T = eltype(grid)
    if accuracy_order == 2
        coefficient_cache = Mattsson2012Cache2(T)
    elseif accuracy_order == 4
        coefficient_cache = Mattsson2012Cache4(T)
    elseif accuracy_order == 6
        coefficient_cache = Mattsson2012Cache6(T)
    elseif accuracy_order == 8
        coefficient_cache = Mattsson2012Cache8(T)
    else
        throw(ArgumentError("Order of accuracy $accuracy_order not implemented/derived."))
    end


    DissipationCoefficients(coefficient_cache, parallel, derivative_order, accuracy_order, source)
end



struct Mattsson2012Cache2{T} <: AbstractCoefficientCache{T}
    half::T

    function Mattsson2012Cache2(::Type{T}) where {T}
        half = one(T) / 2

        new{T}(half)
    end
end

lower_bandwidth(cache::Mattsson2012Cache2) = 1
upper_bandwidth(cache::Mattsson2012Cache2) = 1

function convolve_boundary_coefficients!(dest::AbstractVector, cache::Mattsson2012Cache2,
                                         u::AbstractVector, b::AbstractVector, α)
    @unpack half = cache

    @inbounds begin
        dest[  1] = α * (
                        (half*b[1] - b[2]) * u[1]
                        + (-b[1] + b[2]) * u[2]
                        + half*b[1] * u[3]
                    )

        dest[end] = α * (
                        (half*b[end] - b[end-2]) * u[end]
                        + (-b[end] + b[end-1]) * u[end-1]
                        + half*b[end] * u[end-2]
                    )
    end
end

@inline function convolve_interior_coefficients_loopbody(i, cache::Mattsson2012Cache2, u, b)
    @unpack half = cache
    @inbounds begin
        b_im1 = b[i-1]
        b_i   = b[i]
        b_ip1 = b[i+1]

        retval = (
                    -half*(b_im1 + b_i) * u[i-1]
                    + (half*b_im1 + b_i + half*b_ip1) * u[i]
                     - half*(b_i + b_ip1) * u[i+1]
                )
    end

    retval
end
