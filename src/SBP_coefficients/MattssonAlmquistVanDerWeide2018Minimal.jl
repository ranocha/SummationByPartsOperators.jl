
"""
    MattssonAlmquistVanDerWeide2018Minimal

Coefficients of the optimal SBP operators with nonuniform grid given in
  Mattsson, Almquist, van der Weide (2018)
  Boundary optimized diagonal-norm SBP operators.
  Journal of Computational Physics 374, pp. 1261-1266.
"""
struct MattssonAlmquistVanDerWeide2018Minimal <: SourceOfCoefficients end

function Base.show(io::IO, ::MattssonAlmquistVanDerWeide2018Minimal)
  print(io,
      "  Mattsson, Almquist, van der Weide (2018) \n",
      "  Boundary optimized diagonal-norm SBP operators. \n",
      "  Journal of Computational Physics 374, pp. 1261-1266. \n")
end



struct BoundaryAdaptedGrid{T,M,Grid} <: AbstractArray{T,1}
  xmin::T
  xmax::T
  xstart::SVector{M,T}
  uniform_grid::Grid
end

function BoundaryAdaptedGrid(xmin::T, xmax::T, _xstart::SVector{M,T}, N::Int) where {M,T}
  @argcheck xmin < xmax
  @argcheck N > 2 * M

  Δx = (xmax - xmin) / (2*_xstart[end] + N - 2*M - 1)
  xstart = Δx .* _xstart
  uniform_grid = range(xmin + xstart[end], xmax - xstart[end], length=N-2*M)
  BoundaryAdaptedGrid{T,M,typeof(uniform_grid)}(xmin, xmax, xstart, uniform_grid)
end

function BoundaryAdaptedGrid(xmin, xmax, xstart, N)
  T = promote_type(typeof(xmin), typeof(xmax), eltype(xstart))
  BoundaryAdaptedGrid(T(xmin), T(xmax), T.(xstart), N)
end

function Base.length(grid::BoundaryAdaptedGrid)
  length(grid.uniform_grid) + 2 * length(grid.xstart)
end

function Base.getindex(grid::BoundaryAdaptedGrid, i::Int)
  N = length(grid)
  M = length(grid.xstart)
  @boundscheck begin
    @argcheck i > 0
    @argcheck i <= N
  end

  if i <= M
    grid.xmin + grid.xstart[i]
  elseif i >= N-M+1
    grid.xmax - grid.xstart[end-(N-i)]
  else
    @inbounds x = getindex(grid.uniform_grid, i-M)
    x
  end
end

Base.size(grid::BoundaryAdaptedGrid) = (length(grid),)
Base.step(grid::BoundaryAdaptedGrid) = step(grid.uniform_grid)


function construct_grid(::MattssonAlmquistVanDerWeide2018Minimal, accuracy_order, xmin, xmax, N)
  @argcheck N > 6
  T = promote_type(typeof(xmin), typeof(xmax))

  if accuracy_order == 4
    xstart = SVector(
      T(0.0000000000000e+00),
      T(7.7122987842562e-01),
      T(1.7712298784256e+00),
      T(2.7712298784256e+00),
    )
  elseif accuracy_order == 6
    xstart = SVector(
      T(0.0000000000000e+00),
      T(4.0842950991998e-01),
      T(1.1968523189207e+00),
      T(2.1968523189207e+00),
      T(3.1968523189207e+00;),
      T(4.1968523189207e+00),
    )
  else
    throw(ArgumentError("Accuracy order $accuracy_order not implemented/derived."))
  end

  BoundaryAdaptedGrid(xmin, xmax, xstart, N)
end


function first_derivative_coefficients(source::MattssonAlmquistVanDerWeide2018Minimal, order::Int, T=Float64, parallel=Val{:serial}())
  if order == 4
    left_boundary = (
        # d1
        DerivativeCoefficientRow{T,1,3}(SVector(T(-1.8612097181862932),
                                                T(2.2966302628676907),
                                                T(-0.4354205446813973) )),
        # d2
        DerivativeCoefficientRow{T,1,4}(SVector(T(-0.611186527373936),
                                                T(0),
                                                T(0.6937383659498486),
                                                T(-0.08255183857591565) )),
        # d3
        DerivativeCoefficientRow{T,1,5}(SVector(T(0.11587549591844003),
                                                T(-0.6937383659498486),
                                                T(0),
                                                T(0.6604147086073278),
                                                T(-0.08255183857591597) )),
    )
    right_boundary = .- left_boundary
    upper_coef = SVector(T(2//3), T(-1//12))
    central_coef = zero(T)
    lower_coef = -upper_coef
    left_weights = SVector( T(0.26864248295847),
                            T(1.00946671535),
                            T(1.00946671535) )
    right_weights = left_weights
    left_boundary_derivatives = Tuple{}()
    right_boundary_derivatives = left_boundary_derivatives

    DerivativeCoefficients(left_boundary, right_boundary,
                            left_boundary_derivatives, right_boundary_derivatives,
                            lower_coef, central_coef, upper_coef,
                            left_weights, right_weights, parallel, 1, order, source)
  elseif order == 6
    left_boundary = (
        # d1
        DerivativeCoefficientRow{T,1,5}(SVector(T(-3.924566448353279),
                                                T(4.962014957077398),
                                                T(-1.2881968205658438),
                                                T(0.28645730739095715),
                                                T(-0.03570899554923223) )),
        # d2
        DerivativeCoefficientRow{T,1,5}(SVector(T(-1.022587534557978),
                                                T(0),
                                                T(1.3023996740158794),
                                                T(-0.3357380168806994),
                                                T(0.05592587742279807) )),
        # d3
        DerivativeCoefficientRow{T,1,6}(SVector(T(0.17402334848464412),
                                                T(-0.8537429915394141),
                                                T(0),
                                                T(0.8419418935035105),
                                                T(-0.179894658468737),
                                                T(0.017672408019993497) )),
        # d4
        DerivativeCoefficientRow{T,1,7}(SVector(T(-0.036159059808532096),
                                                T(0.20564365765539105),
                                                T(-0.7867088686169533),
                                                T(0),
                                                T(0.749328780370319),
                                                T(-0.1486175733002554),
                                                T(0.01651306370002838) )),
        # d5
        DerivativeCoefficientRow{T,1,8}(SVector(T(0.004554664962414826),
                                                T(-0.03461379228073425),
                                                T(0.16985243300141972),
                                                T(-0.7571710333787599),
                                                T(0),
                                                T(0.7508648039541874),
                                                T(-0.15017296079083747),
                                                T(0.016685884532315277) )),
    )
    right_boundary = .- left_boundary
    upper_coef = SVector(T(3//4), T(-3//20), T(1//60))
    central_coef = zero(T)
    lower_coef = -upper_coef
    left_weights = SVector( T(0.12740260779883),
                            T(0.61820981002054),
                            T(0.94308973897679),
                            T(1.0093019060199),
                            T(0.99884825610465) )
    right_weights = left_weights
    left_boundary_derivatives = Tuple{}()
    right_boundary_derivatives = left_boundary_derivatives

    DerivativeCoefficients(left_boundary, right_boundary,
                            left_boundary_derivatives, right_boundary_derivatives,
                            lower_coef, central_coef, upper_coef,
                            left_weights, right_weights, parallel, 1, order, source)
  # elseif order == 6
  #   left_boundary = (
  #       # d1
  #       DerivativeCoefficientRow{T,1,5}(SVector(T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T() )),
  #       # d2
  #       DerivativeCoefficientRow{T,1,5}(SVector(T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T() )),
  #       # d3
  #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T() )),
  #       # d4
  #       DerivativeCoefficientRow{T,1,7}(SVector(T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T() )),
  #       # d5
  #       DerivativeCoefficientRow{T,1,8}(SVector(T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T(),
  #                                               T() )),
  #   )
  #   right_boundary = .- left_boundary
  #   upper_coef = SVector(T(2//3), T(-1//12))
  #   central_coef = zero(T)
  #   lower_coef = -upper_coef
  #   left_weights = SVector( T(),
  #                           T(),
  #                           T() )
  #   right_weights = left_weights
  #   left_boundary_derivatives = Tuple{}()
  #   right_boundary_derivatives = left_boundary_derivatives

  #   DerivativeCoefficients(left_boundary, right_boundary,
  #                           left_boundary_derivatives, right_boundary_derivatives,
  #                           lower_coef, central_coef, upper_coef,
  #                           left_weights, right_weights, parallel, 1, order, source)
  else
    throw(ArgumentError("Order $order not implemented/derived."))
  end
end
