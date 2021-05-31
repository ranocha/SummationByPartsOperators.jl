
"""
    Mattsson2017(version::Symbol)

Coefficients of the upwind SBP operators given in
- Mattsson (2017)
  Diagonal-norm upwind SBP operators.
  Journal of Computational Physics 335, pp. 283-310.

You can choose between the different versions `:central`, `:plus`, and `:minus`.
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
  if get(io, :compact, false)
    summary(io, source)
  else
    print(io,
        "Mattsson (2017) \n",
        "  Diagonal-norm upwind SBP operators. \n",
        "  Journal of Computational Physics 335, pp. 283-310. \n",
        "  (upwind coefficients ", source.kind, ")")
  end
end


function first_derivative_coefficients(source::Mattsson2017, order::Int,
                                       T=Float64, mode=FastMode())
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
                                left_weights, right_weights, mode, 1, order, source)
    elseif order == 3
      left_boundary_plus = (
          DerivativeCoefficientRow{T,1,3}(SVector(T(-7//5),
                                                  T(9//5),
                                                  T(-2//5), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(-5//13),
                                                  T(-5//13),
                                                  T(12//13),
                                                  T(-2//13), )),
      )
      right_boundary_plus = (
          DerivativeCoefficientRow{T,1,2}(SVector(T(1),
                                                  T(-1), )),
          DerivativeCoefficientRow{T,1,3}(SVector(T(9//13),
                                                  T(-5//13),
                                                  T(-4//13), )),
      )
      upper_coef_plus = SVector( T(1),
                                 T(-1//6), )
      central_coef_plus = T(-1//2)
      lower_coef_plus = SVector( T(-1//3), )
      left_weights = SVector( T(5//12),
                              T(13//12), )
      right_weights = left_weights
      left_boundary_derivatives = Tuple{}()
      right_boundary_derivatives = left_boundary_derivatives

      left_boundary_minus = (
          DerivativeCoefficientRow{T,1,2}(SVector(T(-1),
                                                  T(1), )),
          DerivativeCoefficientRow{T,1,3}(SVector(T(-9//13),
                                                  T(5//13),
                                                  T(4//13), )),
      )
      right_boundary_minus = (
          DerivativeCoefficientRow{T,1,3}(SVector(T(7//5),
                                                  T(-9//5),
                                                  T(2//5), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(5//13),
                                                  T(5//13),
                                                  T(-12//13),
                                                  T(2//13), )),
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
    elseif order == 4
      left_boundary_plus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(-75//49),
                                                  T(205//98),
                                                  T(-29//49),
                                                  T(3//98), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(-169//366),
                                                  T(-11//61),
                                                  T(99//122),
                                                  T(-43//183),
                                                  T(4//61), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(11//123),
                                                  T(-39//82),
                                                  T(-29//41),
                                                  T(389//246),
                                                  T(-24//41),
                                                  T(4//41), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(9//298),
                                                  T(-11//149),
                                                  T(-65//298),
                                                  T(-117//149),
                                                  T(216//149),
                                                  T(-72//149),
                                                  T(12//149), )),
      )
      right_boundary_plus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(69//49),
                                                  T(-169//98),
                                                  T(11//49),
                                                  T(9//98), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(205//366),
                                                  T(-11//61),
                                                  T(-39//122),
                                                  T(-11//183), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(-29//123),
                                                  T(99//82),
                                                  T(-29//41),
                                                  T(-65//246), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(3//298),
                                                  T(-43//149),
                                                  T(389//298),
                                                  T(-117//149),
                                                  T(-36//149), )),
      )
      upper_coef_plus = SVector( T(3//2),
                                 T(-1//2),
                                 T(1//12), )
      central_coef_plus = T(-5//6)
      lower_coef_plus = SVector( T(-1//4), )
      left_weights = SVector( T(49//144),
                              T(61//48),
                              T(41//48),
                              T(149//144), )
      right_weights = left_weights
      left_boundary_derivatives = Tuple{}()
      right_boundary_derivatives = left_boundary_derivatives

      left_boundary_minus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(-69//49),
                                                  T(169//98),
                                                  T(-11//49),
                                                  T(-9//98), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(-205//366),
                                                  T(11//61),
                                                  T(39//122),
                                                  T(11//183), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(29//123),
                                                  T(-99//82),
                                                  T(29//41),
                                                  T(65//246), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(-3//298),
                                                  T(43//149),
                                                  T(-389//298),
                                                  T(117//149),
                                                  T(36//149), )),
      )
      right_boundary_minus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(75//49),
                                                  T(-205//98),
                                                  T(29//49),
                                                  T(-3//98), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(169//366),
                                                  T(11//61),
                                                  T(-99//122),
                                                  T(43//183),
                                                  T(-4//61), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-11//123),
                                                  T(39//82),
                                                  T(29//41),
                                                  T(-389//246),
                                                  T(24//41),
                                                  T(-4//41), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(-9//298),
                                                  T(11//149),
                                                  T(65//298),
                                                  T(117//149),
                                                  T(-216//149),
                                                  T(72//149),
                                                  T(-12//149), )),
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
    elseif order == 5
      left_boundary_plus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(-366//251),
                                                  T(941//502),
                                                  T(-94//251),
                                                  T(-21//502), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(-869//1794),
                                                  T(-22//299),
                                                  T(375//598),
                                                  T(-86//897),
                                                  T(8//299), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(58//633),
                                                  T(-255//422),
                                                  T(-58//211),
                                                  T(1309//1266),
                                                  T(-60//211),
                                                  T(8//211), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(45//1478),
                                                  T(-22//739),
                                                  T(-661//1478),
                                                  T(-234//739),
                                                  T(720//739),
                                                  T(-180//739),
                                                  T(24//739), )),
      )
      right_boundary_plus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(354//251),
                                                  T(-869//502),
                                                  T(58//251),
                                                  T(45//502), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(941//1794),
                                                  T(-22//299),
                                                  T(-255//598),
                                                  T(-22//897), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(-94//633),
                                                  T(375//422),
                                                  T(-58//211),
                                                  T(-661//1266),
                                                  T(12//211), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-21//1478),
                                                  T(-86//739),
                                                  T(1309//1478),
                                                  T(-234//739),
                                                  T(-360//739),
                                                  T(36//739), )),
      )
      upper_coef_plus = SVector( T(1),
                                 T(-1//4),
                                 T(1//30), )
      central_coef_plus = T(-1//3)
      lower_coef_plus = SVector( T(-1//2),
                                 T(1//20), )
      left_weights = SVector( T(251//720),
                              T(299//240),
                              T(211//240),
                              T(739//720), )
      right_weights = left_weights
      left_boundary_derivatives = Tuple{}()
      right_boundary_derivatives = left_boundary_derivatives

      left_boundary_minus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(-354//251),
                                                  T(869//502),
                                                  T(-58//251),
                                                  T(-45//502), )),
          DerivativeCoefficientRow{T,1,4}(SVector(T(-941//1794),
                                                  T(22//299),
                                                  T(255//598),
                                                  T(22//897), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(94//633),
                                                  T(-375//422),
                                                  T(58//211),
                                                  T(661//1266),
                                                  T(-12//211), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(21//1478),
                                                  T(86//739),
                                                  T(-1309//1478),
                                                  T(234//739),
                                                  T(360//739),
                                                  T(-36//739), )),
      )
      right_boundary_minus = (
          DerivativeCoefficientRow{T,1,4}(SVector(T(366//251),
                                                  T(-941//502),
                                                  T(94//251),
                                                  T(21//502), )),
          DerivativeCoefficientRow{T,1,5}(SVector(T(869//1794),
                                                  T(22//299),
                                                  T(-375//598),
                                                  T(86//897),
                                                  T(-8//299), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-58//633),
                                                  T(255//422),
                                                  T(58//211),
                                                  T(-1309//1266),
                                                  T(60//211),
                                                  T(-8//211), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(-45//1478),
                                                  T(22//739),
                                                  T(661//1478),
                                                  T(234//739),
                                                  T(-720//739),
                                                  T(180//739),
                                                  T(-24//739), )),
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
    elseif order == 6
      left_boundary_plus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(-58148100//36496453),
                                                  T(1146190567//547446795),
                                                  T(-14369571//52137790),
                                                  T(-55265831//182482265),
                                                  T(26269819//1094893590),
                                                  T(9858004//182482265), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-1116490567//2422752675),
                                                  T(-954612//32303369),
                                                  T(190538869//484550535),
                                                  T(102705469//969101070),
                                                  T(4964892//161516845),
                                                  T(-191689861//4845505350), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(9869571//102452500),
                                                  T(-135385429//215150250),
                                                  T(-2198412//7171675),
                                                  T(45137333//35858375),
                                                  T(-253641811//430300500),
                                                  T(70665929//358583750),
                                                  T(-72//2675), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(66965831//723199750),
                                                  T(-208765789//867839700),
                                                  T(-17623253//72319975),
                                                  T(-657684//2066285),
                                                  T(410905829//433919850),
                                                  T(-477953317//1446399500),
                                                  T(576//5395),
                                                  T(-72//5395), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(-49219819//3153258150),
                                                  T(3519588//105108605),
                                                  T(26422771//630651630),
                                                  T(-141938309//315325815),
                                                  T(-12476988//21021721),
                                                  T(2217185207//1576629075),
                                                  T(-4320//7841),
                                                  T(1152//7841),
                                                  T(-144//7841), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(-9498004//587634985),
                                                  T(142906261//3525809910),
                                                  T(-3137129//587634985),
                                                  T(-29884283//1175269970),
                                                  T(-630168407//1762904955),
                                                  T(-9609300//16789571),
                                                  T(57600//43837),
                                                  T(-21600//43837),
                                                  T(5760//43837),
                                                  T(-720//43837), )),
      )
      right_boundary_plus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(57671100//36496453),
                                                  T(-1116490567//547446795),
                                                  T(9869571//52137790),
                                                  T(66965831//182482265),
                                                  T(-49219819//1094893590),
                                                  T(-9498004//182482265), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(1146190567//2422752675),
                                                  T(-954612//32303369),
                                                  T(-135385429//484550535),
                                                  T(-208765789//969101070),
                                                  T(3519588//161516845),
                                                  T(142906261//4845505350), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-14369571//102452500),
                                                  T(190538869//215150250),
                                                  T(-2198412//7171675),
                                                  T(-17623253//35858375),
                                                  T(26422771//430300500),
                                                  T(-3137129//358583750), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-55265831//723199750),
                                                  T(102705469//867839700),
                                                  T(45137333//72319975),
                                                  T(-657684//2066285),
                                                  T(-141938309//433919850),
                                                  T(-2298791//111261500), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(26269819//3153258150),
                                                  T(4964892//105108605),
                                                  T(-253641811//630651630),
                                                  T(410905829//315325815),
                                                  T(-12476988//21021721),
                                                  T(-630168407//1576629075),
                                                  T(288//7841), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(9858004//587634985),
                                                  T(-191689861//3525809910),
                                                  T(70665929//587634985),
                                                  T(-477953317//1175269970),
                                                  T(2217185207//1762904955),
                                                  T(-9609300//16789571),
                                                  T(-17280//43837),
                                                  T(1440//43837), )),
      )
      upper_coef_plus = SVector( T(4//3),
                                 T(-1//2),
                                 T(2//15),
                                 T(-1//60), )
      central_coef_plus = T(-7//12)
      lower_coef_plus = SVector( T(-2//5),
                                 T(1//30), )
      left_weights = SVector( T(13613//43200),
                              T(12049//8640),
                              T(535//864),
                              T(1079//864),
                              T(7841//8640),
                              T(43837//43200), )
      right_weights = left_weights
      left_boundary_derivatives = Tuple{}()
      right_boundary_derivatives = left_boundary_derivatives

      left_boundary_minus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(-57671100//36496453),
                                                  T(1116490567//547446795),
                                                  T(-9869571//52137790),
                                                  T(-66965831//182482265),
                                                  T(49219819//1094893590),
                                                  T(9498004//182482265), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-1146190567//2422752675),
                                                  T(954612//32303369),
                                                  T(135385429//484550535),
                                                  T(208765789//969101070),
                                                  T(-3519588//161516845),
                                                  T(-142906261//4845505350), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(14369571//102452500),
                                                  T(-190538869//215150250),
                                                  T(2198412//7171675),
                                                  T(17623253//35858375),
                                                  T(-26422771//430300500),
                                                  T(3137129//358583750), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(55265831//723199750),
                                                  T(-102705469//867839700),
                                                  T(-45137333//72319975),
                                                  T(657684//2066285),
                                                  T(141938309//433919850),
                                                  T(2298791//111261500), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(-26269819//3153258150),
                                                  T(-4964892//105108605),
                                                  T(253641811//630651630),
                                                  T(-410905829//315325815),
                                                  T(12476988//21021721),
                                                  T(630168407//1576629075),
                                                  T(-288//7841), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-9858004//587634985),
                                                  T(191689861//3525809910),
                                                  T(-70665929//587634985),
                                                  T(477953317//1175269970),
                                                  T(-2217185207//1762904955),
                                                  T(9609300//16789571),
                                                  T(17280//43837),
                                                  T(-1440//43837), )),
      )
      right_boundary_minus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(58148100//36496453),
                                                  T(-1146190567//547446795),
                                                  T(14369571//52137790),
                                                  T(55265831//182482265),
                                                  T(-26269819//1094893590),
                                                  T(-9858004//182482265), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(1116490567//2422752675),
                                                  T(954612//32303369),
                                                  T(-190538869//484550535),
                                                  T(-102705469//969101070),
                                                  T(-4964892//161516845),
                                                  T(191689861//4845505350), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(-9869571//102452500),
                                                  T(135385429//215150250),
                                                  T(2198412//7171675),
                                                  T(-45137333//35858375),
                                                  T(253641811//430300500),
                                                  T(-70665929//358583750),
                                                  T(72//2675), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-66965831//723199750),
                                                  T(208765789//867839700),
                                                  T(17623253//72319975),
                                                  T(657684//2066285),
                                                  T(-410905829//433919850),
                                                  T(477953317//1446399500),
                                                  T(-576//5395),
                                                  T(72//5395), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(49219819//3153258150),
                                                  T(-3519588//105108605),
                                                  T(-26422771//630651630),
                                                  T(141938309//315325815),
                                                  T(12476988//21021721),
                                                  T(-2217185207//1576629075),
                                                  T(4320//7841),
                                                  T(-1152//7841),
                                                  T(144//7841), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(9498004//587634985),
                                                  T(-142906261//3525809910),
                                                  T(3137129//587634985),
                                                  T(29884283//1175269970),
                                                  T(630168407//1762904955),
                                                  T(9609300//16789571),
                                                  T(-57600//43837),
                                                  T(21600//43837),
                                                  T(-5760//43837),
                                                  T(720//43837), )),
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
    elseif order == 7
      left_boundary_plus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(-81216540//51172247),
                                                  T(1587945773//767583705),
                                                  T(-17337249//73103210),
                                                  T(-84398989//255861235),
                                                  T(48781961//1535167410),
                                                  T(13716476//255861235), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-1570125773//3386062785),
                                                  T(-2863836//225737519),
                                                  T(240029831//677212557),
                                                  T(202934303//1354425114),
                                                  T(1418484//225737519),
                                                  T(-231357719//6772125570), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(14637249//144536540),
                                                  T(-206937767//303526734),
                                                  T(-6595236//50587789),
                                                  T(49602727//50587789),
                                                  T(-218919665//607053468),
                                                  T(51815011//505877890),
                                                  T(-216//18869), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(91418989//1008619010),
                                                  T(-266570495//1210342812),
                                                  T(-33094279//100861901),
                                                  T(-1973052//14408843),
                                                  T(440626231//605171406),
                                                  T(-365711063//2017238020),
                                                  T(2016//37621),
                                                  T(-216//37621), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(-62551961//4426143330),
                                                  T(9588//385217),
                                                  T(82588241//885228666),
                                                  T(-279245719//442614333),
                                                  T(-37430964//147538111),
                                                  T(2312302333//2213071665),
                                                  T(-18144//55031),
                                                  T(4032//55031),
                                                  T(-432//55031), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(-13500476//822302915),
                                                  T(202087559//4933817490),
                                                  T(-11297731//822302915),
                                                  T(61008503//1644605830),
                                                  T(-1360092253//2466908745),
                                                  T(-5765580//23494369),
                                                  T(60480//61343),
                                                  T(-18144//61343),
                                                  T(4032//61343),
                                                  T(-432//61343), )),
      )
      right_boundary_plus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(80930340//51172247),
                                                  T(-1570125773//767583705),
                                                  T(14637249//73103210),
                                                  T(91418989//255861235),
                                                  T(-62551961//1535167410),
                                                  T(-13500476//255861235), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(1587945773//3386062785),
                                                  T(-2863836//225737519),
                                                  T(-206937767//677212557),
                                                  T(-266570495//1354425114),
                                                  T(9588//589393),
                                                  T(202087559//6772125570), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-17337249//144536540),
                                                  T(240029831//303526734),
                                                  T(-6595236//50587789),
                                                  T(-33094279//50587789),
                                                  T(82588241//607053468),
                                                  T(-11297731//505877890), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(-84398989//1008619010),
                                                  T(202934303//1210342812),
                                                  T(49602727//100861901),
                                                  T(-1973052//14408843),
                                                  T(-279245719//605171406),
                                                  T(61008503//2017238020),
                                                  T(-288//37621), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(48781961//4426143330),
                                                  T(1418484//147538111),
                                                  T(-218919665//885228666),
                                                  T(440626231//442614333),
                                                  T(-37430964//147538111),
                                                  T(-1360092253//2213071665),
                                                  T(6048//55031),
                                                  T(-576//55031), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(13716476//822302915),
                                                  T(-231357719//4933817490),
                                                  T(51815011//822302915),
                                                  T(-365711063//1644605830),
                                                  T(2312302333//2466908745),
                                                  T(-5765580//23494369),
                                                  T(-36288//61343),
                                                  T(6048//61343),
                                                  T(-576//61343), )),
      )
      upper_coef_plus = SVector( T(1),
                                 T(-3//10),
                                 T(1//15),
                                 T(-1//140), )
      central_coef_plus = T(-1//4)
      lower_coef_plus = SVector( T(-3//5),
                                 T(1//10),
                                 T(-1//105), )
      left_weights = SVector( T(19087//60480),
                              T(84199//60480),
                              T(18869//30240),
                              T(37621//30240),
                              T(55031//60480),
                              T(61343//60480), )
      right_weights = left_weights
      left_boundary_derivatives = Tuple{}()
      right_boundary_derivatives = left_boundary_derivatives

      left_boundary_minus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(-80930340//51172247),
                                                  T(1570125773//767583705),
                                                  T(-14637249//73103210),
                                                  T(-91418989//255861235),
                                                  T(62551961//1535167410),
                                                  T(13500476//255861235), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(-1587945773//3386062785),
                                                  T(2863836//225737519),
                                                  T(206937767//677212557),
                                                  T(266570495//1354425114),
                                                  T(-9588//589393),
                                                  T(-202087559//6772125570), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(17337249//144536540),
                                                  T(-240029831//303526734),
                                                  T(6595236//50587789),
                                                  T(33094279//50587789),
                                                  T(-82588241//607053468),
                                                  T(11297731//505877890), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(84398989//1008619010),
                                                  T(-202934303//1210342812),
                                                  T(-49602727//100861901),
                                                  T(1973052//14408843),
                                                  T(279245719//605171406),
                                                  T(-61008503//2017238020),
                                                  T(288//37621), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-48781961//4426143330),
                                                  T(-1418484//147538111),
                                                  T(218919665//885228666),
                                                  T(-440626231//442614333),
                                                  T(37430964//147538111),
                                                  T(1360092253//2213071665),
                                                  T(-6048//55031),
                                                  T(576//55031), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(-13716476//822302915),
                                                  T(231357719//4933817490),
                                                  T(-51815011//822302915),
                                                  T(365711063//1644605830),
                                                  T(-2312302333//2466908745),
                                                  T(5765580//23494369),
                                                  T(36288//61343),
                                                  T(-6048//61343),
                                                  T(576//61343), )),
      )
      right_boundary_minus = (
          DerivativeCoefficientRow{T,1,6}(SVector(T(81216540//51172247),
                                                  T(-1587945773//767583705),
                                                  T(17337249//73103210),
                                                  T(84398989//255861235),
                                                  T(-48781961//1535167410),
                                                  T(-13716476//255861235), )),
          DerivativeCoefficientRow{T,1,6}(SVector(T(1570125773//3386062785),
                                                  T(2863836//225737519),
                                                  T(-240029831//677212557),
                                                  T(-202934303//1354425114),
                                                  T(-1418484//225737519),
                                                  T(231357719//6772125570), )),
          DerivativeCoefficientRow{T,1,7}(SVector(T(-14637249//144536540),
                                                  T(206937767//303526734),
                                                  T(6595236//50587789),
                                                  T(-49602727//50587789),
                                                  T(218919665//607053468),
                                                  T(-51815011//505877890),
                                                  T(216//18869), )),
          DerivativeCoefficientRow{T,1,8}(SVector(T(-91418989//1008619010),
                                                  T(266570495//1210342812),
                                                  T(33094279//100861901),
                                                  T(1973052//14408843),
                                                  T(-440626231//605171406),
                                                  T(365711063//2017238020),
                                                  T(-2016//37621),
                                                  T(216//37621), )),
          DerivativeCoefficientRow{T,1,9}(SVector(T(62551961//4426143330),
                                                  T(-9588//385217),
                                                  T(-82588241//885228666),
                                                  T(279245719//442614333),
                                                  T(37430964//147538111),
                                                  T(-2312302333//2213071665),
                                                  T(18144//55031),
                                                  T(-4032//55031),
                                                  T(432//55031), )),
          DerivativeCoefficientRow{T,1,10}(SVector(T(13500476//822302915),
                                                  T(-202087559//4933817490),
                                                  T(11297731//822302915),
                                                  T(-61008503//1644605830),
                                                  T(1360092253//2466908745),
                                                  T(5765580//23494369),
                                                  T(-60480//61343),
                                                  T(18144//61343),
                                                  T(-4032//61343),
                                                  T(432//61343), )),
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
    # elseif order == 6
    #   left_boundary_plus = (
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,7}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,8}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,9}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,10}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #   )
    #   right_boundary_plus = (
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,7}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,8}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #   )
    #   upper_coef_plus = SVector( T(),
    #                              T(),
    #                              T(), )
    #   central_coef_plus = T()
    #   lower_coef_plus = SVector( T(),
    #                              T(),
    #                              T(), )
    #   left_weights = SVector( T(),
    #                           T(),
    #                           T(),
    #                           T(),
    #                           T(),
    #                           T(), )
    #   right_weights = left_weights
    #   left_boundary_derivatives = Tuple{}()
    #   right_boundary_derivatives = left_boundary_derivatives

    #   left_boundary_minus = (
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,7}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,8}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #   )
    #   right_boundary_minus = (
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,6}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,7}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,8}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,9}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #       DerivativeCoefficientRow{T,1,10}(SVector(T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(),
    #                                               T(), )),
    #   )
    #   upper_coef_minus     = .- lower_coef_plus
    #   central_coef_minus   = .- central_coef_plus
    #   lower_coef_minus     = .- upper_coef_plus

    #   left_boundary_central  = (left_boundary_plus  .+ left_boundary_minus)  ./ 2
    #   right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
    #   upper_coef_central     = widening_plus(upper_coef_plus, upper_coef_minus) / 2
    #   central_coef_central   = (central_coef_plus   + central_coef_minus) / 2
    #   lower_coef_central     = widening_plus(lower_coef_plus, lower_coef_minus) / 2

    #   if source.kind === :plus
    #     left_boundary  = left_boundary_plus
    #     right_boundary = right_boundary_plus
    #     upper_coef     = upper_coef_plus
    #     central_coef   = central_coef_plus
    #     lower_coef     = lower_coef_plus
    #   elseif source.kind === :minus
    #     left_boundary  = left_boundary_minus
    #     right_boundary = right_boundary_minus
    #     upper_coef     = upper_coef_minus
    #     central_coef   = central_coef_minus
    #     lower_coef     = lower_coef_minus
    #   elseif source.kind === :central
    #     left_boundary  = left_boundary_central
    #     right_boundary = right_boundary_central
    #     upper_coef     = upper_coef_central
    #     central_coef   = central_coef_central
    #     lower_coef     = lower_coef_central
    #   end
    #   DerivativeCoefficients(left_boundary, right_boundary,
    #                           left_boundary_derivatives, right_boundary_derivatives,
    #                           lower_coef, central_coef, upper_coef,
    #                           left_weights, right_weights, mode, 1, order, source)
  else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end


# function second_derivative_coefficients(source::Mattsson2017, order::Int, T=Float64, mode=FastMode())
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
#                                 left_weights, right_weights, mode, 2, order, source)
#     else
#         throw(ArgumentError("Order $order not implemented/derived."))
#     end
# end
