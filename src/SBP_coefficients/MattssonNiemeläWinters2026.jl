
"""
    MattssonNiemeläWinters2026(version::Symbol())

Coefficients of the optimized upwind SBP operators with nonuniform grid given in
- Mattsson, Niemelä, Winters (2026)
  Optimal boundary closures for diagonal-norm upwind SBP operators.
  [arXiv:2602.05727](https://arxiv.org/abs/2602.05727)

You can choose between the different versions `:central`, `:plus`, and `:minus`.
"""
struct MattssonNiemeläWinters2026 <: SourceOfCoefficients
    kind::Symbol

    function MattssonNiemeläWinters2026(kind::Symbol)
        if (kind !== :plus) && (kind !== :minus) && (kind !== :central)
            throw(ArgumentError("The only choices are :plus, :minus, and :central, not :$kind."))
        end
        new(kind)
    end
end

function Base.show(io::IO, source::MattssonNiemeläWinters2026)
    if get(io, :compact, false)
        summary(io, source)
    else
        print(io,
              "Mattsson, Niemelä, Winters (2026) \n",
              "  Optimal boundary closures for diagonal-norm upwind SBP operators. \n",
              "  arXiv e-print 2602.05727.")
    end
end

function construct_grid(source::MattssonNiemeläWinters2026, accuracy_order, xmin, xmax, N)
    T = promote_type(typeof(xmin), typeof(xmax))

    if accuracy_order == 2
        xstart = SVector(T(0.0000000000000e+00),
                         T(1.0158066888926694105),
                         T(2.0158066888926694105))
    elseif accuracy_order == 3
        xstart = SVector(T(0.0000000000000e+00),
                         T(0.78643399666738163986),
                         T(1.78643399666738163986))
    elseif accuracy_order == 4
        xstart = SVector(T(0.0000000000000e+00),
                         T(0.67776056835278144950),
                         T(1.67776056835278144950))
    elseif accuracy_order == 5
        xstart = SVector(T(0.0000000000000e+00),
                         T(0.62048159180300444152),
                         T(1.62048159180300444152))
    elseif accuracy_order == 6
        xstart = SVector(T(0.0000000000000e+00),
                         T(0.46488458913890803995),
                         T(1.33571258906181027793))
    elseif accuracy_order == 7
        xstart = SVector(T(0.0000000000000e+00),
                         T(0.47968913918690932952),
                         T(1.37532312518273333739))
    elseif accuracy_order == 8
        xstart = SVector(T(0.0000000000000e+00),
                         T(0.46931054871854486577),
                         T(1.37099127623205908804))
    elseif accuracy_order == 9
        xstart = SVector(T(0.0000000000000e+00),
                         T(0.46457087199875349857),
                         T(1.35931258539346483716))
    else
        throw(ArgumentError("Accuracy order $accuracy_order not implemented/derived."))
    end

    BoundaryAdaptedGrid(xmin, xmax, xstart, N)
end

function first_derivative_coefficients(source::MattssonNiemeläWinters2026,
                                       order::Int,
                                       T = Float64,
                                       mode = FastMode())
    if order == 2
        left_boundary_plus = (DerivativeCoefficientRow{T, 1, 3}(SVector(T(-1.672364656310962),
                                                                        T(2.371163860459328),
                                                                        T(-0.698799204148366))),
                              DerivativeCoefficientRow{T, 1, 4}(SVector(T(-0.525917296230532),
                                                                        T(-0.384308750879223),
                                                                        T(1.354682401534832),
                                                                        T(-0.444456354425077))),
                              DerivativeCoefficientRow{T, 1, 5}(SVector(T(0.128177867047527),
                                                                        T(-0.319438992586323),
                                                                        T(-1.217153088225320),
                                                                        T(1.877885618352155),
                                                                        T(-0.469471404588039))))
        right_boundary_plus = (DerivativeCoefficientRow{T, 1, 3}(SVector(T(1.396915985042026),
                                                                         T(-1.815912586468807),
                                                                         T(0.418996601426782))),
                               DerivativeCoefficientRow{T, 1, 3}(SVector(T(0.686726935924425),
                                                                         T(-0.384308750879223),
                                                                         T(-0.302418185045202))),
                               DerivativeCoefficientRow{T, 1, 3}(SVector(T(-0.213774028660944),
                                                                         T(1.430927116886264),
                                                                         T(-1.217153088225320))))
        upper_coef_plus = SVector(T(2),
                                  T(-1 // 2))
        central_coef_plus = T(-3 // 2)
        lower_coef_plus = SVector{0, T}()
        left_weights = SVector(T(0.325809242246152),
                               T(1.124969853669369),
                               T(1.065027592977149))
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        left_boundary_minus = (DerivativeCoefficientRow{T, 1, 3}(SVector(T(-1.396915985042026),
                                                                         T(1.815912586468807),
                                                                         T(-0.418996601426782))),
                               DerivativeCoefficientRow{T, 1, 3}(SVector(T(-0.686726935924425),
                                                                         T(0.384308750879223),
                                                                         T(0.302418185045202))),
                               DerivativeCoefficientRow{T, 1, 3}(SVector(T(0.213774028660944),
                                                                         T(-1.430927116886264),
                                                                         T(1.217153088225320))))
        right_boundary_minus = (DerivativeCoefficientRow{T, 1, 3}(SVector(T(1.672364656310962),
                                                                          T(-2.371163860459328),
                                                                          T(0.698799204148366))),
                                DerivativeCoefficientRow{T, 1, 4}(SVector(T(0.525917296230532),
                                                                          T(0.384308750879223),
                                                                          T(-1.354682401534832),
                                                                          T(0.444456354425077))),
                                DerivativeCoefficientRow{T, 1, 5}(SVector(T(-0.128177867047527),
                                                                          T(0.319438992586323),
                                                                          T(1.217153088225320),
                                                                          T(-1.877885618352155),
                                                                          T(0.469471404588039))))
        upper_coef_minus = .-lower_coef_plus
        central_coef_minus = .-central_coef_plus
        lower_coef_minus = .-upper_coef_plus

        left_boundary_central = (left_boundary_plus .+ left_boundary_minus) ./ 2
        right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
        upper_coef_central = widening_plus(upper_coef_plus, upper_coef_minus) / 2
        central_coef_central = (central_coef_plus + central_coef_minus) / 2
        lower_coef_central = widening_plus(lower_coef_plus, lower_coef_minus) / 2

        if source.kind === :plus
            left_boundary = left_boundary_plus
            right_boundary = right_boundary_plus
            upper_coef = upper_coef_plus
            central_coef = central_coef_plus
            lower_coef = lower_coef_plus
        elseif source.kind === :minus
            left_boundary = left_boundary_minus
            right_boundary = right_boundary_minus
            upper_coef = upper_coef_minus
            central_coef = central_coef_minus
            lower_coef = lower_coef_minus
        elseif source.kind === :central
            left_boundary = left_boundary_central
            right_boundary = right_boundary_central
            upper_coef = upper_coef_central
            central_coef = central_coef_central
            lower_coef = lower_coef_central
        end
        DerivativeCoefficients(left_boundary, right_boundary,
                               left_boundary_derivatives, right_boundary_derivatives,
                               lower_coef, central_coef, upper_coef,
                               left_weights, right_weights, mode, 1, order, source)
    elseif order == 3
        left_boundary_plus = (DerivativeCoefficientRow{T, 1, 3}(SVector(T(-1.891809941491256),
                                                                        T(2.379593594713310),
                                                                        T(-0.487783653222054))),
                              DerivativeCoefficientRow{T, 1, 4}(SVector(T(-0.574873939656967),
                                                                        T(-0.135833228671310),
                                                                        T(0.873514746600907),
                                                                        T(-0.162807578272630))),
                              DerivativeCoefficientRow{T, 1, 5}(SVector(T(0.107116362952890),
                                                                        T(-0.517414050013166),
                                                                        T(-0.432130140896273),
                                                                        T(1.010913393547859),
                                                                        T(-0.168485565591310))))
        right_boundary_plus = (DerivativeCoefficientRow{T, 1, 3}(SVector(T(1.764146924242998),
                                                                         T(-2.151532040583887),
                                                                         T(0.387385116340889))),
                               DerivativeCoefficientRow{T, 1, 3}(SVector(T(0.635810352238159),
                                                                         T(-0.135833228671310),
                                                                         T(-0.499977123566849))),
                               DerivativeCoefficientRow{T, 1, 4}(SVector(T(-0.134877693119841),
                                                                         T(0.903978965198734),
                                                                         T(-0.432130140896273),
                                                                         T(-0.336971131182620))))
        upper_coef_plus = SVector(T(1),
                                  T(-1 // 6))
        central_coef_plus = T(-1 // 2)
        lower_coef_plus = SVector(T(-1 // 3))
        left_weights = SVector(T(0.273526203050309),
                               T(1.023703370782746),
                               T(0.989204422834326))
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        left_boundary_minus = (DerivativeCoefficientRow{T, 1, 3}(SVector(T(-1.764146924242998),
                                                                         T(2.151532040583887),
                                                                         T(-0.387385116340889))),
                               DerivativeCoefficientRow{T, 1, 3}(SVector(T(-0.635810352238159),
                                                                         T(0.135833228671310),
                                                                         T(0.499977123566849))),
                               DerivativeCoefficientRow{T, 1, 4}(SVector(T(0.134877693119841),
                                                                         T(-0.903978965198734),
                                                                         T(0.432130140896273),
                                                                         T(0.336971131182620))))
        right_boundary_minus = (DerivativeCoefficientRow{T, 1, 3}(SVector(T(1.891809941491256),
                                                                          T(-2.379593594713310),
                                                                          T(0.487783653222054))),
                                DerivativeCoefficientRow{T, 1, 4}(SVector(T(0.574873939656967),
                                                                          T(0.135833228671310),
                                                                          T(-0.873514746600907),
                                                                          T(0.162807578272630))),
                                DerivativeCoefficientRow{T, 1, 5}(SVector(T(-0.107116362952890),
                                                                          T(0.517414050013166),
                                                                          T(0.432130140896273),
                                                                          T(-1.010913393547859),
                                                                          T(0.168485565591310))))
        upper_coef_minus = .-lower_coef_plus
        central_coef_minus = .-central_coef_plus
        lower_coef_minus = .-upper_coef_plus

        left_boundary_central = (left_boundary_plus .+ left_boundary_minus) ./ 2
        right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
        upper_coef_central = widening_plus(upper_coef_plus, upper_coef_minus) / 2
        central_coef_central = (central_coef_plus + central_coef_minus) / 2
        lower_coef_central = widening_plus(lower_coef_plus, lower_coef_minus) / 2

        if source.kind === :plus
            left_boundary = left_boundary_plus
            right_boundary = right_boundary_plus
            upper_coef = upper_coef_plus
            central_coef = central_coef_plus
            lower_coef = lower_coef_plus
        elseif source.kind === :minus
            left_boundary = left_boundary_minus
            right_boundary = right_boundary_minus
            upper_coef = upper_coef_minus
            central_coef = central_coef_minus
            lower_coef = lower_coef_minus
        elseif source.kind === :central
            left_boundary = left_boundary_central
            right_boundary = right_boundary_central
            upper_coef = upper_coef_central
            central_coef = central_coef_central
            lower_coef = lower_coef_central
        end
        DerivativeCoefficients(left_boundary, right_boundary,
                               left_boundary_derivatives, right_boundary_derivatives,
                               lower_coef, central_coef, upper_coef,
                               left_weights, right_weights, mode, 1, order, source)
    elseif order == 4
        left_boundary_plus = (DerivativeCoefficientRow{T, 1, 4}(SVector(T(-2.556358917106399),
                                                                        T(3.564640991464135),
                                                                        T(-1.283964876143738),
                                                                        T(0.275682801786002))),
                              DerivativeCoefficientRow{T, 1, 5}(SVector(T(-0.649346730903335),
                                                                        T(-0.125078806540877),
                                                                        T(1.072672394278584),
                                                                        T(-0.381966566829468),
                                                                        T(0.083719709995096))),
                              DerivativeCoefficientRow{T, 1, 6}(SVector(T(0.205936021419767),
                                                                        T(-0.704908069660006),
                                                                        T(-0.485218213810784),
                                                                        T(1.413674006262887),
                                                                        T(-0.515380493054236),
                                                                        T(0.085896748842373))),
                              DerivativeCoefficientRow{T, 1, 7}(SVector(T(-0.036650201684175),
                                                                        T(0.087379701967370),
                                                                        T(-0.334196340376952),
                                                                        T(-0.788921357978960),
                                                                        T(1.484845197331454),
                                                                        T(-0.494948399110485),
                                                                        T(0.082491399851747))))
        right_boundary_plus = (DerivativeCoefficientRow{T, 1, 4}(SVector(T(2.393833280243093),
                                                                         T(-3.199556314039223),
                                                                         T(0.989000263032899),
                                                                         T(-0.183277229236769))),
                               DerivativeCoefficientRow{T, 1, 4}(SVector(T(0.723440298423478),
                                                                         T(-0.125078806540877),
                                                                         T(-0.687042291594004),
                                                                         T(0.088680799711403))),
                               DerivativeCoefficientRow{T, 1, 4}(SVector(T(-0.267355457949934),
                                                                         T(1.100566058421508),
                                                                         T(-0.485218213810784),
                                                                         T(-0.347992386660790))),
                               DerivativeCoefficientRow{T, 1, 5}(SVector(T(0.055128672167248),
                                                                         T(-0.376362469437301),
                                                                         T(1.357629354804255),
                                                                         T(-0.788921357978960),
                                                                         T(-0.247474199555242))))
        upper_coef_plus = SVector(T(3 // 2),
                                  T(-1 // 2),
                                  T(1 // 12))
        central_coef_plus = T(-5 // 6)
        lower_coef_plus = SVector(T(-1 // 4))
        left_weights = SVector(T(0.202012358335387),
                               T(0.995384878163280),
                               T(0.970157013582162),
                               T(1.010206318271953))
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        left_boundary_minus = (DerivativeCoefficientRow{T, 1, 4}(SVector(T(-2.393833280243093),
                                                                         T(3.199556314039223),
                                                                         T(-0.989000263032899),
                                                                         T(0.183277229236769))),
                               DerivativeCoefficientRow{T, 1, 4}(SVector(T(-0.723440298423478),
                                                                         T(0.125078806540877),
                                                                         T(0.687042291594004),
                                                                         T(-0.088680799711403))),
                               DerivativeCoefficientRow{T, 1, 4}(SVector(T(0.267355457949934),
                                                                         T(-1.100566058421508),
                                                                         T(0.485218213810784),
                                                                         T(0.347992386660790))),
                               DerivativeCoefficientRow{T, 1, 5}(SVector(T(-0.055128672167248),
                                                                         T(0.376362469437301),
                                                                         T(-1.357629354804255),
                                                                         T(0.788921357978960),
                                                                         T(0.247474199555242))))
        right_boundary_minus = (DerivativeCoefficientRow{T, 1, 4}(SVector(T(2.556358917106399),
                                                                          T(-3.564640991464135),
                                                                          T(1.283964876143738),
                                                                          T(-0.275682801786002))),
                                DerivativeCoefficientRow{T, 1, 5}(SVector(T(0.649346730903335),
                                                                          T(0.125078806540877),
                                                                          T(-1.072672394278584),
                                                                          T(0.381966566829468),
                                                                          T(-0.083719709995096))),
                                DerivativeCoefficientRow{T, 1, 6}(SVector(T(-0.205936021419767),
                                                                          T(0.704908069660006),
                                                                          T(0.485218213810784),
                                                                          T(-1.413674006262887),
                                                                          T(0.515380493054236),
                                                                          T(-0.085896748842373))),
                                DerivativeCoefficientRow{T, 1, 7}(SVector(T(0.036650201684175),
                                                                          T(-0.087379701967370),
                                                                          T(0.334196340376952),
                                                                          T(0.788921357978960),
                                                                          T(-1.484845197331454),
                                                                          T(0.494948399110485),
                                                                          T(-0.082491399851747))))
        upper_coef_minus = .-lower_coef_plus
        central_coef_minus = .-central_coef_plus
        lower_coef_minus = .-upper_coef_plus

        left_boundary_central = (left_boundary_plus .+ left_boundary_minus) ./ 2
        right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
        upper_coef_central = widening_plus(upper_coef_plus, upper_coef_minus) / 2
        central_coef_central = (central_coef_plus + central_coef_minus) / 2
        lower_coef_central = widening_plus(lower_coef_plus, lower_coef_minus) / 2

        if source.kind === :plus
            left_boundary = left_boundary_plus
            right_boundary = right_boundary_plus
            upper_coef = upper_coef_plus
            central_coef = central_coef_plus
            lower_coef = lower_coef_plus
        elseif source.kind === :minus
            left_boundary = left_boundary_minus
            right_boundary = right_boundary_minus
            upper_coef = upper_coef_minus
            central_coef = central_coef_minus
            lower_coef = lower_coef_minus
        elseif source.kind === :central
            left_boundary = left_boundary_central
            right_boundary = right_boundary_central
            upper_coef = upper_coef_central
            central_coef = central_coef_central
            lower_coef = lower_coef_central
        end
        DerivativeCoefficients(left_boundary, right_boundary,
                               left_boundary_derivatives, right_boundary_derivatives,
                               lower_coef, central_coef, upper_coef,
                               left_weights, right_weights, mode, 1, order, source)
    elseif order == 5
        left_boundary_plus = (DerivativeCoefficientRow{T, 1, 4}(SVector(T(-2.654846183131157),
                                                                        T(3.516343816487329),
                                                                        T(-1.075712081010993),
                                                                        T(0.214214447654821))),
                              DerivativeCoefficientRow{T, 1, 5}(SVector(T(-0.700183010215649),
                                                                        T(-0.049626684071346),
                                                                        T(0.970340075938412),
                                                                        T(-0.256800400283813),
                                                                        T(0.036270018632396))),
                              DerivativeCoefficientRow{T, 1, 6}(SVector(T(0.182096705039648),
                                                                        T(-0.738439255266119),
                                                                        T(-0.181426486663972),
                                                                        T(0.951824594913516),
                                                                        T(-0.246987182334315),
                                                                        T(0.032931624311242))),
                              DerivativeCoefficientRow{T, 1, 7}(SVector(T(-0.034458638496284),
                                                                        T(0.122486191317302),
                                                                        T(-0.553771332521901),
                                                                        T(-0.318768238242615),
                                                                        T(1.001504703757659),
                                                                        T(-0.250376175939415),
                                                                        T(0.033383490125255))))
        right_boundary_plus = (DerivativeCoefficientRow{T, 1, 4}(SVector(T(2.587533828965011),
                                                                         T(-3.373424806464659),
                                                                         T(0.966264845958962),
                                                                         T(-0.180373868459314))),
                               DerivativeCoefficientRow{T, 1, 4}(SVector(T(0.729847066299824),
                                                                         T(-0.049626684071346),
                                                                         T(-0.813297433927422),
                                                                         T(0.133077051698945))),
                               DerivativeCoefficientRow{T, 1, 5}(SVector(T(-0.202722500298602),
                                                                         T(0.881027251703800),
                                                                         T(-0.181426486663972),
                                                                         T(-0.546275701208089),
                                                                         T(0.049397436466863))),
                               DerivativeCoefficientRow{T, 1, 6}(SVector(T(0.040923545497300),
                                                                         T(-0.236363088586311),
                                                                         T(0.964884898022573),
                                                                         T(-0.318768238242615),
                                                                         T(-0.500752351878829),
                                                                         T(0.050075235187883))))
        upper_coef_plus = SVector(T(1),
                                  T(-1 // 4),
                                  T(1 // 30))
        central_coef_plus = T(-1 // 3)
        lower_coef_plus = SVector(T(-1 // 2),
                                  T(1 // 20))
        left_weights = SVector(T(0.190753054470034),
                               T(0.919032704977988),
                               T(1.012198275380973),
                               T(0.998497556974009))
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        left_boundary_minus = (DerivativeCoefficientRow{T, 1, 4}(SVector(T(-2.587533828965011),
                                                                         T(3.373424806464659),
                                                                         T(-0.966264845958962),
                                                                         T(0.180373868459314))),
                               DerivativeCoefficientRow{T, 1, 4}(SVector(T(-0.729847066299824),
                                                                         T(0.049626684071346),
                                                                         T(0.813297433927422),
                                                                         T(-0.133077051698945))),
                               DerivativeCoefficientRow{T, 1, 5}(SVector(T(0.202722500298602),
                                                                         T(-0.881027251703800),
                                                                         T(0.181426486663972),
                                                                         T(0.546275701208089),
                                                                         T(-0.049397436466863))),
                               DerivativeCoefficientRow{T, 1, 6}(SVector(T(-0.040923545497300),
                                                                         T(0.236363088586311),
                                                                         T(-0.964884898022573),
                                                                         T(0.318768238242615),
                                                                         T(0.500752351878829),
                                                                         T(-0.050075235187883))))
        right_boundary_minus = (DerivativeCoefficientRow{T, 1, 4}(SVector(T(2.654846183131157),
                                                                          T(-3.516343816487329),
                                                                          T(1.075712081010993),
                                                                          T(-0.214214447654821))),
                                DerivativeCoefficientRow{T, 1, 5}(SVector(T(0.700183010215649),
                                                                          T(0.049626684071346),
                                                                          T(-0.970340075938412),
                                                                          T(0.256800400283813),
                                                                          T(-0.036270018632396))),
                                DerivativeCoefficientRow{T, 1, 6}(SVector(T(-0.182096705039648),
                                                                          T(0.738439255266119),
                                                                          T(0.181426486663972),
                                                                          T(-0.951824594913516),
                                                                          T(0.246987182334315),
                                                                          T(-0.032931624311242))),
                                DerivativeCoefficientRow{T, 1, 7}(SVector(T(0.034458638496284),
                                                                          T(-0.122486191317302),
                                                                          T(0.553771332521901),
                                                                          T(0.318768238242615),
                                                                          T(-1.001504703757659),
                                                                          T(0.250376175939415),
                                                                          T(-0.033383490125255))))
        upper_coef_minus = .-lower_coef_plus
        central_coef_minus = .-central_coef_plus
        lower_coef_minus = .-upper_coef_plus

        left_boundary_central = (left_boundary_plus .+ left_boundary_minus) ./ 2
        right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
        upper_coef_central = widening_plus(upper_coef_plus, upper_coef_minus) / 2
        central_coef_central = (central_coef_plus + central_coef_minus) / 2
        lower_coef_central = widening_plus(lower_coef_plus, lower_coef_minus) / 2

        if source.kind === :plus
            left_boundary = left_boundary_plus
            right_boundary = right_boundary_plus
            upper_coef = upper_coef_plus
            central_coef = central_coef_plus
            lower_coef = lower_coef_plus
        elseif source.kind === :minus
            left_boundary = left_boundary_minus
            right_boundary = right_boundary_minus
            upper_coef = upper_coef_minus
            central_coef = central_coef_minus
            lower_coef = lower_coef_minus
        elseif source.kind === :central
            left_boundary = left_boundary_central
            right_boundary = right_boundary_central
            upper_coef = upper_coef_central
            central_coef = central_coef_central
            lower_coef = lower_coef_central
        end
        DerivativeCoefficients(left_boundary, right_boundary,
                               left_boundary_derivatives, right_boundary_derivatives,
                               lower_coef, central_coef, upper_coef,
                               left_weights, right_weights, mode, 1, order, source)
    elseif order == 6
        left_boundary_plus = (DerivativeCoefficientRow{T, 1, 6}(SVector(T(-3.755677079270491),
                                                                        T(5.094616935299430),
                                                                        T(-1.904387215859687),
                                                                        T(0.729129749158562),
                                                                        T(-0.181947339732544),
                                                                        T(0.018264950404730))),
                              DerivativeCoefficientRow{T, 1, 6}(SVector(T(-0.906195140897546),
                                                                        T(-0.049610433237538),
                                                                        T(1.319719117473809),
                                                                        T(-0.501327659814975),
                                                                        T(0.164533101808823),
                                                                        T(-0.027118985332573))),
                              DerivativeCoefficientRow{T, 1, 7}(SVector(T(0.237732085727226),
                                                                        T(-0.883548218164051),
                                                                        T(-0.127655240294382),
                                                                        T(1.090339416028211),
                                                                        T(-0.425677904826528),
                                                                        T(0.126099840422713),
                                                                        T(-0.017289978893190))),
                              DerivativeCoefficientRow{T, 1, 8}(SVector(T(-0.078345651107586),
                                                                        T(0.253363947264913),
                                                                        T(-0.704237731113124),
                                                                        T(-0.310852524115995),
                                                                        T(1.192205144457065),
                                                                        T(-0.467435946763637),
                                                                        T(0.131774584432416),
                                                                        T(-0.016471823054052))),
                              DerivativeCoefficientRow{T, 1, 9}(SVector(T(0.017228501873662),
                                                                        T(-0.060815316062831),
                                                                        T(0.146499870422555),
                                                                        T(-0.525347168840596),
                                                                        T(-0.512611118374119),
                                                                        T(1.319975450567247),
                                                                        T(-0.502082895112066),
                                                                        T(0.133888772029884),
                                                                        T(-0.016736096503736))),
                              DerivativeCoefficientRow{T, 1, 10}(SVector(T(-0.001304824701429),
                                                                         T(0.006644310538497),
                                                                         T(-0.015912211601629),
                                                                         T(0.053919312875515),
                                                                         T(-0.413969955689860),
                                                                         T(-0.578702616364331),
                                                                         T(1.332387347288753),
                                                                         T(-0.499645255233282),
                                                                         T(0.133238734728875),
                                                                         T(-0.016654841841109))))
        right_boundary_plus = (DerivativeCoefficientRow{T, 1, 6}(SVector(T(3.677305363325555),
                                                                         T(-4.909252416542798),
                                                                         T(1.703355095812326),
                                                                         T(-0.589230314929821),
                                                                         T(0.127527897470572),
                                                                         T(-0.009705625135835))),
                               DerivativeCoefficientRow{T, 1, 6}(SVector(T(0.940411435343120),
                                                                         T(-0.049610433237538),
                                                                         T(-1.168568069849456),
                                                                         T(0.351739685798358),
                                                                         T(-0.083095407099030),
                                                                         T(0.009122789044547))),
                               DerivativeCoefficientRow{T, 1, 6}(SVector(T(-0.265789526782542),
                                                                         T(0.997832736325951),
                                                                         T(-0.127655240294382),
                                                                         T(-0.739217235807965),
                                                                         T(0.151348294800747),
                                                                         T(-0.016519028241809))),
                               DerivativeCoefficientRow{T, 1, 6}(SVector(T(0.096947057020551),
                                                                         T(-0.361114653512880),
                                                                         T(1.038744930842483),
                                                                         T(-0.310852524115995),
                                                                         T(-0.517051607891842),
                                                                         T(0.053326797657683))),
                               DerivativeCoefficientRow{T, 1, 7}(SVector(T(-0.024580347874183),
                                                                         T(0.120417396564127),
                                                                         T(-0.412041364462887),
                                                                         T(1.211332849096812),
                                                                         T(-0.512611118374119),
                                                                         T(-0.415989607957220),
                                                                         T(0.033472193007471))),
                               DerivativeCoefficientRow{T, 1, 8}(SVector(T(0.002455540794633),
                                                                         T(-0.019751301839679),
                                                                         T(0.121467638069620),
                                                                         T(-0.472629637815500),
                                                                         T(1.313566897659665),
                                                                         T(-0.578702616364331),
                                                                         T(-0.399716204186626),
                                                                         T(0.033309683682219))))
        upper_coef_plus = SVector(T(4 // 3),
                                  T(-1 // 2),
                                  T(2 // 15),
                                  T(-1 // 60))
        central_coef_plus = T(-7 // 12)
        lower_coef_plus = SVector(T(-2 // 5),
                                  T(1 // 30))
        left_weights = SVector(T(0.134535498734575),
                               T(0.728837192416789),
                               T(0.963949509113091),
                               T(1.011828903939489),
                               T(0.995851491591640),
                               T(1.000709993266226))
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        left_boundary_minus = (DerivativeCoefficientRow{T, 1, 6}(SVector(T(-3.677305363325555),
                                                                         T(4.909252416542798),
                                                                         T(-1.703355095812326),
                                                                         T(0.589230314929821),
                                                                         T(-0.127527897470572),
                                                                         T(0.009705625135835))),
                               DerivativeCoefficientRow{T, 1, 6}(SVector(T(-0.940411435343120),
                                                                         T(0.049610433237538),
                                                                         T(1.168568069849456),
                                                                         T(-0.351739685798358),
                                                                         T(0.083095407099030),
                                                                         T(-0.009122789044547))),
                               DerivativeCoefficientRow{T, 1, 6}(SVector(T(0.265789526782542),
                                                                         T(-0.997832736325951),
                                                                         T(0.127655240294382),
                                                                         T(0.739217235807965),
                                                                         T(-0.151348294800747),
                                                                         T(0.016519028241809))),
                               DerivativeCoefficientRow{T, 1, 6}(SVector(T(-0.096947057020551),
                                                                         T(0.361114653512880),
                                                                         T(-1.038744930842483),
                                                                         T(0.310852524115995),
                                                                         T(0.517051607891842),
                                                                         T(-0.053326797657683))),
                               DerivativeCoefficientRow{T, 1, 7}(SVector(T(0.024580347874183),
                                                                         T(-0.120417396564127),
                                                                         T(0.412041364462887),
                                                                         T(-1.211332849096812),
                                                                         T(0.512611118374119),
                                                                         T(0.415989607957220),
                                                                         T(-0.033472193007471))),
                               DerivativeCoefficientRow{T, 1, 8}(SVector(T(-0.002455540794633),
                                                                         T(0.019751301839679),
                                                                         T(-0.121467638069620),
                                                                         T(0.472629637815500),
                                                                         T(-1.313566897659665),
                                                                         T(0.578702616364331),
                                                                         T(0.399716204186626),
                                                                         T(-0.033309683682219))))
        right_boundary_minus = (DerivativeCoefficientRow{T, 1, 6}(SVector(T(3.755677079270491),
                                                                          T(-5.094616935299430),
                                                                          T(1.904387215859687),
                                                                          T(-0.729129749158562),
                                                                          T(0.181947339732544),
                                                                          T(-0.018264950404730))),
                                DerivativeCoefficientRow{T, 1, 6}(SVector(T(0.906195140897546),
                                                                          T(0.049610433237538),
                                                                          T(-1.319719117473809),
                                                                          T(0.501327659814975),
                                                                          T(-0.164533101808823),
                                                                          T(0.027118985332573))),
                                DerivativeCoefficientRow{T, 1, 7}(SVector(T(-0.237732085727226),
                                                                          T(0.883548218164051),
                                                                          T(0.127655240294382),
                                                                          T(-1.090339416028211),
                                                                          T(0.425677904826528),
                                                                          T(-0.126099840422713),
                                                                          T(0.017289978893190))),
                                DerivativeCoefficientRow{T, 1, 8}(SVector(T(0.078345651107586),
                                                                          T(-0.253363947264913),
                                                                          T(0.704237731113124),
                                                                          T(0.310852524115995),
                                                                          T(-1.192205144457065),
                                                                          T(0.467435946763637),
                                                                          T(-0.131774584432416),
                                                                          T(0.016471823054052))),
                                DerivativeCoefficientRow{T, 1, 9}(SVector(T(-0.017228501873662),
                                                                          T(0.060815316062831),
                                                                          T(-0.146499870422555),
                                                                          T(0.525347168840596),
                                                                          T(0.512611118374119),
                                                                          T(-1.319975450567247),
                                                                          T(0.502082895112066),
                                                                          T(-0.133888772029884),
                                                                          T(0.016736096503736))),
                                DerivativeCoefficientRow{T, 1, 10}(SVector(T(0.001304824701429),
                                                                           T(-0.006644310538497),
                                                                           T(0.015912211601629),
                                                                           T(-0.053919312875515),
                                                                           T(0.413969955689860),
                                                                           T(0.578702616364331),
                                                                           T(-1.332387347288753),
                                                                           T(0.499645255233282),
                                                                           T(-0.133238734728875),
                                                                           T(0.016654841841109))))
        upper_coef_minus = .-lower_coef_plus
        central_coef_minus = .-central_coef_plus
        lower_coef_minus = .-upper_coef_plus

        left_boundary_central = (left_boundary_plus .+ left_boundary_minus) ./ 2
        right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
        upper_coef_central = widening_plus(upper_coef_plus, upper_coef_minus) / 2
        central_coef_central = (central_coef_plus + central_coef_minus) / 2
        lower_coef_central = widening_plus(lower_coef_plus, lower_coef_minus) / 2

        if source.kind === :plus
            left_boundary = left_boundary_plus
            right_boundary = right_boundary_plus
            upper_coef = upper_coef_plus
            central_coef = central_coef_plus
            lower_coef = lower_coef_plus
        elseif source.kind === :minus
            left_boundary = left_boundary_minus
            right_boundary = right_boundary_minus
            upper_coef = upper_coef_minus
            central_coef = central_coef_minus
            lower_coef = lower_coef_minus
        elseif source.kind === :central
            left_boundary = left_boundary_central
            right_boundary = right_boundary_central
            upper_coef = upper_coef_central
            central_coef = central_coef_central
            lower_coef = lower_coef_central
        end
        DerivativeCoefficients(left_boundary, right_boundary,
                               left_boundary_derivatives, right_boundary_derivatives,
                               lower_coef, central_coef, upper_coef,
                               left_weights, right_weights, mode, 1, order, source)
    elseif order == 7
        left_boundary_plus = (DerivativeCoefficientRow{T, 1, 6}(SVector(T(-3.574597921569969),
                                                                        T(4.786324846066822),
                                                                        T(-1.688767520049126),
                                                                        T(0.590379354802345),
                                                                        T(-0.120204938336236),
                                                                        T(0.006866179086163))),
                              DerivativeCoefficientRow{T, 1, 6}(SVector(T(-0.883923454973464),
                                                                        T(-0.021028659504377),
                                                                        T(1.217132668682104),
                                                                        T(-0.401478932917210),
                                                                        T(0.100930553870111),
                                                                        T(-0.011632175157163))),
                              DerivativeCoefficientRow{T, 1, 7}(SVector(T(0.228437011813396),
                                                                        T(-0.875921961544600),
                                                                        T(-0.054450212538879),
                                                                        T(0.924153268224264),
                                                                        T(-0.279422110582068),
                                                                        T(0.064455838416543),
                                                                        T(-0.007251833788656))),
                              DerivativeCoefficientRow{T, 1, 8}(SVector(T(-0.073951206233434),
                                                                        T(0.252161137373666),
                                                                        T(-0.761614155953779),
                                                                        T(-0.136498318688696),
                                                                        T(0.951375100888152),
                                                                        T(-0.290907863728331),
                                                                        T(0.066567543103514),
                                                                        T(-0.007132236761091))),
                              DerivativeCoefficientRow{T, 1, 9}(SVector(T(0.013343770541662),
                                                                        T(-0.049601736079086),
                                                                        T(0.160649297156179),
                                                                        T(-0.657837705860311),
                                                                        T(-0.219550214526392),
                                                                        T(0.993510363631668),
                                                                        T(-0.300046887453750),
                                                                        T(0.066677086100833),
                                                                        T(-0.007143973510804))),
                              DerivativeCoefficientRow{T, 1, 10}(SVector(T(-0.000391462486334),
                                                                         T(0.003009678529462),
                                                                         T(-0.017964428161274),
                                                                         T(0.111241495585250),
                                                                         T(-0.607178485175455),
                                                                         T(-0.248233272383674),
                                                                         T(0.999990342064734),
                                                                         T(-0.299997102619420),
                                                                         T(0.066666022804316),
                                                                         T(-0.007142788157605))))
        right_boundary_plus = (DerivativeCoefficientRow{T, 1, 6}(SVector(T(3.542121045611527),
                                                                         T(-4.708502023508660),
                                                                         T(1.601291541420409),
                                                                         T(-0.527073633296073),
                                                                         T(0.094949025180805),
                                                                         T(-0.002785955408008))),
                               DerivativeCoefficientRow{T, 1, 6}(SVector(T(0.898533073456783),
                                                                         T(-0.021028659504377),
                                                                         T(-1.152660076088575),
                                                                         T(0.337393009740077),
                                                                         T(-0.066258365606224),
                                                                         T(0.004021018002315))),
                               DerivativeCoefficientRow{T, 1, 6}(SVector(T(-0.240916157956684),
                                                                         T(0.924915555529416),
                                                                         T(-0.054450212538879),
                                                                         T(-0.774385295254804),
                                                                         T(0.163074792967703),
                                                                         T(-0.018238682746751))),
                               DerivativeCoefficientRow{T, 1, 7}(SVector(T(0.082833332318153),
                                                                         T(-0.300057741071639),
                                                                         T(0.908912159959023),
                                                                         T(-0.136498318688696),
                                                                         T(-0.656756952062214),
                                                                         T(0.111077168560161),
                                                                         T(-0.009509649014788))),
                               DerivativeCoefficientRow{T, 1, 8}(SVector(T(-0.016893139366928),
                                                                         T(0.075557714857230),
                                                                         T(-0.275266120888442),
                                                                         T(0.952940675870603),
                                                                         T(-0.219550214526392),
                                                                         T(-0.607279247082916),
                                                                         T(0.100015629151250),
                                                                         T(-0.009525298014405))),
                               DerivativeCoefficientRow{T, 1, 9}(SVector(T(0.000964786273662),
                                                                         T(-0.008706528496340),
                                                                         T(0.063486617695290),
                                                                         T(-0.291338231412718),
                                                                         T(0.993345516899620),
                                                                         T(-0.248233272383674),
                                                                         T(-0.599994205238841),
                                                                         T(0.099999034206473),
                                                                         T(-0.009523717543474))))
        upper_coef_plus = SVector(T(1),
                                  T(-3 // 10),
                                  T(1 // 15),
                                  T(-1 // 140))
        central_coef_plus = T(-1 // 4)
        lower_coef_plus = SVector(T(-3 // 5),
                                  T(1 // 10),
                                  T(-1 // 105))
        left_weights = SVector(T(0.140514189841058),
                               T(0.748493937428286),
                               T(0.984972539501733),
                               T(1.001489067472393),
                               T(0.999843732910721),
                               T(1.000009658028542))
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        left_boundary_minus = (DerivativeCoefficientRow{T, 1, 6}(SVector(T(-3.542121045611527),
                                                                         T(4.708502023508660),
                                                                         T(-1.601291541420409),
                                                                         T(0.527073633296073),
                                                                         T(-0.094949025180805),
                                                                         T(0.002785955408008))),
                               DerivativeCoefficientRow{T, 1, 6}(SVector(T(-0.898533073456783),
                                                                         T(0.021028659504377),
                                                                         T(1.152660076088575),
                                                                         T(-0.337393009740077),
                                                                         T(0.066258365606224),
                                                                         T(-0.004021018002315))),
                               DerivativeCoefficientRow{T, 1, 6}(SVector(T(0.240916157956684),
                                                                         T(-0.924915555529416),
                                                                         T(0.054450212538879),
                                                                         T(0.774385295254804),
                                                                         T(-0.163074792967703),
                                                                         T(0.018238682746751))),
                               DerivativeCoefficientRow{T, 1, 7}(SVector(T(-0.082833332318153),
                                                                         T(0.300057741071639),
                                                                         T(-0.908912159959023),
                                                                         T(0.136498318688696),
                                                                         T(0.656756952062214),
                                                                         T(-0.111077168560161),
                                                                         T(0.009509649014788))),
                               DerivativeCoefficientRow{T, 1, 8}(SVector(T(0.016893139366928),
                                                                         T(-0.075557714857230),
                                                                         T(0.275266120888442),
                                                                         T(-0.952940675870603),
                                                                         T(0.219550214526392),
                                                                         T(0.607279247082916),
                                                                         T(-0.100015629151250),
                                                                         T(0.009525298014405))),
                               DerivativeCoefficientRow{T, 1, 9}(SVector(T(-0.000964786273662),
                                                                         T(0.008706528496340),
                                                                         T(-0.063486617695290),
                                                                         T(0.291338231412718),
                                                                         T(-0.993345516899620),
                                                                         T(0.248233272383674),
                                                                         T(0.599994205238841),
                                                                         T(-0.099999034206473),
                                                                         T(0.009523717543474))))
        right_boundary_minus = (DerivativeCoefficientRow{T, 1, 6}(SVector(T(3.574597921569969),
                                                                          T(-4.786324846066822),
                                                                          T(1.688767520049126),
                                                                          T(-0.590379354802345),
                                                                          T(0.120204938336236),
                                                                          T(-0.006866179086163))),
                                DerivativeCoefficientRow{T, 1, 6}(SVector(T(0.883923454973464),
                                                                          T(0.021028659504377),
                                                                          T(-1.217132668682104),
                                                                          T(0.401478932917210),
                                                                          T(-0.100930553870111),
                                                                          T(0.011632175157163))),
                                DerivativeCoefficientRow{T, 1, 7}(SVector(T(-0.228437011813396),
                                                                          T(0.875921961544600),
                                                                          T(0.054450212538879),
                                                                          T(-0.924153268224264),
                                                                          T(0.279422110582068),
                                                                          T(-0.064455838416543),
                                                                          T(0.007251833788656))),
                                DerivativeCoefficientRow{T, 1, 8}(SVector(T(0.073951206233434),
                                                                          T(-0.252161137373666),
                                                                          T(0.761614155953779),
                                                                          T(0.136498318688696),
                                                                          T(-0.951375100888152),
                                                                          T(0.290907863728331),
                                                                          T(-0.066567543103514),
                                                                          T(0.007132236761091))),
                                DerivativeCoefficientRow{T, 1, 9}(SVector(T(-0.013343770541662),
                                                                          T(0.049601736079086),
                                                                          T(-0.160649297156179),
                                                                          T(0.657837705860311),
                                                                          T(0.219550214526392),
                                                                          T(-0.993510363631668),
                                                                          T(0.300046887453750),
                                                                          T(-0.066677086100833),
                                                                          T(0.007143973510804))),
                                DerivativeCoefficientRow{T, 1, 10}(SVector(T(0.000391462486334),
                                                                           T(-0.003009678529462),
                                                                           T(0.017964428161274),
                                                                           T(-0.111241495585250),
                                                                           T(0.607178485175455),
                                                                           T(0.248233272383674),
                                                                           T(-0.999990342064734),
                                                                           T(0.299997102619420),
                                                                           T(-0.066666022804316),
                                                                           T(0.007142788157605))))
        upper_coef_minus = .-lower_coef_plus
        central_coef_minus = .-central_coef_plus
        lower_coef_minus = .-upper_coef_plus

        left_boundary_central = (left_boundary_plus .+ left_boundary_minus) ./ 2
        right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
        upper_coef_central = widening_plus(upper_coef_plus, upper_coef_minus) / 2
        central_coef_central = (central_coef_plus + central_coef_minus) / 2
        lower_coef_central = widening_plus(lower_coef_plus, lower_coef_minus) / 2

        if source.kind === :plus
            left_boundary = left_boundary_plus
            right_boundary = right_boundary_plus
            upper_coef = upper_coef_plus
            central_coef = central_coef_plus
            lower_coef = lower_coef_plus
        elseif source.kind === :minus
            left_boundary = left_boundary_minus
            right_boundary = right_boundary_minus
            upper_coef = upper_coef_minus
            central_coef = central_coef_minus
            lower_coef = lower_coef_minus
        elseif source.kind === :central
            left_boundary = left_boundary_central
            right_boundary = right_boundary_central
            upper_coef = upper_coef_central
            central_coef = central_coef_central
            lower_coef = lower_coef_central
        end
        DerivativeCoefficients(left_boundary, right_boundary,
                               left_boundary_derivatives, right_boundary_derivatives,
                               lower_coef, central_coef, upper_coef,
                               left_weights, right_weights, mode, 1, order, source)

    elseif order == 8
        left_boundary_plus = (DerivativeCoefficientRow{T, 1, 8}(SVector(T(-4.543535228767764),
                                                                        T(7.407530637111112),
                                                                        T(-6.039882962893306),
                                                                        T(6.319056442849687),
                                                                        T(-5.218600967009760),
                                                                        T(2.862774280523586),
                                                                        T(-0.916614110722202),
                                                                        T(0.129271908908647))),
                              DerivativeCoefficientRow{T, 1, 8}(SVector(T(-0.453549041115923),
                                                                        T(-1.302492113991826),
                                                                        T(3.489999172042541),
                                                                        T(-3.444330809240609),
                                                                        T(2.838607549200516),
                                                                        T(-1.552270889875925),
                                                                        T(0.492496100076883),
                                                                        T(-0.068459967095658))),
                              DerivativeCoefficientRow{T, 1, 8}(SVector(T(-0.363628825639525),
                                                                        T(0.853654212904449),
                                                                        T(-3.088407204888069),
                                                                        T(5.015425547907850),
                                                                        T(-4.004801751671350),
                                                                        T(2.180026955384745),
                                                                        T(-0.686633363269485),
                                                                        T(0.094364429271383))),
                              DerivativeCoefficientRow{T, 1, 9}(SVector(T(0.762082954127772),
                                                                        T(-2.170557710810797),
                                                                        T(3.478053889004435),
                                                                        T(-5.867335588289929),
                                                                        T(6.192970273593883),
                                                                        T(-3.302588795031592),
                                                                        T(1.063639172959331),
                                                                        T(-0.159848933758619),
                                                                        T(0.003584738205516))),
                              DerivativeCoefficientRow{T, 1, 10}(SVector(T(-0.779310457381702),
                                                                         T(2.239677793384923),
                                                                         T(-3.827100500011441),
                                                                         T(4.712159306750459),
                                                                         T(-5.137055265951099),
                                                                         T(3.872716058460845),
                                                                         T(-1.328788707540612),
                                                                         T(0.279769354952040),
                                                                         T(-0.035630647403792),
                                                                         T(0.003563064740379))),
                              DerivativeCoefficientRow{T, 1, 11}(SVector(T(0.465743750422248),
                                                                         T(-1.344731351052584),
                                                                         T(2.331043287406878),
                                                                         T(-3.047561141011216),
                                                                         T(2.292200936606628),
                                                                         T(-2.000921871702452),
                                                                         T(1.735230507732654),
                                                                         T(-0.565633740229944),
                                                                         T(0.166797761556549),
                                                                         T(-0.035742377476403),
                                                                         T(0.003574237747640))),
                              DerivativeCoefficientRow{T, 1, 12}(SVector(T(-0.152129575964355),
                                                                         T(0.442423175108889),
                                                                         T(-0.777066846787855),
                                                                         T(1.037949230819469),
                                                                         T(-0.863784398125473),
                                                                         T(0.016456074826725),
                                                                         T(-0.609333986755292),
                                                                         T(1.270912394669787),
                                                                         T(-0.499931428249499),
                                                                         T(0.166643809416500),
                                                                         T(-0.035709387732107),
                                                                         T(0.003570938773211))),
                              DerivativeCoefficientRow{T, 1, 13}(SVector(T(0.021191153042340),
                                                                         T(-0.062255226314800),
                                                                         T(0.110992139857950),
                                                                         T(-0.150604165918153),
                                                                         T(0.129397484811924),
                                                                         T(-0.003021916733589),
                                                                         T(-0.477285585226883),
                                                                         T(-0.452945127565490),
                                                                         T(1.250010506391703),
                                                                         T(-0.500004202556681),
                                                                         T(0.166668067518894),
                                                                         T(-0.035714585896906),
                                                                         T(0.003571458589691))))
        right_boundary_plus = (DerivativeCoefficientRow{T, 1, 8}(SVector(T(2.859874032009999),
                                                                         T(-2.496614758534824),
                                                                         T(-2.677060761960096),
                                                                         T(5.621063994683933),
                                                                         T(-5.783097535404512),
                                                                         T(3.445381563525398),
                                                                         T(-1.126431994360043),
                                                                         T(0.156885460040145))),
                               DerivativeCoefficientRow{T, 1, 8}(SVector(T(1.345693566063129),
                                                                         T(-1.302492113991826),
                                                                         T(1.141706918569862),
                                                                         T(-2.908441540909868),
                                                                         T(3.019313971148636),
                                                                         T(-1.807168047258096),
                                                                         T(0.595116388874535),
                                                                         T(-0.083729142496372))),
                               DerivativeCoefficientRow{T, 1, 8}(SVector(T(-0.820405565688017),
                                                                         T(2.609472227757943),
                                                                         T(-3.088407204888069),
                                                                         T(3.484597843136334),
                                                                         T(-3.857624531277929),
                                                                         T(2.342290206887081),
                                                                         T(-0.781537421615602),
                                                                         T(0.111614445688256))),
                               DerivativeCoefficientRow{T, 1, 9}(SVector(T(0.856714174722325),
                                                                         T(-2.570489621820486),
                                                                         T(5.006006752335181),
                                                                         T(-5.867335588289929),
                                                                         T(4.740822502032001),
                                                                         T(-3.056514319183563),
                                                                         T(1.041960251746113),
                                                                         T(-0.151164151541642),
                                                                         T(0 // 1))),
                               DerivativeCoefficientRow{T, 1, 10}(SVector(T(-0.703240829260628),
                                                                          T(2.105632720819886),
                                                                          T(-3.973113158628417),
                                                                          T(6.155527337004563),
                                                                          T(-5.137055265951099),
                                                                          T(2.285035554917741),
                                                                          T(-0.861879726233252),
                                                                          T(0.129093367331207),
                                                                          T(0 // 1),
                                                                          T(0 // 1))),
                               DerivativeCoefficientRow{T, 1, 11}(SVector(T(0.386987393250901),
                                                                          T(-1.155059892802726),
                                                                          T(2.169559171520989),
                                                                          T(-3.292914812571790),
                                                                          T(3.884860065879826),
                                                                          T(-2.000921871702452),
                                                                          T(0.016471277599305),
                                                                          T(-0.003024268261320),
                                                                          T(-0.005957062912734),
                                                                          T(0 // 1),
                                                                          T(0 // 1))),
                               DerivativeCoefficientRow{T, 1, 12}(SVector(T(-0.123792751524547),
                                                                          T(0.366132898367712),
                                                                          T(-0.682705661607576),
                                                                          T(1.059544699130806),
                                                                          T(-1.331725203695376),
                                                                          T(1.733628912797258),
                                                                          T(-0.609333986755292),
                                                                          T(-0.477216117555181),
                                                                          T(0.071418775464214),
                                                                          T(-0.005951564622018),
                                                                          T(0 // 1),
                                                                          T(0 // 1))),
                               DerivativeCoefficientRow{T, 1, 13}(SVector(T(0.017461279108068),
                                                                          T(-0.050902118640814),
                                                                          T(0.093838300828537),
                                                                          T(-0.159256775472928),
                                                                          T(0.280428434137662),
                                                                          T(-0.565193930229378),
                                                                          T(1.271097399580030),
                                                                          T(-0.452945127565490),
                                                                          T(-0.500004202556681),
                                                                          T(0.071429171793812),
                                                                          T(-0.005952430982818),
                                                                          T(0 // 1),
                                                                          T(0 // 1))))
        upper_coef_plus = SVector(T(5 // 4),
                                  T(-1 // 2),
                                  T(1 // 6),
                                  T(-1 // 28),
                                  T(1 // 280))
        central_coef_plus = T(-9 // 20)
        lower_coef_plus = SVector(T(-1 // 2),
                                  T(1 // 14),
                                  T(-1 // 168))
        left_weights = SVector(T(0.135072905573093),
                               T(0.743524909020311),
                               T(0.994416146403454),
                               T(0.996287139164853),
                               T(1.002347369935375),
                               T(0.999214048865769),
                               T(1.000137162311921),
                               T(0.999991594957283))
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        left_boundary_minus = (DerivativeCoefficientRow{T, 1, 8}(SVector(T(-2.859874032009999),
                                                                         T(2.496614758534824),
                                                                         T(2.677060761960096),
                                                                         T(-5.621063994683933),
                                                                         T(5.783097535404512),
                                                                         T(-3.445381563525398),
                                                                         T(1.126431994360043),
                                                                         T(-0.156885460040145))),
                               DerivativeCoefficientRow{T, 1, 8}(SVector(T(-1.345693566063129),
                                                                         T(1.302492113991826),
                                                                         T(-1.141706918569862),
                                                                         T(2.908441540909868),
                                                                         T(-3.019313971148636),
                                                                         T(1.807168047258096),
                                                                         T(-0.595116388874535),
                                                                         T(0.083729142496372))),
                               DerivativeCoefficientRow{T, 1, 8}(SVector(T(0.820405565688017),
                                                                         T(-2.609472227757943),
                                                                         T(3.088407204888069),
                                                                         T(-3.484597843136334),
                                                                         T(3.857624531277929),
                                                                         T(-2.342290206887081),
                                                                         T(0.781537421615602),
                                                                         T(-0.111614445688256))),
                               DerivativeCoefficientRow{T, 1, 9}(SVector(T(-0.856714174722325),
                                                                         T(2.570489621820486),
                                                                         T(-5.006006752335181),
                                                                         T(5.867335588289929),
                                                                         T(-4.740822502032001),
                                                                         T(3.056514319183563),
                                                                         T(-1.041960251746113),
                                                                         T(0.151164151541642),
                                                                         T(0 // 1))),
                               DerivativeCoefficientRow{T, 1, 10}(SVector(T(0.703240829260628),
                                                                          T(-2.105632720819886),
                                                                          T(3.973113158628417),
                                                                          T(-6.155527337004563),
                                                                          T(5.137055265951099),
                                                                          T(-2.285035554917741),
                                                                          T(0.861879726233252),
                                                                          T(-0.129093367331207),
                                                                          T(0 // 1),
                                                                          T(0 // 1))),
                               DerivativeCoefficientRow{T, 1, 11}(SVector(T(-0.386987393250901),
                                                                          T(1.155059892802726),
                                                                          T(-2.169559171520989),
                                                                          T(3.292914812571790),
                                                                          T(-3.884860065879826),
                                                                          T(2.000921871702452),
                                                                          T(-0.016471277599305),
                                                                          T(0.003024268261320),
                                                                          T(0.005957062912734),
                                                                          T(0 // 1),
                                                                          T(0 // 1))),
                               DerivativeCoefficientRow{T, 1, 12}(SVector(T(0.123792751524547),
                                                                          T(-0.366132898367712),
                                                                          T(0.682705661607576),
                                                                          T(-1.059544699130806),
                                                                          T(1.331725203695376),
                                                                          T(-1.733628912797258),
                                                                          T(0.609333986755292),
                                                                          T(0.477216117555181),
                                                                          T(-0.071418775464214),
                                                                          T(0.005951564622018),
                                                                          T(0 // 1),
                                                                          T(0 // 1))),
                               DerivativeCoefficientRow{T, 1, 13}(SVector(T(-0.017461279108068),
                                                                          T(0.050902118640814),
                                                                          T(-0.093838300828537),
                                                                          T(0.159256775472928),
                                                                          T(-0.280428434137662),
                                                                          T(0.565193930229378),
                                                                          T(-1.271097399580030),
                                                                          T(0.452945127565490),
                                                                          T(0.500004202556681),
                                                                          T(-0.071429171793812),
                                                                          T(0.005952430982818),
                                                                          T(0 // 1),
                                                                          T(0 // 1))))
        right_boundary_minus = (DerivativeCoefficientRow{T, 1, 8}(SVector(T(4.543535228767764),
                                                                          T(-7.407530637111112),
                                                                          T(6.039882962893306),
                                                                          T(-6.319056442849687),
                                                                          T(5.218600967009760),
                                                                          T(-2.862774280523586),
                                                                          T(0.916614110722202),
                                                                          T(-0.129271908908647))),
                                DerivativeCoefficientRow{T, 1, 8}(SVector(T(0.453549041115923),
                                                                          T(1.302492113991826),
                                                                          T(-3.489999172042541),
                                                                          T(3.444330809240609),
                                                                          T(-2.838607549200516),
                                                                          T(1.552270889875925),
                                                                          T(-0.492496100076883),
                                                                          T(0.068459967095658))),
                                DerivativeCoefficientRow{T, 1, 8}(SVector(T(0.363628825639525),
                                                                          T(-0.853654212904449),
                                                                          T(3.088407204888069),
                                                                          T(-5.015425547907850),
                                                                          T(4.004801751671350),
                                                                          T(-2.180026955384745),
                                                                          T(0.686633363269485),
                                                                          T(-0.094364429271383))),
                                DerivativeCoefficientRow{T, 1, 9}(SVector(T(-0.762082954127772),
                                                                          T(2.170557710810797),
                                                                          T(-3.478053889004435),
                                                                          T(5.867335588289929),
                                                                          T(-6.192970273593883),
                                                                          T(3.302588795031592),
                                                                          T(-1.063639172959331),
                                                                          T(0.159848933758619),
                                                                          T(-0.003584738205516))),
                                DerivativeCoefficientRow{T, 1, 10}(SVector(T(0.779310457381702),
                                                                           T(-2.239677793384923),
                                                                           T(3.827100500011441),
                                                                           T(-4.712159306750459),
                                                                           T(5.137055265951099),
                                                                           T(-3.872716058460845),
                                                                           T(1.328788707540612),
                                                                           T(-0.279769354952040),
                                                                           T(0.035630647403792),
                                                                           T(-0.003563064740379))),
                                DerivativeCoefficientRow{T, 1, 11}(SVector(T(-0.465743750422248),
                                                                           T(1.344731351052584),
                                                                           T(-2.331043287406878),
                                                                           T(3.047561141011216),
                                                                           T(-2.292200936606628),
                                                                           T(2.000921871702452),
                                                                           T(-1.735230507732654),
                                                                           T(0.565633740229944),
                                                                           T(-0.166797761556549),
                                                                           T(0.035742377476403),
                                                                           T(-0.003574237747640))),
                                DerivativeCoefficientRow{T, 1, 12}(SVector(T(0.152129575964355),
                                                                           T(-0.442423175108889),
                                                                           T(0.777066846787855),
                                                                           T(-1.037949230819469),
                                                                           T(0.863784398125473),
                                                                           T(-0.016456074826725),
                                                                           T(0.609333986755292),
                                                                           T(-1.270912394669787),
                                                                           T(0.499931428249499),
                                                                           T(-0.166643809416500),
                                                                           T(0.035709387732107),
                                                                           T(-0.003570938773211))),
                                DerivativeCoefficientRow{T, 1, 13}(SVector(T(-0.021191153042340),
                                                                           T(0.062255226314800),
                                                                           T(-0.110992139857950),
                                                                           T(0.150604165918153),
                                                                           T(-0.129397484811924),
                                                                           T(0.003021916733589),
                                                                           T(0.477285585226883),
                                                                           T(0.452945127565490),
                                                                           T(-1.250010506391703),
                                                                           T(0.500004202556681),
                                                                           T(-0.166668067518894),
                                                                           T(0.035714585896906),
                                                                           T(-0.003571458589691))))
        upper_coef_minus = .-lower_coef_plus
        central_coef_minus = .-central_coef_plus
        lower_coef_minus = .-upper_coef_plus

        left_boundary_central = (left_boundary_plus .+ left_boundary_minus) ./ 2
        right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
        upper_coef_central = widening_plus(upper_coef_plus, upper_coef_minus) / 2
        central_coef_central = (central_coef_plus + central_coef_minus) / 2
        lower_coef_central = widening_plus(lower_coef_plus, lower_coef_minus) / 2

        if source.kind === :plus
            left_boundary = left_boundary_plus
            right_boundary = right_boundary_plus
            upper_coef = upper_coef_plus
            central_coef = central_coef_plus
            lower_coef = lower_coef_plus
        elseif source.kind === :minus
            left_boundary = left_boundary_minus
            right_boundary = right_boundary_minus
            upper_coef = upper_coef_minus
            central_coef = central_coef_minus
            lower_coef = lower_coef_minus
        elseif source.kind === :central
            left_boundary = left_boundary_central
            right_boundary = right_boundary_central
            upper_coef = upper_coef_central
            central_coef = central_coef_central
            lower_coef = lower_coef_central
        end
        DerivativeCoefficients(left_boundary, right_boundary,
                               left_boundary_derivatives, right_boundary_derivatives,
                               lower_coef, central_coef, upper_coef,
                               left_weights, right_weights, mode, 1, order, source)
    elseif order == 9
        left_boundary_plus = (DerivativeCoefficientRow{T, 1, 8}(SVector(T(-3.847347650663411),
                                                                        T(5.305832094510697),
                                                                        T(-2.206406993775266),
                                                                        T(1.013492019588300),
                                                                        T(-0.290876732169925),
                                                                        T(0.012990985098595),
                                                                        T(0.014691680302891),
                                                                        T(-0.002375402891880))),
                              DerivativeCoefficientRow{T, 1, 8}(SVector(T(-0.849765846168166),
                                                                        T(-0.164536922777836),
                                                                        T(1.458999779350386),
                                                                        T(-0.624172228642107),
                                                                        T(0.217690528814413),
                                                                        T(-0.038400185835984),
                                                                        T(-0.000751894687916),
                                                                        T(0.000936769947210))),
                              DerivativeCoefficientRow{T, 1, 8}(SVector(T(0.156101537218983),
                                                                        T(-0.673899820896248),
                                                                        T(-0.353171152479185),
                                                                        T(1.237055655554625),
                                                                        T(-0.493549096557796),
                                                                        T(0.152030190829188),
                                                                        T(-0.026013037157170),
                                                                        T(0.001445723487603))),
                              DerivativeCoefficientRow{T, 1, 9}(SVector(T(0.042022234353631),
                                                                        T(-0.061708708869846),
                                                                        T(-0.293462615788667),
                                                                        T(-0.655092084194954),
                                                                        T(1.338460288188073),
                                                                        T(-0.477637843725215),
                                                                        T(0.125692979905497),
                                                                        T(-0.019866614701244),
                                                                        T(0.001592364832725))),
                              DerivativeCoefficientRow{T, 1, 10}(SVector(T(-0.111400431425707),
                                                                         T(0.288830641698750),
                                                                         T(-0.349486849379579),
                                                                         T(-0.090855368055662),
                                                                         T(-0.637622496465684),
                                                                         T(1.197428491898929),
                                                                         T(-0.379929338644352),
                                                                         T(0.099262149995093),
                                                                         T(-0.017809902023915),
                                                                         T(0.001583102402126))),
                              DerivativeCoefficientRow{T, 1, 11}(SVector(T(0.077432608323769),
                                                                         T(-0.212549757428126),
                                                                         T(0.318942068005413),
                                                                         T(-0.270520964395790),
                                                                         T(-0.337417719864865),
                                                                         T(-0.358427819360282),
                                                                         T(1.041349345910250),
                                                                         T(-0.337866369883873),
                                                                         T(0.095347065761010),
                                                                         T(-0.017877574830189),
                                                                         T(0.001589117762683))),
                              DerivativeCoefficientRow{T, 1, 12}(SVector(T(-0.025236556853842),
                                                                         T(0.072726393490041),
                                                                         T(-0.122002325006802),
                                                                         T(0.135175121239691),
                                                                         T(0.009603432702573),
                                                                         T(-0.600133750769375),
                                                                         T(-0.217473389128111),
                                                                         T(1.001634179423633),
                                                                         T(-0.333239014480747),
                                                                         T(0.095211146994499),
                                                                         T(-0.017852090061469),
                                                                         T(0.001586852449908))),
                              DerivativeCoefficientRow{T, 1, 13}(SVector(T(0.003253562097921),
                                                                         T(-0.009982358570495),
                                                                         T(0.018651197497656),
                                                                         T(-0.023637956771815),
                                                                         T(-0.001425199511631),
                                                                         T(0.131480108475396),
                                                                         T(-0.663662797821441),
                                                                         T(-0.200334853975450),
                                                                         T(1.000031353071446),
                                                                         T(-0.333343784357149),
                                                                         T(0.095241081244900),
                                                                         T(-0.017857702733419),
                                                                         T(0.001587351354082))))
        right_boundary_plus = (DerivativeCoefficientRow{T, 1, 8}(SVector(T(3.620708549227146),
                                                                         T(-4.668297867288131),
                                                                         T(1.155528368173997),
                                                                         T(0.312826539744191),
                                                                         T(-0.834151419031349),
                                                                         T(0.577610174596973),
                                                                         T(-0.188521368233643),
                                                                         T(0.024297022810815))),
                               DerivativeCoefficientRow{T, 1, 8}(SVector(T(0.965815598660001),
                                                                         T(-0.164536922777836),
                                                                         T(-0.908049421078713),
                                                                         T(-0.083620281085714),
                                                                         T(0.393678772789026),
                                                                         T(-0.288610596572745),
                                                                         T(0.098892474879078),
                                                                         T(-0.013569624813097))),
                               DerivativeCoefficientRow{T, 1, 8}(SVector(T(-0.298065831134292),
                                                                         T(1.082782134064773),
                                                                         T(-0.353171152479185),
                                                                         T(-0.295123506766509),
                                                                         T(-0.353521170706653),
                                                                         T(0.321402551225085),
                                                                         T(-0.123119021850349),
                                                                         T(0.018815997647129))),
                               DerivativeCoefficientRow{T, 1, 8}(SVector(T(0.136143177613706),
                                                                         T(-0.460616274447069),
                                                                         T(1.230093775086648),
                                                                         T(-0.655092084194954),
                                                                         T(-0.091386945507643),
                                                                         T(-0.271073724260165),
                                                                         T(0.135644690427143),
                                                                         T(-0.023712614717665))),
                               DerivativeCoefficientRow{T, 1, 9}(SVector(T(-0.038846416509197),
                                                                         T(0.159713196329492),
                                                                         T(-0.487916801206936),
                                                                         T(1.330674763618661),
                                                                         T(-0.637622496465684),
                                                                         T(-0.336140477050505),
                                                                         T(0.009580737882073),
                                                                         T(-0.001421384600561),
                                                                         T(0.001978878002657))),
                               DerivativeCoefficientRow{T, 1, 10}(SVector(T(0.001741530715904),
                                                                          T(-0.028280147304211),
                                                                          T(0.150866330330903),
                                                                          T(-0.476663868728418),
                                                                          T(1.201978396005712),
                                                                          T(-0.358427819360282),
                                                                          T(-0.600990472295380),
                                                                          T(0.131626419872670),
                                                                          T(-0.023836766440252),
                                                                          T(0.001986397203354))),
                               DerivativeCoefficientRow{T, 1, 11}(SVector(T(0.001966712997663),
                                                                          T(-0.000552949948955),
                                                                          T(-0.025777097364538),
                                                                          T(0.125257861138543),
                                                                          T(-0.380829314016763),
                                                                          T(1.039864885770100),
                                                                          T(-0.217473389128111),
                                                                          T(-0.663454208753449),
                                                                          T(0.142816720491749),
                                                                          T(-0.023802786748625),
                                                                          T(0.001983565562385))),
                               DerivativeCoefficientRow{T, 1, 12}(SVector(T(-0.000318085095301),
                                                                          T(0.000689125428294),
                                                                          T(0.001433061100451),
                                                                          T(-0.019804065688314),
                                                                          T(0.099528563655892),
                                                                          T(-0.337490809257689),
                                                                          T(1.001949091797688),
                                                                          T(-0.200334853975450),
                                                                          T(-0.666687568714297),
                                                                          T(0.142861621867349),
                                                                          T(-0.023810270311225),
                                                                          T(0.001984189192602))))
        upper_coef_plus = SVector(T(1 // 1),
                                  T(-1 // 3),
                                  T(2 // 21),
                                  T(-1 // 56),
                                  T(1 // 630))
        central_coef_plus = T(-1 // 5)
        lower_coef_plus = SVector(T(-2 // 3),
                                  T(1 // 7),
                                  T(-1 // 42),
                                  T(1 // 504))
        left_weights = SVector(T(0.133903652199973),
                               T(0.735616919420790),
                               T(0.991210410068642),
                               T(0.996820298138376),
                               T(1.002652503824264),
                               T(0.998857117185048),
                               T(1.000283036644833),
                               T(0.999968647911538))
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        # upper left block of the `D^-` derivative matrix
        # Equivalent to the nonzero entries in the first eight rows of `h * Dm` from Ken's code
        left_boundary_minus = (DerivativeCoefficientRow{T, 1, 8}(SVector(T(-3.620708549227146),
                                                                         T(4.668297867288131),
                                                                         T(-1.155528368173997),
                                                                         T(-0.312826539744191),
                                                                         T(0.834151419031349),
                                                                         T(-0.577610174596973),
                                                                         T(0.188521368233643),
                                                                         T(-0.024297022810815))),
                               DerivativeCoefficientRow{T, 1, 8}(SVector(T(-0.965815598660001),
                                                                         T(0.164536922777836),
                                                                         T(0.908049421078713),
                                                                         T(0.083620281085714),
                                                                         T(-0.393678772789026),
                                                                         T(0.288610596572745),
                                                                         T(-0.098892474879078),
                                                                         T(0.013569624813097))),
                               DerivativeCoefficientRow{T, 1, 8}(SVector(T(0.298065831134292),
                                                                         T(-1.082782134064773),
                                                                         T(0.353171152479185),
                                                                         T(0.295123506766509),
                                                                         T(0.353521170706653),
                                                                         T(-0.321402551225085),
                                                                         T(0.123119021850349),
                                                                         T(-0.018815997647129))),
                               DerivativeCoefficientRow{T, 1, 9}(SVector(T(-0.136143177613706),
                                                                         T(0.460616274447069),
                                                                         T(-1.230093775086648),
                                                                         T(0.655092084194954),
                                                                         T(0.091386945507643),
                                                                         T(0.271073724260165),
                                                                         T(-0.135644690427143),
                                                                         T(0.023712614717665),
                                                                         T(0 // 1))),
                               DerivativeCoefficientRow{T, 1, 10}(SVector(T(0.038846416509197),
                                                                          T(-0.159713196329492),
                                                                          T(0.487916801206936),
                                                                          T(-1.330674763618661),
                                                                          T(0.637622496465684),
                                                                          T(0.336140477050505),
                                                                          T(-0.009580737882073),
                                                                          T(0.001421384600561),
                                                                          T(-0.001978878002657),
                                                                          T(0 // 1))),
                               DerivativeCoefficientRow{T, 1, 11}(SVector(T(-0.001741530715904),
                                                                          T(0.028280147304211),
                                                                          T(-0.150866330330903),
                                                                          T(0.476663868728418),
                                                                          T(-1.201978396005712),
                                                                          T(0.358427819360282),
                                                                          T(0.600990472295380),
                                                                          T(-0.131626419872670),
                                                                          T(0.023836766440252),
                                                                          T(-0.001986397203354),
                                                                          T(0 // 1))),
                               DerivativeCoefficientRow{T, 1, 12}(SVector(T(-0.001966712997663),
                                                                          T(0.000552949948955),
                                                                          T(0.025777097364538),
                                                                          T(-0.125257861138543),
                                                                          T(0.380829314016763),
                                                                          T(-1.039864885770100),
                                                                          T(0.217473389128111),
                                                                          T(0.663454208753449),
                                                                          T(-0.142816720491749),
                                                                          T(0.023802786748625),
                                                                          T(-0.001983565562385),
                                                                          T(0 // 1))),
                               DerivativeCoefficientRow{T, 1, 13}(SVector(T(0.000318085095301),
                                                                          T(-0.000689125428294),
                                                                          T(-0.001433061100451),
                                                                          T(0.019804065688314),
                                                                          T(-0.099528563655892),
                                                                          T(0.337490809257689),
                                                                          T(-1.001949091797688),
                                                                          T(0.200334853975450),
                                                                          T(0.666687568714297),
                                                                          T(-0.142861621867349),
                                                                          T(0.023810270311225),
                                                                          T(-0.001984189192602),
                                                                          T(0 // 1))))
        right_boundary_minus = (DerivativeCoefficientRow{T, 1, 8}(SVector(T(3.847347650663411),
                                                                          T(-5.305832094510697),
                                                                          T(2.206406993775266),
                                                                          T(-1.013492019588300),
                                                                          T(0.290876732169925),
                                                                          T(-0.012990985098595),
                                                                          T(-0.014691680302891),
                                                                          T(0.002375402891880))),
                                DerivativeCoefficientRow{T, 1, 8}(SVector(T(0.849765846168166),
                                                                          T(0.164536922777836),
                                                                          T(-1.458999779350386),
                                                                          T(0.624172228642107),
                                                                          T(-0.217690528814413),
                                                                          T(0.038400185835984),
                                                                          T(0.000751894687916),
                                                                          T(-0.000936769947210))),
                                DerivativeCoefficientRow{T, 1, 8}(SVector(T(-0.156101537218983),
                                                                          T(0.673899820896248),
                                                                          T(0.353171152479185),
                                                                          T(-1.237055655554625),
                                                                          T(0.493549096557796),
                                                                          T(-0.152030190829188),
                                                                          T(0.026013037157170),
                                                                          T(-0.001445723487603))),
                                DerivativeCoefficientRow{T, 1, 9}(SVector(T(-0.042022234353631),
                                                                          T(0.061708708869846),
                                                                          T(0.293462615788667),
                                                                          T(0.655092084194954),
                                                                          T(-1.338460288188073),
                                                                          T(0.477637843725215),
                                                                          T(-0.125692979905497),
                                                                          T(0.019866614701244),
                                                                          T(-0.001592364832725))),
                                DerivativeCoefficientRow{T, 1, 10}(SVector(T(0.111400431425707),
                                                                           T(-0.288830641698750),
                                                                           T(0.349486849379579),
                                                                           T(0.090855368055662),
                                                                           T(0.637622496465684),
                                                                           T(-1.197428491898929),
                                                                           T(0.379929338644352),
                                                                           T(-0.099262149995093),
                                                                           T(0.017809902023915),
                                                                           T(-0.001583102402126))),
                                DerivativeCoefficientRow{T, 1, 11}(SVector(T(-0.077432608323769),
                                                                           T(0.212549757428126),
                                                                           T(-0.318942068005413),
                                                                           T(0.270520964395790),
                                                                           T(0.337417719864865),
                                                                           T(0.358427819360282),
                                                                           T(-1.041349345910250),
                                                                           T(0.337866369883873),
                                                                           T(-0.095347065761010),
                                                                           T(0.017877574830189),
                                                                           T(-0.001589117762683))),
                                DerivativeCoefficientRow{T, 1, 12}(SVector(T(0.025236556853842),
                                                                           T(-0.072726393490041),
                                                                           T(0.122002325006802),
                                                                           T(-0.135175121239691),
                                                                           T(-0.009603432702573),
                                                                           T(0.600133750769375),
                                                                           T(0.217473389128111),
                                                                           T(-1.001634179423633),
                                                                           T(0.333239014480747),
                                                                           T(-0.095211146994499),
                                                                           T(0.017852090061469),
                                                                           T(-0.001586852449908))),
                                DerivativeCoefficientRow{T, 1, 13}(SVector(T(-0.003253562097921),
                                                                           T(0.009982358570495),
                                                                           T(-0.018651197497656),
                                                                           T(0.023637956771815),
                                                                           T(0.001425199511631),
                                                                           T(-0.131480108475396),
                                                                           T(0.663662797821441),
                                                                           T(0.200334853975450),
                                                                           T(-1.000031353071446),
                                                                           T(0.333343784357149),
                                                                           T(-0.095241081244900),
                                                                           T(0.017857702733419),
                                                                           T(-0.001587351354082))))
        upper_coef_minus = .-lower_coef_plus
        central_coef_minus = .-central_coef_plus
        lower_coef_minus = .-upper_coef_plus

        left_boundary_central = (left_boundary_plus .+ left_boundary_minus) ./ 2
        right_boundary_central = (right_boundary_plus .+ right_boundary_minus) ./ 2
        upper_coef_central = widening_plus(upper_coef_plus, upper_coef_minus) / 2
        central_coef_central = (central_coef_plus + central_coef_minus) / 2
        lower_coef_central = widening_plus(lower_coef_plus, lower_coef_minus) / 2

        if source.kind === :plus
            left_boundary = left_boundary_plus
            right_boundary = right_boundary_plus
            upper_coef = upper_coef_plus
            central_coef = central_coef_plus
            lower_coef = lower_coef_plus
        elseif source.kind === :minus
            left_boundary = left_boundary_minus
            right_boundary = right_boundary_minus
            upper_coef = upper_coef_minus
            central_coef = central_coef_minus
            lower_coef = lower_coef_minus
        elseif source.kind === :central
            left_boundary = left_boundary_central
            right_boundary = right_boundary_central
            upper_coef = upper_coef_central
            central_coef = central_coef_central
            lower_coef = lower_coef_central
        end
        DerivativeCoefficients(left_boundary, right_boundary,
                               left_boundary_derivatives, right_boundary_derivatives,
                               lower_coef, central_coef, upper_coef,
                               left_weights, right_weights, mode, 1, order, source)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end
