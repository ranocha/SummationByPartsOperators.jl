
"""
    DienerDorbandSchnetterTiglio2007()

Coefficients of the SBP operators given in
- Diener, Dorband, Schnetter, Tiglio (2007)
  Optimized high-order derivative and dissipation operators satisfying
  summation by parts, and applications in three-dimensional multi-block
  evolutions.
  Journal of Scientific Computing 32.1, pp. 109-145.

See also (second- and fourth-order operators)
- Mattsson, Nordström (2004)
  Summation by parts operators for finite difference approximations of second
  derivatives.
  Journal of Computational Physics 199, pp. 503-540.

The dissipation operators proposed by Diener, Dorband, Schnetter, Tiglio (2007)
for the diagonal-norm operators are the same as the ones of
- Mattsson, Svärd, Nordström (2004)
  Stable and Accurate Artificial Dissipation.
  Journal of Scientific Computing 21.1, pp. 57-79.
"""
struct DienerDorbandSchnetterTiglio2007 <: SourceOfCoefficients end

function Base.show(io::IO, source::DienerDorbandSchnetterTiglio2007)
    if get(io, :compact, false)
      summary(io, source)
    else
        print(io,
            "Diener, Dorband, Schnetter, Tiglio (2007) \n",
            "  Optimized high-order derivative and dissipation operators satisfying \n",
            "  summation by parts, and applications in three-dimensional multi-block \n",
            "  evolutions. \n",
            "  Journal of Scientific Computing 32.1, pp. 109-145. \n",
            "See also (second- and fourth-order operators) \n",
            "  Mattsson, Nordström (2004) \n",
            "  Summation by parts operators for finite difference approximations of second \n",
            "    derivatives. \n",
            "  Journal of Computational Physics 199, pp. 503-540.\n",
            "The dissipation operators proposed by Diener, Dorband, Schnetter, Tiglio (2007)\n",
            "for the diagonal-norm operators are the same as the ones of\n",
            "  Mattsson, Svärd, Nordström (2004)\n",
            "  Stable and Accurate Artificial Dissipation.\n",
            "  Journal of Scientific Computing 21.1, pp. 57-79.")
    end
end


function first_derivative_coefficients(source::DienerDorbandSchnetterTiglio2007, order::Int,
                                       T=Float64, mode=FastMode())
    if order == 2
        left_boundary = (
            DerivativeCoefficientRow{T,1,2}(SVector(-one(T), one(T))),
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
                                left_weights, right_weights, mode, 1, order, source)
    elseif order == 4
        left_boundary = (
            # q1
            DerivativeCoefficientRow{T,1,4}(SVector(T(-24//17),
                                                    T(59//34),
                                                    T(-4//17),
                                                    T(-3//34) )),
            # q2
            DerivativeCoefficientRow{T,1,3}(SVector(T(-1//2),
                                                    T(0),
                                                    T(1//2))),
            # q3
            DerivativeCoefficientRow{T,1,5}(SVector(T(4//43),
                                                    T(-59//86),
                                                    T(0),
                                                    T(59//86),
                                                    T(-4//43))),
            # q4
            DerivativeCoefficientRow{T,1,6}(SVector(T(3//98),
                                                    T(0),
                                                    T(-59//98),
                                                    T(0),
                                                    T(32//49),
                                                    T(-4//49))),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(2//3), T(-1//12))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(17//48),
                                T(59//48),
                                T(43//48),
                                T(49//48) )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, mode, 1, order, source)
    elseif order == 6
        left_boundary = (
            # q1
            DerivativeCoefficientRow{T,1,6}(SVector(convert(T, big"-1.582533518939116418785258993332844897062"),
                                                    convert(T, big"2.033426786468126253898161347360808173712"),
                                                    convert(T, big"-0.1417052898146741610733887894481170575600"),
                                                    convert(T, big"-0.4501096599735708523162117824920488989702"),
                                                    convert(T, big"0.1042956382142412661862395105494407610836"),
                                                    convert(T, big"0.03662604404499391209045870736276191879693") )),
            # q2
            DerivativeCoefficientRow{T,1,6}(SVector(convert(T, big"-0.4620701275035953590186631853846278325646"),
                                                    convert(T, 0),
                                                    convert(T, big"0.2873679417026202568532985205129449923126"),
                                                    convert(T, big"0.2585974499280928196267362923074433487080"),
                                                    convert(T, big"-0.06894808744606961472005221923058251153103"),
                                                    convert(T, big"-0.01494717668104810274131940820517799692506") )),
            # q3
            DerivativeCoefficientRow{T,1,6}(SVector(convert(T, big"0.07134398748360337973038301686379010397038"),
                                                    convert(T, big"-0.6366933020423417826592908754928085932593"),
                                                    convert(T, 0),
                                                    convert(T, big"0.6067199374180168986519150843189505198519"),
                                                    convert(T, big"-0.02338660408468356531858175098561718651857"),
                                                    convert(T, big"-0.01798401877459493040442547470431484404443") )),
            # q4
            DerivativeCoefficientRow{T,1,7}(SVector(convert(T, big"0.1146397975178068401430112823144985150596"),
                                                    convert(T, big"-0.2898424301162697370942324201800071793273"),
                                                    convert(T, big"-0.3069262456316931913128086944558079603132"),
                                                    convert(T, 0),
                                                    convert(T, big"0.5203848121857539166740071338174418292578"),
                                                    convert(T, big"-0.05169127637022742348368508279860701098408"),
                                                    convert(T, big"0.01343534241462959507370778130248180630715") )),
            # q5
            DerivativeCoefficientRow{T,1,8}(SVector(convert(T, big"-0.03614399304268576976452921364705641609825"),
                                                    convert(T, big"0.1051508663818248421520867474440761344449"),
                                                    convert(T, big"0.01609777419666805778308369351834662756172"),
                                                    convert(T, big"-0.7080721616106272031118456849378369336023"),
                                                    convert(T, 0),
                                                    convert(T, big"0.7692160858661111736140494493705980473867"),
                                                    convert(T, big"-0.1645296432652024882569506157166433921544"),
                                                    convert(T, big"0.01828107147391138758410562396851593246160") )),
            # q6
            DerivativeCoefficientRow{T,1,9}(SVector(convert(T, big"-0.01141318406360863692889821914555232596651"),
                                                    convert(T, big"0.02049729840293952857599941220163960606616"),
                                                    convert(T, big"0.01113095018331244864875173213474522093204"),
                                                    convert(T, big"0.06324365883611076515355091406993789453750"),
                                                    convert(T, big"-0.6916640154753724474963890679085181638850"),
                                                    convert(T, 0),
                                                    convert(T, big"0.7397091390607520376247117645715851236273"),
                                                    convert(T, big"-0.1479418278121504075249423529143170247255"),
                                                    convert(T, big"0.01643798086801671194721581699047966941394") ))
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(3//4), T(-3//20), T(1//60))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(13649//43200),
                                T(12013//8640),
                                T(2711//4320),
                                T(5359//4320),
                                T(7877//8640),
                                T(43801//43200) )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, mode, 1, order, source)
    elseif order == 8
        left_boundary = (
            # q1
            DerivativeCoefficientRow{T,1,8}(SVector(convert(T, big"-1.695543604431898508749855654248370812054"),
                                                    convert(T, big"2.244525109969933844775276252674718713114"),
                                                    convert(T, big"-0.002163206607021876385102025416828692597556"),
                                                    convert(T, big"-0.8963677947307838167920190976769921451579"),
                                                    convert(T, big"0.2272938673777916084251861567415263328923"),
                                                    convert(T, big"0.1564962758424078514454013191006671546757"),
                                                    convert(T, big"0.003067629343740613524557473056482141287158"),
                                                    convert(T, big"-0.03730827676416971624344442423120269215938") )),
            # q2
            DerivativeCoefficientRow{T,1,8}(SVector(convert(T, big"-0.4338209217401506177055540526837828066711"),
                                                    convert(T, 0),
                                                    convert(T, big"0.06631961192627293470305188523936353462439"),
                                                    convert(T, big"0.6021342012739069461791324177355805603656"),
                                                    convert(T, big"-0.2213339967864188453859900127368926875520"),
                                                    convert(T, big"-0.001759980003751516130708214237001705688807"),
                                                    convert(T, big"-0.03586646217917856587352426166797879385310"),
                                                    convert(T, big"0.02432754750931966421359223835071189877503") )),
            # q3
            DerivativeCoefficientRow{T,1,8}(SVector(convert(T, big"0.002477771724790106958560398469919805908522"),
                                                    convert(T, big"-0.3930241559935857537756812928554282732813"),
                                                    convert(T, 0),
                                                    convert(T, big"-0.6080263371305341035810612202824357742461"),
                                                    convert(T, big"2.372533124692864412954923458457782022478"),
                                                    convert(T, big"-2.126112217612194556872952028715623246817"),
                                                    convert(T, big"0.9075309435559160631495435310307745215131"),
                                                    convert(T, big"-0.1553791292372561688333328461049890555549") )),
            # q4
            DerivativeCoefficientRow{T,1,8}(SVector(convert(T, big"0.1470043328582935680963279007411074464266"),
                                                    convert(T, big"-0.5109179516689331400658365102850902594996"),
                                                    convert(T, big"0.08705685833207777685654755900085198684605"),
                                                    convert(T, 0),
                                                    convert(T, big"-0.1306945333157231839535552435788086367322"),
                                                    convert(T, big"0.5981933016362239542598518794410035592720"),
                                                    convert(T, big"-0.2031099149801233772210520666270517105660"),
                                                    convert(T, big"0.01246790713818440202771648130798761425319") )),
            # q5
            DerivativeCoefficientRow{T,1,9}(SVector(convert(T, big"-0.1624073990846984662267508265053107632238"),
                                                    convert(T, big"0.8182390368133059538132603839842499379861"),
                                                    convert(T, big"-1.480018301574606037840376638130713134769"),
                                                    convert(T, big"0.5694185675497882973361521309100738568707"),
                                                    convert(T, 0),
                                                    convert(T, big"0.02929508293831017534661096502884570435894"),
                                                    convert(T, big"0.2825910005984090924748504462528740123097"),
                                                    convert(T, big"-0.04846434332860791750685068920563907354793"),
                                                    convert(T, big"-0.008653643911901097396895772334380539984709") )),
            # q6
            DerivativeCoefficientRow{T,1,10}(SVector(convert(T, big"-0.03609686950604370828405582087760384548297"),
                                                     convert(T, big"0.002100328577309696555612805397001149641337"),
                                                     convert(T, big"0.4281425817419204360479874008765051769418"),
                                                     convert(T, big"-0.8413238238875046736839550701643761633826"),
                                                     convert(T, big"-0.009456755727630000971085851751473603919416"),
                                                     convert(T, 0),
                                                     convert(T, big"0.5172389789041733771234394082779034212042"),
                                                     convert(T, big"-0.08760813565107717873901807278223424202774"),
                                                     convert(T, big"0.02979718129528502284256573906127239395909"),
                                                     convert(T, big"-0.002793485746432970891490538036994286933664") )),
             # q7
            DerivativeCoefficientRow{T,1,11}(SVector(convert(T, big"-0.0009797678134978722516935350416937185004514"),
                                                     convert(T, big"0.05926834509975463070197111976550139351983"),
                                                     convert(T, big"-0.2530570463899371286637621744353665227114"),
                                                     convert(T, big"0.3955555826583941989223787901885502942688"),
                                                     convert(T, big"-0.1263166266018192756531792392597193430222"),
                                                     convert(T, big"-0.7162192643577544899896537844517283192029"),
                                                     convert(T, 0),
                                                     convert(T, big"0.8209721963136350137518635528607308559584"),
                                                     convert(T, big"-0.2166153552278720352907291696202456084323"),
                                                     convert(T, big"0.04126006766245181624585317516576106827282"),
                                                     convert(T, big"-0.003868131343354857773048735171790100150577") )),
              # q8
            DerivativeCoefficientRow{T,1,12}(SVector(convert(T, big"0.01090012273307913185972171872891927027272"),
                                                     convert(T, big"-0.03677379943661633440187210478144113487019"),
                                                     convert(T, big"0.03963287609450569641558898819403594845401"),
                                                     convert(T, big"-0.02221139656912423686782339404035534965598"),
                                                     convert(T, big"0.01981665906734246925389947980969432794418"),
                                                     convert(T, big"0.1109698768905366566973791789571759073485"),
                                                     convert(T, big"-0.7509903604688148129224205834189298636300"),
                                                     convert(T, 0),
                                                     convert(T, big"0.7926019635554773751160111698352821779736"),
                                                     convert(T, big"-0.1981504908888693437790027924588205444934"),
                                                     convert(T, big"0.03774295064549892262457196046834677037969"),
                                                     convert(T, big"-0.003538401623015523996053621293907509723096") ))
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(4//5), T(-1//5), T(4//105), T(-1//280))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(1498139//5080320),
                                T(1107307//725760),
                                T(20761//80640),
                                T(1304999//725760),
                                T(299527//725760),
                                T(103097//80640),
                                T(670091//725760),
                                T(5127739//5080320) )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, mode, 1, order, source)
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end
