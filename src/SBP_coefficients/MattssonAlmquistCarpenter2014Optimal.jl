
"""
    MattssonAlmquistCarpenter2014Optimal

Coefficients of the optimal SBP operators with nonuniform grid given in
  Mattsson, Almquist, Carpenter (2014)
  Optimal diagonal-norm SBP operators.
  Journal of Computational Physics 264, pp. 91-111.
"""
struct MattssonAlmquistCarpenter2014Optimal <: SourceOfCoefficients end

function Base.show(io::IO, ::MattssonAlmquistCarpenter2014Optimal)
    print(io,
        "  Mattsson, Almquist, Carpenter (2014) \n",
        "  Optimal diagonal-norm SBP operators. \n",
        "  Journal of Computational Physics 264, pp. 91-111. \n")
end



struct MattssonAlmquistCarpenter2014OptimalGrid{T,Grid} <: AbstractArray{T,1}
    xmin::T
    xmax::T
    d1::T
    d2::T
    d3::T
    uniform_grid::Grid
end

function MattssonAlmquistCarpenter2014OptimalGrid(xmin::T, xmax::T, d1::T, d2::T, d3::T, N::Int) where T
    @argcheck xmin < xmax
    @argcheck N > 6
    @argcheck d1 > 0
    @argcheck d2 > 0
    @argcheck d3 > 0

    uniform_grid = range(xmin+d1+d2+d3, stop=xmax-d1-d2-d3, length=N-6)
    MattssonAlmquistCarpenter2014OptimalGrid{T,typeof(uniform_grid)}(xmin, xmax, d1, d2, d3, uniform_grid)
end

function MattssonAlmquistCarpenter2014OptimalGrid(xmin, xmax, d1, d2, d3, N)
    MattssonAlmquistCarpenter2014OptimalGrid(promote(xmin, xmax, d1, d2, d3)..., N)
end

function Base.length(grid::MattssonAlmquistCarpenter2014OptimalGrid)
    length(grid.uniform_grid) + 6
end

function Base.getindex(grid::MattssonAlmquistCarpenter2014OptimalGrid, i::Int)
    N = length(grid)
    @boundscheck begin
        @argcheck i > 0
        @argcheck i <= N
    end

    if i == 1
        grid.xmin
    elseif i == 2
        grid.xmin + grid.d1
    elseif i == 3
        grid.xmin + grid.d1 + grid.d2
    elseif i == N
        grid.xmax
    elseif i == N-1
        grid.xmax - grid.d1
    elseif i == N-2
        grid.xmax - grid.d1 - grid.d2
    else
        @inbounds x = getindex(grid.uniform_grid, i-3)
        x
    end
end

Base.size(grid::MattssonAlmquistCarpenter2014OptimalGrid) = (length(grid),)
Base.step(grid::MattssonAlmquistCarpenter2014OptimalGrid) = step(grid.uniform_grid)


function construct_grid(::MattssonAlmquistCarpenter2014Optimal, accuracy_order, xmin, xmax, N)
    @argcheck N > 6
    T = promote_type(typeof(xmin), typeof(xmax))

    if accuracy_order == 2
        d1 = T(78866488858096586513//Int128(10)^20)
        d2 = T(95915098594220826013//Int128(10)^20)
        d3 = T(1)
    elseif accuracy_order == 4
        d1 = T(72181367003646814327//Int128(10)^20)
        d2 = T(13409118421582217252//Int128(10)^19)
        d3 = T(12898797485951900258//Int128(10)^19)
    elseif accuracy_order == 6
        d1 = T(51670081689316731234//Int128(10)^20)
        d2 = T(98190527037374634269//Int128(10)^20)
        d3 = T(10868393364992957832//Int128(10)^19)
    elseif accuracy_order == 8
        d1 = T(41669687672575697416//Int128(10)^20)
        d2 = T(78703773886730090312//Int128(10)^20)
        d3 = T(92685925671601406028//Int128(10)^20)
    else
        throw(ArgumentError("Accuracy order $accuracy_order not implemented/derived."))
    end

    h = (xmax-xmin) / (2*(d1+d2+d3) + N - 7)

    MattssonAlmquistCarpenter2014OptimalGrid(xmin, xmax, d1*h, d2*h, d3*h, N)
end


function first_derivative_coefficients(source::MattssonAlmquistCarpenter2014Optimal, order::Int, T=Float64, parallel=Val{:serial}())
    if order == 2
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,3}(SVector(T(-50000000000000000000//33743097329453577701),
                                                    T(55932483188770411252//33743097329453577701),
                                                    T(-5932483188770411252//33743097329453577701) )),
            # d2
            DerivativeCoefficientRow{T,1,3}(SVector(T(-13983120797192602813//24439920504708372824),
                                                    T(0),
                                                    T(13983120797192602813//24439920504708372824) )),
            # d3
            DerivativeCoefficientRow{T,1,4}(SVector(T(2966241594385205626//46639404052015171765),
                                                    T(-27966241594385205626//46639404052015171765),
                                                    T(0),
                                                    T(5000000000000000000//9327880810403034353) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(1//2))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(33743097329453577701//Int128(10)^20),
                                T(97759682018833491296//Int128(10)^20),
                                T(93278808104030343530//Int128(10)^20) )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 1, order, source)
    elseif order == 4
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,5}(SVector(T(-50000000000000000000//21427296612044126417),
                                                    T(66884898686930380508//21427296612044126417),
                                                    T(-25171531878753856238//21427296612044126417),
                                                    T(10997619816825822803//21427296612044126417),
                                                    T(-13554933125011735296//107136483060220632085) )),
            # d2
            DerivativeCoefficientRow{T,1,5}(SVector(T(-16721224671732595127//28093989712218483700),
                                                    T(0),
                                                    T(92214436948640491071//112375958848873934800),
                                                    T(-16206184326771260201//56187979424436967400),
                                                    T(17707075979581024571//280939897122184837000) )),
            # d3
            DerivativeCoefficientRow{T,1,5}(SVector(T(12585765939376928119//71722939624706300000),
                                                    T(-92214436948640491071//143445879249412600000),
                                                    T(0),
                                                    T(1636075617843355867//2868917584988252000),
                                                    T(-14760875822281158529//143445879249412600000) )),
            # d4
            DerivativeCoefficientRow{T,1,6}(SVector(T(-10997619816825822803//109173230217361308360),
                                                    T(16206184326771260201//54586615108680654180),
                                                    T(-8180378089216779335//10917323021736130836),
                                                    T(0),
                                                    T(17180591347196107273//27293307554340327090),
                                                    T(-625000000000000000//8187992266302098127) )),
            # d5
            DerivativeCoefficientRow{T,1,7}(SVector(T(2259155520835289216//82365134292746683125),
                                                    T(-17707075979581024571//247095402878240049375),
                                                    T(14760875822281158529//98838161151296019750),
                                                    T(-34361182694392214546//49419080575648009875),
                                                    T(0),
                                                    T(800000000000000000//1186057933815552237),
                                                    T(-100000000000000000//1186057933815552237) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(2//3), T(-1//12))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(21427296612044126417//100000000000000000000),
                                T(280939897122184837//250000000000000000),
                                T(717229396247063//500000000000000),
                                T(2729330755434032709//2500000000000000000),
                                T(395352644605184079//400000000000000000) )
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
            DerivativeCoefficientRow{T,1,7}(SVector(T(-3125000000000000000//944357158252257333),
                                                    T(22223596967296279011//5036571510678705776),
                                                    T(-5854697895099786871//3777428633009029332),
                                                    T(28083754862953532289//50365715106787057760),
                                                    T(-1592329083817967435//15109714532036117328),
                                                    T(-15653772860347171721//1510971453203611732800),
                                                    T(23913677348563239189//5036571510678705776000) )),
            # d2
            DerivativeCoefficientRow{T,1,7}(SVector(T(-22223596967296279011//26989195119035671001),
                                                    T(0),
                                                    T(89405599296515541581//80967585357107013003),
                                                    T(-28597427787314667440//80967585357107013003),
                                                    T(57056178538117177397//809675853571070130030),
                                                    T(41320613074890940489//8096758535710701300300),
                                                    T(-5124091837453295329//1619351707142140260060) )),
            # d3
            DerivativeCoefficientRow{T,1,7}(SVector(T(5854697895099786871//27278567870198137125),
                                                    T(-89405599296515541581//109114271480792548500),
                                                    T(0),
                                                    T(27653905086569037761//36371423826930849500),
                                                    T(-6077915680998146409//36371423826930849500),
                                                    T(83784382166533084621//10911427148079254850000),
                                                    T(42099567773838744673//10911427148079254850000) )),
            # d4
            DerivativeCoefficientRow{T,1,7}(SVector(T(-84251264588860596867//1043526904157157775600),
                                                    T(714935694682866686//2608817260392894439),
                                                    T(-82961715259707113283//104352690415715777560),
                                                    T(0),
                                                    T(75419218459746682761//104352690415715777560),
                                                    T(-14034899831339963049//104352690415715777560),
                                                    T(3512738257179466361//260881726039289443900) )),
            # d5
            DerivativeCoefficientRow{T,1,8}(SVector(T(1592329083817967435//98680905919946100728),
                                                    T(-57056178538117177397//986809059199461007280),
                                                    T(18233747042994439227//98680905919946100728),
                                                    T(-75419218459746682761//98680905919946100728),
                                                    T(0),
                                                    T(18687868497479752801//24670226479986525182),
                                                    T(-15119380469839682133//98680905919946100728),
                                                    T(625000000000000000//37005339719979787773) )),
            # d6
            DerivativeCoefficientRow{T,1,9}(SVector(T(15653772860347171721//10037581831426163456000),
                                                    T(-41320613074890940489//10037581831426163456000),
                                                    T(-29162680879405877//3493763254934272000),
                                                    T(14034899831339963049//100375818314261634560),
                                                    T(-18687868497479752801//25093954578565408640),
                                                    T(0),
                                                    T(4696526232232696867//6273488644641352160),
                                                    T(-5859375000000000//39209304029008451),
                                                    T(1953125000000000//117627912087025353) )),
            # d7
            DerivativeCoefficientRow{T,1,10}(SVector(T(-71741032045689717567//99943556356761752125000),
                                                    T(5124091837453295329//1998871127135235042500),
                                                    T(-42099567773838744673//9994355635676175212500),
                                                    T(-7025476514358932722//499717781783808760625),
                                                    T(15119380469839682133//99943556356761752125),
                                                    T(-75144419715723149872//99943556356761752125),
                                                    T(0),
                                                    T(600000000000000000//799548450854094017),
                                                    T(-120000000000000000//799548450854094017),
                                                    T(40000000000000000//2398645352562282051) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(3//4), T(-3//20), T(1//60))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(944357158252257333//6250000000000000000),
                                T(80967585357107013003//100000000000000000000),
                                T(218228542961585097//200000000000000000),
                                T(2608817260392894439//2500000000000000000),
                                T(12335113239993262591//12500000000000000000),
                                T(39209304029008451//39062500000000000),
                                T(799548450854094017//800000000000000000) )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 1, order, source)
    #= FIXME: There seem to be some errors in the coefficients...
    elseif order == 8
        left_boundary = (
            # d1
            DerivativeCoefficientRow{T,1,8}(SVector(T(-25000000000000000000//6081611055353751439),
                                                    T(66670790901888837033//12163222110707502878),
                                                    T(-10997015095317519523//6081611055353751439),
                                                    T(61752567584332553851//121632221107075028780),
                                                    T(-32312350944133128873//1216322211070750287800),
                                                    T(-16967490160001675093//608161105535375143900),
                                                    T(3031405594112644741//24326444221415005756000),
                                                    T(32350133072942893419//12163222110707502878000) )),
            # d2
            DerivativeCoefficientRow{T,1,8}(SVector(T(-1257939450979034661//1230864766727295094),
                                                    T(0),
                                                    T(86688767821045233147//65235832636546639982),
                                                    T(-24298087640343350527//65235832636546639982),
                                                    T(39549469619698650847//652358326365466399820),
                                                    T(2076352837148473751//652358326365466399820),
                                                    T(-4065343087310119561//4077239539784164998875),
                                                    T(-8167205603296830489//13047166527309327996400) )),
            # d3
            DerivativeCoefficientRow{T,1,8}(SVector(T(10997015095317519523//43865207099050505477),
                                                    T(-86688767821045233147//87730414198101010954),
                                                    T(0),
                                                    T(41032546292236417573//43865207099050505477),
                                                    T(-21014872891771196683//87730414198101010954),
                                                    T(38572177503610408523//877304141981010109540),
                                                    T(-16637807199883547459//87730414198101010954000),
                                                    T(-11580435163412724923//219326035495252527385000) )),
            # d4
            DerivativeCoefficientRow{T,1,8}(SVector(T(-61752567584332553851//973889517710795427990),
                                                    T(24298087640343350527//97388951771079542799),
                                                    T(-82065092584472835146//97388951771079542799),
                                                    T(0),
                                                    T(81102837324727866266//97388951771079542799),
                                                    T(-20423896795865859484//97388951771079542799),
                                                    T(776602705809652689//21641989282462120622),
                                                    T(1285505569126131521//541049732061553015550) )),
            # d5
            DerivativeCoefficientRow{T,1,9}(SVector(T(32312350944133128873//10072514376844677230000),
                                                    T(-39549469619698650847//1007251437684467723000),
                                                    T(21014872891771196683//100725143768446772300),
                                                    T(-40551418662363933133//50362571884223386150),
                                                    T(0),
                                                    T(3202162603743798001//4029005750737870892),
                                                    T(-1969916799269047253//10072514376844677230),
                                                    T(37220336417235830241//1007251437684467723000),
                                                    T(-25000000000000000//7050760063791274061) )),
            # d6
            DerivativeCoefficientRow{T,1,10}(SVector(T(16967490160001675093//4988436332888823941700),
                                                    T(-296621833878353393//142526752368252112620),
                                                    T(-38572177503610408523//997687266577764788340),
                                                    T(10211948397932929742//49884363328888239417),
                                                    T(-80054065093594950025//99768726657776478834),
                                                    T(0),
                                                    T(39951539656523293130//49884363328888239417),
                                                    T(-2222198748535843633//11085414073086275426),
                                                    T(40000000000000000000//1047571629906653027757),
                                                    T(-1250000000000000000//349190543302217675919) )),
            # d7
            DerivativeCoefficientRow{T,1,11}(SVector(T(-3031405594112644741//200106059975821710280000),
                                                    T(4065343087310119561//6253314374244428446250),
                                                    T(16637807199883547459//100053029987910855140000),
                                                    T(-6989424352286874201//200106059975821710280),
                                                    T(1969916799269047253//10005302998791085514),
                                                    T(-3995153965652329313//5002651499395542757),
                                                    T(0),
                                                    T(40008167342759928887//50026514993955427570),
                                                    T(-1000000000000000000//5002651499395542757),
                                                    T(4000000000000000000//105055681487306397897),
                                                    T(-125000000000000000//35018560495768799299) )),
            # d8
            DerivativeCoefficientRow{T,1,12}(SVector(T(-32350133072942893419//99994066100338390832000),
                                                    T(8167205603296830489//19998813220067678166400),
                                                    T(11580435163412724923//249985165250845977080000),
                                                    T(-11569550122135183689//4999703305016919541600),
                                                    T(-37220336417235830241//999940661003383908320),
                                                    T(19999788736822592697//99994066100338390832),
                                                    T(-40008167342759928887//49997033050169195416),
                                                    T(0),
                                                    T(5000000000000000000//6249629131271149427),
                                                    T(-1250000000000000000//6249629131271149427),
                                                    T(5000000000000000000//131242211756694137967),
                                                    T(-156250000000000000//43747403918898045989) )),
        )
        right_boundary = .- left_boundary
        upper_coef = SVector(T(4//5), T(-1//5), T(4//105), T(-1//280))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_weights = SVector( T(6081611055353751439//50000000000000000000),
                                T(32617916318273319991//50000000000000000000),
                                T(43865207099050505477//50000000000000000000),
                                T(97388951771079542799//100000000000000000000),
                                T(1007251437684467723//1000000000000000000),
                                T(49884363328888239417//50000000000000000000),
                                T(5002651499395542757//5000000000000000000),
                                T(6249629131271149427//6250000000000000000) )
        right_weights = left_weights
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, parallel, 1, order, source)
    =#
    else
        throw(ArgumentError("Order $order not implemented/derived."))
    end
end
