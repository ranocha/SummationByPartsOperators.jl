
"""
    SharanBradyLivescu2022(alpha_left, alpha_right)

Coefficients of the cut-cell SBP operators given in
- Sharan, Brady, Livescu (2022)
  High-order dimensionally-split Cartesian embedded boundary method for non-dissipative schemes.
  Journal of Computational Physics 464, 111341.

Here, `alpha_left` * Δx is the spacing between the left endpoint and the second node, while
`alpha_right` * Δx is the spacing between the right endpoint and the last node.
"""
struct SharanBradyLivescu2022{T} <: SourceOfCoefficients
    alpha_left::T
    alpha_right::T
end

function Base.show(io::IO, source::SharanBradyLivescu2022)
    if get(io, :compact, false)
      summary(io, source)
    else
        print(io,
        "SharanBradyLivescu2022($(source.alpha_left), $(source.alpha_right)) \n",
        "   Coefficients of the cut-cell SBP operators given in \n",
        "    Sharan, Brady, Livescu (2022) \n",
        "    High-order dimensionally-split Cartesian embedded boundary method for non-dissipative schemes. \n",
        "    Journal of Computational Physics 464, 111341. \n",
        "Here, $(source.alpha_left) * Δx is the spacing between the left endpoint and the second node, while \n",
        "$(source.alpha_right) * Δx is the spacing between the right endpoint and the last node. ")
    end
end

# A grid with unequal spacing between the the first and second nodes and the
# second-to-last and last nodes.
struct CutCellGrid{T,Grid} <: AbstractArray{T,1}
    xmin::T
    xmax::T
    uniform_interior_grid::Grid # the interior nodes
end

function CutCellGrid(xmin, xmax, alpha_left, alpha_right, N::Int)
    @argcheck xmin < xmax
    @argcheck N > 3

    # build the grid on [-1 + α_left * Δx, 1 - α_right * Δx].
    Δx = (xmax - xmin) / (N - 3) / (1 + (alpha_left + alpha_right) / (N - 3))

    return CutCellGrid(xmin, xmax, range(xmin + alpha_left * Δx, xmax - alpha_right * Δx, length=N-2))
end

Base.length(grid::CutCellGrid) = length(grid.uniform_interior_grid) + 2

function Base.getindex(grid::CutCellGrid, i::Int)
    N = length(grid)
    @boundscheck begin
      @argcheck i > 0
      @argcheck i <= N
    end

    if i == 1
        return grid.xmin
    elseif i == N
        return grid.xmax
    else
      @inbounds x = getindex(grid.uniform_interior_grid, i - 1)
      return x
    end
end

Base.size(grid::CutCellGrid) = (length(grid),)
Base.step(grid::CutCellGrid) = step(grid.uniform_interior_grid)

function boundary_coefficients_weights_order_2(alpha, T=Float64; sign=true)
    coeffs = (
            DerivativeCoefficientRow{T, 1, 3}(sign * SVector(-(2 / (1 + alpha)), 1, (1 - alpha) / (1 + alpha))),
            DerivativeCoefficientRow{T, 1, 3}(sign * SVector(-(1 / (1 + alpha)), 0, 1 / (1 + alpha))),
            DerivativeCoefficientRow{T, 1, 4}(sign * SVector((1 - alpha) / (-4 - alpha + alpha^2),
                                                      (1 + alpha) / (-4 - alpha + alpha^2),
                                                      0, 2 / (4 + alpha - alpha^2))),
           )
    weights = SVector((1 + alpha) / 4, (1 // 4) * (1 + alpha)^2, (1 // 4) * (4 + alpha - alpha^2))
    return coeffs, weights
end


function boundary_coefficients_weights_order_4(alpha, T=Float64; sign=true)
    b1 = (-59 + 72 * alpha^2 + 40 * alpha^3 + 6 * alpha^4) /
            (144 * alpha + 120 * alpha^2 + 24 * alpha^3)
    b2 = (43 + 36 * alpha^2 + 32 * alpha^3 + 6 * alpha^4) /
            (72 * alpha + 96 * alpha^2 + 24 * alpha^3)
    b3 = (17 + 48 * alpha + 44 * alpha^2 + 16 * alpha^3 + 2 * alpha^4) /
            (48 + 88 * alpha + 48 * alpha^2 + 8 * alpha^3)
    s = 147 // 233

    h11 = s * max(zero(T), b1) + (1 - s) * min(b2, b3)

    coeffs = (
        DerivativeCoefficientRow{T, 1, 4}(sign * SVector(-1 / (2h11), (24 + 12(alpha^2) + 4alpha * (9 - 12h11) - 72h11) / (48h11),
                                                    (32h11 - 8(alpha^2) - 16alpha * (1 - 2h11))/(16h11),
                                                    (4(alpha^2) + 4alpha * (1 - 4h11) - 8h11)/(16h11))),
        DerivativeCoefficientRow{T, 1, 5}(sign * SVector((72h11 - 24 - 12(alpha^2) - 4alpha * (9 - 12h11))/(17 + 2(alpha^4) + alpha * (48 - 88h11) + (alpha^2) * (44 - 48h11) - 48h11 - 8(alpha^3) * (h11 - 2)),
                                                    0, (59 + 2(alpha^4) + alpha * (96 - 208h11) + (alpha^2) * (56 - 72h11) + 8(alpha^3) * (2 - h11) - 192h11)/(2(17 + 2(alpha^4) + alpha * (48 - 88h11) + (alpha^2) * (44 - 48h11) - 48h11 - 8(alpha^3) * (h11 - 2))),
                                                    (24h11 - 4 - 4(alpha^4) - 4alpha * (3 - 20h11) - (alpha^2) * (28 - 72h11) - 4(alpha^3) * (5 - 4h11))/(17 + 2(alpha^4) + alpha * (48 - 88h11) + (alpha^2) * (44 - 48h11) - 48h11 - 8(alpha^3) * (h11 - 2)),
                                                    (6(alpha^4) + 8(alpha^3) * (3 - 3h11) + 3(alpha^2) * (8 - 24h11) - 3 - 48alpha * h11)/(2(17 + 2(alpha^4) + alpha * (48 - 88h11) + (alpha^2) * (44 - 48h11) - 48h11 - 8(alpha^3) * (h11 - 2))))),
        DerivativeCoefficientRow{T, 1, 5}(sign * SVector((-3(32h11 - 8(alpha^2) - 16alpha * (1 - 2h11)))/(59 + 8(alpha^3) * (3h11 - 5) + 24(alpha^2) * (5h11 - 3) + 144alpha * h11 - 6(alpha^4)),
                                                    (192h11 - 59 - 2(alpha^4) - alpha * (96 - 208h11) - (alpha^2) * (56 - 72h11) - 8(alpha^3) * (2 - h11))/(2(59 + 8(alpha^3) * (3h11 - 5) + 24(alpha^2) * (5h11 - 3) + 144alpha * h11 - 6(alpha^4))),
                                                    0, (59 + 18(alpha^4) + 3(alpha^2) * (24 - 88h11) + 8(alpha^3) * (10 - 9h11) - 144alpha * h11)/(2(59 + 8(alpha^3) * (3h11 - 5) + 24(alpha^2) * (5h11 - 3) + 144alpha * h11 - 6(alpha^4))),
                                                    (alpha * (2 + alpha) * (32h11 - 8(alpha^2) - 16alpha * (1 - 2h11)))/(59 + 8(alpha^3) * (3h11 - 5) + 24(alpha^2) * (5h11 - 3) + 144alpha * h11 - 6(alpha^4)))),
        DerivativeCoefficientRow{T, 1, 6}(sign * SVector((3(8h11 - 4(alpha^2) - 4alpha * (1 - 4h11)))/(43 + 6(alpha^4) + (alpha^2) * (36 - 96h11) - 8(alpha^3) * (3h11 - 4) - 72alpha * h11),
                                                    (4 + 4(alpha^4) + (alpha^2) * (28 - 72h11) + 4alpha * (3 - 20h11) + 4(alpha^3) * (5 - 4h11) - 24h11)/(43 + 6(alpha^4) + (alpha^2) * (36 - 96h11) - 8(alpha^3) * (3h11 - 4) - 72alpha * h11),
                                                    (144alpha * h11 - 59 - 18(alpha^4) - 3(alpha^2) * (24 - 88h11) - 8(alpha^3) * (10 - 9h11))/(2(43 + 6(alpha^4) + (alpha^2) * (36 - 96h11) - 8(alpha^3) * (3h11 - 4) - 72alpha * h11)),
                                                    0, (59 + 10(alpha^4) + (alpha^2) * (40 - 120h11) + 8(alpha^3) * (5 - 5h11) - 80alpha * h11)/(2(43 + 6(alpha^4) + (alpha^2) * (36 - 96h11) - 8(alpha^3) * (3h11 - 4) - 72alpha * h11)),
                                                    -4/(43 + 6(alpha^4) + (alpha^2) * (36 - 96h11) - 8(alpha^3) * (3h11 - 4) - 72alpha * h11))),
        DerivativeCoefficientRow{T, 1, 7}(sign * SVector(0, (6(alpha^4) + 8(alpha^3) * (3 - 3h11) + 3(alpha^2) * (8 - 24h11) - 3 - 48alpha * h11)/(2(2(alpha^4) + (alpha^2) * (8 - 24h11) - 49 - 8(alpha^3) * (h11 - 1) - 16alpha * h11)),
                                                    (alpha * (2 + alpha) * (32h11 - 8(alpha^2) - 16alpha * (1 - 2h11)))/(2(alpha^4) + (alpha^2) * (8 - 24h11) - 49 - 8(alpha^3) * (h11 - 1) - 16alpha * h11),
                                                    (59 + 10(alpha^4) + (alpha^2) * (40 - 120h11) + 8(alpha^3) * (5 - 5h11) - 80alpha * h11)/(2(2(alpha^4) + (alpha^2) * (8 - 24h11) - 49 - 8(alpha^3) * (h11 - 1) - 16alpha * h11)),
                                                    0, -32/(2(alpha^4) + (alpha^2) * (8 - 24h11) - 49 - 8(alpha^3) * (h11 - 1) - 16alpha * h11),
                                                    4/(2(alpha^4) + (alpha^2) * (8 - 24h11) - 49 - 8(alpha^3) * (h11 - 1) - 16alpha * h11)))
       )

    weights = SVector(h11,
                     (1 // 48) * (17 + 2 * alpha^4 + alpha * (48 - 88 * h11) + alpha^2 * (44 - 48 * h11) - 8 * alpha^3 * (-2 + h11) - 48 * h11),
                     (1 // 48) * (59 - 6 * alpha^4 + 144 * alpha * h11 + 8 * alpha^3 * (-5 + 3 * h11) + 24 * alpha^2 * (-3 + 5 * h11)),
                     (1 // 48) * (43 + 6 * alpha^4 + alpha^2 * (36 - 96 * h11) - 72 * alpha * h11 - 8 * alpha^3 * (-4 + 3 * h11)),
                     (1 // 48) * (49 - 2 * alpha^4 + 8 * alpha^3 * (-1 + h11) + 16 * alpha * h11 + 8 * alpha^2 * (-1 + 3 * h11))
                    )

    return coeffs, weights
end

function boundary_coefficients_weights_order_6(alpha, T=Float64; sign=true)
    b1 = (-12013 +
          21600 * alpha^2 +
          18480 * alpha^3 +
          6390 * alpha^4 +
          1008 * alpha^5 +
          60 * alpha^6) /
         (360 * alpha * (120 + 154 * alpha + 71 * alpha^2 + 14 * alpha^3 + alpha^4))
    b2 = (2711 +
          10800 * alpha^2 +
          12840 * alpha^3 +
          5310 * alpha^4 +
          936 * alpha^5 +
          60 * alpha^6) /
         (360 * alpha * (60 + 107 * alpha + 59 * alpha^2 + 13 * alpha^3 + alpha^4))
    b3 = (13649 +
          43200 * alpha +
          49320 * alpha^2 +
          27000 * alpha^3 +
          7650 * alpha^4 +
          1080 * alpha^5 +
          60 * alpha^6) /
         (360 * (120 + 274 * alpha + 225 * alpha^2 + 85 * alpha^3 + 15 * alpha^4 + alpha^5))
    s = 16077 // 26921

    h11 = s * max(0, b1) + (1 - s) * min(b2, b3)
    d61 = 15025 // 525612

    coeffs = (
        DerivativeCoefficientRow{T, 1, 6}(sign * SVector(-1/(2h11), (21600 + 720alpha * (55 + 3h11 * (25d61 - 40)) + 600(alpha^3) * (6 + d61 * (123h11 - 61)) + 360(alpha^2) * (5d61 * (61h11 - 15) - 60(h11 - 1)) + 450d61 * (alpha^4) * (44h11 - 41) + 360d61 * (alpha^5) * (5h11 - 11) - 39385d61 - 79200h11 - 300d61 * (alpha^6))/(43200h11),
                                                    (25920h11 + 31508d61 + 240d61 * (alpha^6) - 720(alpha^2) * (3(5 - 6h11) + 2d61 * (61h11 - 15)) - 4320alpha * (3 + 2h11 * (5d61 - 5)) - 120(alpha^3) * (18 + 4d61 * (123h11 - 61)) - 360d61 * (alpha^4) * (44h11 - 41) - 288d61 * (alpha^5) * (5h11 - 11))/(8640h11),
                                                    (360(alpha^2) * (2(6 - 9h11) + 3d61 * (61h11 - 15)) + 120(alpha^3) * (9 + 3d61 * (123h11 - 61)) + 1080alpha * (3 + 2h11 * (15d61 - 8)) + 270d61 * (alpha^4) * (44h11 - 41) + 216d61 * (alpha^5) * (5h11 - 11) - 6480h11 - 23631d61 - 180d61 * (alpha^6))/(4320h11),
                                                    (1440h11 + 15754d61 + 120d61 * (alpha^6) - 720alpha * (1 + 6h11 * (5d61 - 1)) - 360(3 + 2d61 * (61h11 - 15) - 6h11) * (alpha^2) - 120(alpha^3) * (3 + 2d61 * (123h11 - 61)) - 90d61 * (alpha^4) * (88h11 - 82) - 144d61 * (alpha^5) * (5h11 - 11))/(4320h11),
                                                    (d61 * (72(alpha^5) * (5h11 - 11) + 120(alpha^3) * (123h11 - 61) + 90(alpha^4) * (44h11 - 41) + 360(alpha^2) * (61h11 - 15) + 10800alpha * h11 - 7877 - 60(alpha^6)))/(8640h11))),
        DerivativeCoefficientRow{T, 1, 6}(sign * SVector((21600 + 720alpha * (55 + 3h11 * (25d61 - 40)) + 600(alpha^3) * (6 + d61 * (123h11 - 61)) + 360(alpha^2) * (5d61 * (61h11 - 15) - 60(h11 - 1)) + 450d61 * (alpha^4) * (44h11 - 41) + 360d61 * (alpha^5) * (5h11 - 11) - 39385d61 - 79200h11 - 300d61 * (alpha^6))/(43200h11 + 360(alpha^2) * (225h11 - 137) + 1800(alpha^3) * (17h11 - 15) + 360(alpha^5) * (h11 - 3) + 450(alpha^4) * (12h11 - 17) + 720alpha * (137h11 - 60) - 13649 - 60(alpha^6)),
                                                    0, ((alpha^3) * (3600(173h11 - 90) - 5d61 * (1022400h11 - 323957)) + 2(388800h11 + 472620d61 - 81763) + 300d61 * (alpha^9) - 30(alpha^6) * (68 + d61 * (9960h11 - 10427)) - 60(2145 + 10d61 * (5727h11 - 2729) - 2340h11) * (alpha^4) - 9(alpha^2) * (5d61 * (89760h11 - 22277) - 240(635h11 - 213)) - 180(alpha^5) * (d61 * (7400h11 - 5173) - 4(17h11 - 36)) - 2alpha * (d61 * (648000h11 - 512005) - 6480(121h11 - 30)) - 180d61 * (alpha^8) * (10h11 - 37) - 30d61 * (alpha^7) * (1200h11 - 2063))/(6(43200h11 + 360(alpha^2) * (225h11 - 137) + 1800(alpha^3) * (17h11 - 15) + 360(alpha^5) * (h11 - 3) + 450(alpha^4) * (12h11 - 17) + 720alpha * (137h11 - 60) - 13649 - 60(alpha^6))),
                                                    (alpha * (-2160(309h11 - 30) - d61 * (748315 - 648000h11)) + (alpha^3) * (-600(999h11 - 465) - d61 * (991585 - 3403800h11)) + 30(996 + d61 * (36600h11 - 23929) - 492h11) * (alpha^5) + 30(alpha^4) * (15(301 - 348h11) + 20d61 * (4287h11 - 1888)) + 60(alpha^6) * (41 + 8d61 * (555h11 - 548)) + 30d61 * (alpha^7) * (1140h11 - 1861) + 120d61 * (alpha^8) * (15h11 - 53) - 131 - 129600h11 - 472620d61 - 300d61 * (alpha^9) - 8(alpha^2) * (45(2765h11 - 687) + d61 * (79885 - 292950h11)))/(2(43200h11 + 360(alpha^2) * (225h11 - 137) + 1800(alpha^3) * (17h11 - 15) + 360(alpha^5) * (h11 - 3) + 450(alpha^4) * (12h11 - 17) + 720alpha * (137h11 - 60) - 13649 - 60(alpha^6))),
                                                    ((alpha^3) * (3600(155h11 - 70) - 5d61 * (501120h11 - 142037)) + (alpha^2) * (720(1115h11 - 261) - 5d61 * (326880h11 - 98339)) + 2(9143 + 14400h11 + 157540d61) + 300d61 * (alpha^9) - 2alpha * (7200 + 5d61 * (43200h11 - 55139) - 209520h11) - 60(528 + d61 * (15300h11 - 9553) - 276h11) * (alpha^5) - 30(alpha^6) * (92 + d61 * (7920h11 - 7453)) - 60(2235 + 10d61 * (3357h11 - 1415) - 2700h11) * (alpha^4) - 30d61 * (alpha^7) * (1080h11 - 1679) - 60d61 * (alpha^8) * (30h11 - 101))/(2(43200h11 + 360(alpha^2) * (225h11 - 137) + 1800(alpha^3) * (17h11 - 15) + 360(alpha^5) * (h11 - 3) + 450(alpha^4) * (12h11 - 17) + 720alpha * (137h11 - 60) - 13649 - 60(alpha^6))),
                                                    (alpha * (5d61 * (64800h11 - 86647) - 304560h11) + (alpha^3) * (-1800(277h11 - 123) - d61 * (555985 - 1974600h11)) + 6(alpha^2) * (5d61 * (41760h11 - 13277) - 180(615h11 - 141)) + 30(4155 + 10d61 * (5478h11 - 2251) - 5220h11) * (alpha^4) + 90(alpha^5) * (d61 * (8680h11 - 5259) - 4(49h11 - 87)) + 60(alpha^6) * (49 + d61 * (3540h11 - 3211)) + 360d61 * (alpha^8) * (5h11 - 16) + 30d61 * (alpha^7) * (1020h11 - 1517) - 20539 - 236310d61 - 300d61 * (alpha^9))/(6(43200h11 + 360(alpha^2) * (225h11 - 137) + 1800(alpha^3) * (17h11 - 15) + 360(alpha^5) * (h11 - 3) + 450(alpha^4) * (12h11 - 17) + 720alpha * (137h11 - 60) - 13649 - 60(alpha^6))))),
        DerivativeCoefficientRow{T, 1, 7}(sign * SVector((720(alpha^2) * (3(5 - 6h11) + 2d61 * (61h11 - 15)) + 120(alpha^3) * (18 + 4d61 * (123h11 - 61)) + 4320alpha * (3 + 2h11 * (5d61 - 5)) + 360d61 * (alpha^4) * (44h11 - 41) + 288d61 * (alpha^5) * (5h11 - 11) - 25920h11 - 31508d61 - 240d61 * (alpha^6))/(12013 + 90(alpha^4) * (56h11 - 71) + 72(alpha^5) * (5h11 - 14) + 720(alpha^2) * (77h11 - 30) + 120(alpha^3) * (213h11 - 154) + 43200alpha * h11 - 60(alpha^6)),
                                                    ((alpha^3) * (3600(173h11 - 90) - 5d61 * (1022400h11 - 323957)) + 2(388800h11 + 472620d61 - 81763) + 300d61 * (alpha^9) - 30(alpha^6) * (68 + d61 * (9960h11 - 10427)) - 60(2145 + 10d61 * (5727h11 - 2729) - 2340h11) * (alpha^4) - 9(alpha^2) * (5d61 * (89760h11 - 22277) - 240(635h11 - 213)) - 180(alpha^5) * (d61 * (7400h11 - 5173) - 4(17h11 - 36)) - 2alpha * (d61 * (648000h11 - 512005) - 6480(121h11 - 30)) - 180d61 * (alpha^8) * (10h11 - 37) - 30d61 * (alpha^7) * (1200h11 - 2063))/(30(12013 + 90(alpha^4) * (56h11 - 71) + 72(alpha^5) * (5h11 - 14) + 720(alpha^2) * (77h11 - 30) + 120(alpha^3) * (213h11 - 154) + 43200alpha * h11 - 60(alpha^6))),
                                                    0, (14714 + alpha * (-449280h11 - 567144d61) + (alpha^3) * (6d61 * (339120h11 - 72677) - 240(3153h11 - 1399)) + 6(alpha^2) * (d61 * (129600h11 - 55139) - 480(370h11 - 78)) + 12(alpha^6) * (290 + 3d61 * (7800h11 - 7109)) + 18(alpha^5) * (40d61 * (1440h11 - 841) - 8(145h11 - 287)) + 60(2991 + 36d61 * (949h11 - 349) - 3552h11) * (alpha^4) + 72d61 * (alpha^8) * (30h11 - 101) + 108d61 * (alpha^7) * (360h11 - 553) - 360d61 * (alpha^9))/(6(12013 + 90(alpha^4) * (56h11 - 71) + 72(alpha^5) * (5h11 - 14) + 720(alpha^2) * (77h11 - 30) + 120(alpha^3) * (213h11 - 154) + 43200alpha * h11 - 60(alpha^6))),
                                                    (30637 + 480d61 * (alpha^9) + 2alpha * (252064d61 + 298080h11) - (alpha^3) * (8d61 * (240480h11 - 51077) - 360(2687h11 - 1198)) - 6(alpha^6) * (850 + 16d61 * (3450h11 - 2983)) - 36(1596 + 80d61 * (395h11 - 219) - 850h11) * (alpha^5) - 3(alpha^2) * (16d61 * (14400h11 - 7877) - 720(607h11 - 138)) - 30(7917 + 64d61 * (1086h11 - 379) - 9720h11) * (alpha^4) - 48d61 * (alpha^7) * (1020h11 - 1487) - 576d61 * (alpha^8) * (5h11 - 16))/(6(12013 + 90(alpha^4) * (56h11 - 71) + 72(alpha^5) * (5h11 - 14) + 720(alpha^2) * (77h11 - 30) + 120(alpha^3) * (213h11 - 154) + 43200alpha * h11 - 60(alpha^6))),
                                                    ((alpha^2) * (-720(223h11 - 51) - d61 * (39385 - 64800h11)) + alpha * (-47262d61 - 73440h11) + (alpha^3) * (-240(504h11 - 223) - d61 * (40277 - 185760h11)) + 60(504 + 6d61 * (581h11 - 197) - 636h11) * (alpha^4) + 6(alpha^6) * (120 + d61 * (6120h11 - 5087)) + 6(alpha^5) * (10d61 * (1992h11 - 1069) - 24(30h11 - 53)) + 90d61 * (alpha^7) * (64h11 - 89) + 12d61 * (alpha^8) * (30h11 - 91) - 4656 - 60d61 * (alpha^9))/(2(12013 + 90(alpha^4) * (56h11 - 71) + 72(alpha^5) * (5h11 - 14) + 720(alpha^2) * (77h11 - 30) + 120(alpha^3) * (213h11 - 154) + 43200alpha * h11 - 60(alpha^6))),
                                                    (6611 + 6600(alpha^3) * (21h11 - 10) + 4950(alpha^4) * (8h11 - 7) + 7920(alpha^2) * (25h11 - 6) + 3960(alpha^5) * (h11 - 2) + 95040alpha * h11 - 660(alpha^6))/(30(12013 + 90(alpha^4) * (56h11 - 71) + 72(alpha^5) * (5h11 - 14) + 720(alpha^2) * (77h11 - 30) + 120(alpha^3) * (213h11 - 154) + 43200alpha * h11 - 60(alpha^6))))),
        DerivativeCoefficientRow{T, 1, 7}(sign * SVector((6480h11 + 23631d61 + 180d61 * (alpha^6) - 360(alpha^2) * (2(6 - 9h11) + 3d61 * (61h11 - 15)) - 120(alpha^3) * (9 + 3d61 * (123h11 - 61)) - 1080alpha * (3 + 2h11 * (15d61 - 8)) - 270d61 * (alpha^4) * (44h11 - 41) - 216d61 * (alpha^5) * (5h11 - 11))/(2711 + 60(alpha^6) - 72(alpha^5) * (5h11 - 13) - 360(alpha^2) * (107h11 - 30) - 120(alpha^3) * (177h11 - 107) - 21600alpha * h11 - 90(alpha^4) * (52h11 - 59)),
                                                    (alpha * (d61 * (648000h11 - 748315) - 2160(309h11 - 30)) + (alpha^3) * (5d61 * (680760h11 - 198317) - 600(999h11 - 465)) + 8(alpha^2) * (5d61 * (58590h11 - 15977) - 45(2765h11 - 687)) + 30(996 + d61 * (36600h11 - 23929) - 492h11) * (alpha^5) + 30(alpha^4) * (15(301 - 348h11) + 20d61 * (4287h11 - 1888)) + 60(alpha^6) * (41 + 8d61 * (555h11 - 548)) + 30d61 * (alpha^7) * (1140h11 - 1861) + 120d61 * (alpha^8) * (15h11 - 53) - 131 - 129600h11 - 472620d61 - 300d61 * (alpha^9))/(20(2711 + 60(alpha^6) - 72(alpha^5) * (5h11 - 13) - 360(alpha^2) * (107h11 - 30) - 120(alpha^3) * (177h11 - 107) - 21600alpha * h11 - 90(alpha^4) * (52h11 - 59))),
                                                    ((alpha^3) * (240(3153h11 - 1399) - 6d61 * (339120h11 - 72677)) + alpha * (449280h11 + 567144d61) + 360d61 * (alpha^9) - 14714 - 6(alpha^2) * (d61 * (129600h11 - 55139) - 480(370h11 - 78)) - 12(alpha^6) * (290 + 3d61 * (7800h11 - 7109)) - 18(alpha^5) * (40d61 * (1440h11 - 841) - 8(145h11 - 287)) - 60(2991 + 36d61 * (949h11 - 349) - 3552h11) * (alpha^4) - 72d61 * (alpha^8) * (30h11 - 101) - 108d61 * (alpha^7) * (360h11 - 553))/(12(2711 + 60(alpha^6) - 72(alpha^5) * (5h11 - 13) - 360(alpha^2) * (107h11 - 30) - 120(alpha^3) * (177h11 - 107) - 21600alpha * h11 - 90(alpha^4) * (52h11 - 59))),
                                                    0, (1290 + (alpha^2) * (5d61 * (8640h11 - 7877) - 720(235h11 - 54)) + alpha * (-31508d61 - 77760h11) + (alpha^3) * (d61 * (141840h11 - 29477) - 120(1056h11 - 469)) + 6(alpha^5) * (20d61 * (930h11 - 473) - 3(240h11 - 436)) + 6(alpha^6) * (120 + d61 * (6000h11 - 4823)) + 60(525 + d61 * (2994h11 - 938) - 660h11) * (alpha^4) + 30d61 * (alpha^7) * (192h11 - 263) + 12d61 * (alpha^8) * (30h11 - 91) - 60d61 * (alpha^9))/(2(2711 + 60(alpha^6) - 72(alpha^5) * (5h11 - 13) - 360(alpha^2) * (107h11 - 30) - 120(alpha^3) * (177h11 - 107) - 21600alpha * h11 - 90(alpha^4) * (52h11 - 59))),
                                                    (11237 + (alpha^3) * (120(3909h11 - 1723) + d61 * (72231 - 327240h11)) + alpha * (70893d61 + 282960h11) + 180d61 * (alpha^9) - 12(alpha^6) * (235 + d61 * (7920h11 - 6114)) - 12(alpha^2) * (d61 * (8100h11 - 7877) - 2(25845h11 - 5895)) - 30(3909 + 36d61 * (397h11 - 121) - 4956h11) * (alpha^4) - 18(1652 + 5d61 * (3096h11 - 1525) - 940h11) * (alpha^5) - 54d61 * (alpha^7) * (300h11 - 391) - 72d61 * (alpha^8) * (15h11 - 43))/(12(2711 + 60(alpha^6) - 72(alpha^5) * (5h11 - 13) - 360(alpha^2) * (107h11 - 30) - 120(alpha^3) * (177h11 - 107) - 21600alpha * h11 - 90(alpha^4) * (52h11 - 59))),
                                                    (840(alpha^6) - 6974 - 8400(alpha^3) * (21h11 - 10) - 6300(alpha^4) * (8h11 - 7) - 10080(alpha^2) * (25h11 - 6) - 5040(alpha^5) * (h11 - 2) - 120960alpha * h11)/(20(2711 + 60(alpha^6) - 72(alpha^5) * (5h11 - 13) - 360(alpha^2) * (107h11 - 30) - 120(alpha^3) * (177h11 - 107) - 21600alpha * h11 - 90(alpha^4) * (52h11 - 59))))),
        DerivativeCoefficientRow{T, 1, 8}(sign * SVector((720alpha * (1 + 6h11 * (5d61 - 1)) + 120(alpha^3) * (3 + 2d61 * (123h11 - 61)) + 360(3 + 2d61 * (61h11 - 15) - 6h11) * (alpha^2) + 90d61 * (alpha^4) * (88h11 - 82) + 144d61 * (alpha^5) * (5h11 - 11) - 1440h11 - 15754d61 - 120d61 * (alpha^6))/(5359 + 720(alpha^2) * (39h11 - 10) + 90(alpha^4) * (48h11 - 49) + 14400alpha * h11 + 360(alpha^3) * (49h11 - 26) + 72(alpha^5) * (5h11 - 12) - 60(alpha^6)),
                                                    ((alpha^3) * (3600(155h11 - 70) - 5d61 * (501120h11 - 142037)) + (alpha^2) * (720(1115h11 - 261) - 5d61 * (326880h11 - 98339)) + 2(9143 + 14400h11 + 157540d61) + 300d61 * (alpha^9) - 2alpha * (7200 + 5d61 * (43200h11 - 55139) - 209520h11) - 60(528 + d61 * (15300h11 - 9553) - 276h11) * (alpha^5) - 30(alpha^6) * (92 + d61 * (7920h11 - 7453)) - 60(2235 + 10d61 * (3357h11 - 1415) - 2700h11) * (alpha^4) - 30d61 * (alpha^7) * (1080h11 - 1679) - 60d61 * (alpha^8) * (30h11 - 101))/(20(5359 + 720(alpha^2) * (39h11 - 10) + 90(alpha^4) * (48h11 - 49) + 14400alpha * h11 + 360(alpha^3) * (49h11 - 26) + 72(alpha^5) * (5h11 - 12) - 60(alpha^6))),
                                                    ((alpha^3) * (8d61 * (240480h11 - 51077) - 360(2687h11 - 1198)) + alpha * (-504128d61 - 596160h11) + 6(alpha^6) * (850 + 16d61 * (3450h11 - 2983)) + 36(1596 + 80d61 * (395h11 - 219) - 850h11) * (alpha^5) + 3(alpha^2) * (16d61 * (14400h11 - 7877) - 720(607h11 - 138)) + 30(7917 + 64d61 * (1086h11 - 379) - 9720h11) * (alpha^4) + 48d61 * (alpha^7) * (1020h11 - 1487) + 576d61 * (alpha^8) * (5h11 - 16) - 30637 - 480d61 * (alpha^9))/(12(5359 + 90(alpha^4) * (48h11 - 49) + 720(alpha^2) * (39h11 - 10) + 360(alpha^3) * (49h11 - 26) + 14400alpha * h11 + 72(alpha^5) * (5h11 - 12) - 60(alpha^6))),
                                                    ((alpha^2) * (720(235h11 - 54) - 5d61 * (8640h11 - 7877)) + (alpha^3) * (120(1056h11 - 469) - d61 * (141840h11 - 29477)) + alpha * (31508d61 + 77760h11) + 60d61 * (alpha^9) - 1290 - 6(alpha^5) * (20d61 * (930h11 - 473) - 3(240h11 - 436)) - 60(525 + d61 * (2994h11 - 938) - 660h11) * (alpha^4) - 6(alpha^6) * (120 + d61 * (6000h11 - 4823)) - 30d61 * (alpha^7) * (192h11 - 263) - 12d61 * (alpha^8) * (30h11 - 91))/(2(5359 + 90(alpha^4) * (48h11 - 49) + 720(alpha^2) * (39h11 - 10) + 360(alpha^3) * (49h11 - 26) + 14400alpha * h11 + 72(alpha^5) * (5h11 - 12) - 60(alpha^6))),
                                                    0, (27466 + (alpha^3) * (-720(377h11 - 169) - d61 * (18677 - 76320h11)) + 6(alpha^6) * (260 + d61 * (4560h11 - 3329)) + 60(1131 + 2d61 * (885h11 - 257) - 1404h11) * (alpha^4) + 36(468 + 5d61 * (412h11 - 193) - 260h11) * (alpha^5) + 6d61 * (alpha^7) * (840h11 - 1031) + 36d61 * (alpha^8) * (10h11 - 27) - 60d61 * (alpha^9) - 3(alpha^2) * (720(169h11 - 39) + d61 * (7877 - 7200h11)) - 2alpha * (7877d61 + 84240h11))/(12(5359 + 90(alpha^4) * (48h11 - 49) + 720(alpha^2) * (39h11 - 10) + 360(alpha^3) * (49h11 - 26) + 14400alpha * h11 + 72(alpha^5) * (5h11 - 12) - 60(alpha^6))),
                                                    (11400(alpha^3) * (21h11 - 10) + 8550(alpha^4) * (8h11 - 7) + 13680(alpha^2) * (25h11 - 6) + 6840(alpha^5) * (h11 - 2) + 164160alpha * h11 - 1541 - 1140(alpha^6))/(20(5359 + 90(alpha^4) * (48h11 - 49) + 720(alpha^2) * (39h11 - 10) + 360(alpha^3) * (49h11 - 26) + 14400alpha * h11 + 72(alpha^5) * (5h11 - 12) - 60(alpha^6))),
                                                    72/(5359 + 90(alpha^4) * (48h11 - 49) + 720(alpha^2) * (39h11 - 10) + 360(alpha^3) * (49h11 - 26) + 14400alpha * h11 + 72(alpha^5) * (5h11 - 12) - 60(alpha^6)))),
        DerivativeCoefficientRow{T, 1, 9}(sign * SVector(d61, (alpha * (5d61 * (64800h11 - 86647) - 304560h11) + (alpha^3) * (5d61 * (394920h11 - 111197) - 1800(277h11 - 123)) + 6(alpha^2) * (5d61 * (41760h11 - 13277) - 180(615h11 - 141)) + 30(4155 + 10d61 * (5478h11 - 2251) - 5220h11) * (alpha^4) + 90(alpha^5) * (d61 * (8680h11 - 5259) - 4(49h11 - 87)) + 60(alpha^6) * (49 + d61 * (3540h11 - 3211)) + 360d61 * (alpha^8) * (5h11 - 16) + 30d61 * (alpha^7) * (1020h11 - 1517) - 20539 - 236310d61 - 300d61 * (alpha^9))/(30(7877 + 60(alpha^6) - 72(alpha^5) * (5h11 - 11) - 10800alpha * h11 - 90(alpha^4) * (44h11 - 41) - 120(alpha^3) * (123h11 - 61) - 360(alpha^2) * (61h11 - 15))),
                                                    (4656 + (alpha^2) * (720(223h11 - 51) + d61 * (39385 - 64800h11)) + (alpha^3) * (240(504h11 - 223) + d61 * (40277 - 185760h11)) + alpha * (47262d61 + 73440h11) + 60d61 * (alpha^9) - 60(504 + 6d61 * (581h11 - 197) - 636h11) * (alpha^4) - 6(alpha^6) * (120 + d61 * (6120h11 - 5087)) - 6(alpha^5) * (10d61 * (1992h11 - 1069) - 24(30h11 - 53)) - 90d61 * (alpha^7) * (64h11 - 89) - 12d61 * (alpha^8) * (30h11 - 91))/(2(7877 + 60(alpha^6) - 72(alpha^5) * (5h11 - 11) - 120(alpha^3) * (123h11 - 61) - 90(alpha^4) * (44h11 - 41) - 360(alpha^2) * (61h11 - 15) - 10800alpha * h11)),
                                                    (alpha * (-70893d61 - 282960h11) + (alpha^3) * (d61 * (327240h11 - 72231) - 120(3909h11 - 1723)) + 12(alpha^6) * (235 + d61 * (7920h11 - 6114)) + 12(alpha^2) * (d61 * (8100h11 - 7877) - 2(25845h11 - 5895)) + 30(3909 + 36d61 * (397h11 - 121) - 4956h11) * (alpha^4) + 18(1652 + 5d61 * (3096h11 - 1525) - 940h11) * (alpha^5) + 54d61 * (alpha^7) * (300h11 - 391) + 72d61 * (alpha^8) * (15h11 - 43) - 11237 - 180d61 * (alpha^9))/(6(7877 + 60(alpha^6) - 72(alpha^5) * (5h11 - 11) - 120(alpha^3) * (123h11 - 61) - 90(alpha^4) * (44h11 - 41) - 360(alpha^2) * (61h11 - 15) - 10800alpha * h11)),
                                                    ((alpha^3) * (720(377h11 - 169) + d61 * (18677 - 76320h11)) + 60d61 * (alpha^9) + 3(alpha^2) * (720(169h11 - 39) + d61 * (7877 - 7200h11)) + 2alpha * (7877d61 + 84240h11) - 27466 - 6(alpha^6) * (260 + d61 * (4560h11 - 3329)) - 60(1131 + 2d61 * (885h11 - 257) - 1404h11) * (alpha^4) - 36(468 + 5d61 * (412h11 - 193) - 260h11) * (alpha^5) - 6d61 * (alpha^7) * (840h11 - 1031) - 36d61 * (alpha^8) * (10h11 - 27))/(6(7877 + 60(alpha^6) - 72(alpha^5) * (5h11 - 11) - 120(alpha^3) * (123h11 - 61) - 90(alpha^4) * (44h11 - 41) - 360(alpha^2) * (61h11 - 15) - 10800alpha * h11)),
                                                    0, (178774 + 1560(alpha^6) - 15600(alpha^3) * (21h11 - 10) - 11700(alpha^4) * (8h11 - 7) - 18720(alpha^2) * (25h11 - 6) - 9360(alpha^5) * (h11 - 2) - 224640alpha * h11)/(30(7877 + 60(alpha^6) - 72(alpha^5) * (5h11 - 11) - 120(alpha^3) * (123h11 - 61) - 90(alpha^4) * (44h11 - 41) - 360(alpha^2) * (61h11 - 15) - 10800alpha * h11)),
                                                    -1296/(7877 + 60(alpha^6) - 72(alpha^5) * (5h11 - 11) - 120(alpha^3) * (123h11 - 61) - 90(alpha^4) * (44h11 - 41) - 360(alpha^2) * (61h11 - 15) - 10800alpha * h11),
                                                    144/(7877 + 60(alpha^6) - 72(alpha^5) * (5h11 - 11) - 10800alpha * h11 - 90(alpha^4) * (44h11 - 41) - 120(alpha^3) * (123h11 - 61) - 360(alpha^2) * (61h11 - 15)))),
        DerivativeCoefficientRow{T, 1, 10}(sign * SVector(0, 0, (660(alpha^6) - 6611 - 6600(alpha^3) * (21h11 - 10) - 4950(alpha^4) * (8h11 - 7) - 7920(alpha^2) * (25h11 - 6) - 3960(alpha^5) * (h11 - 2) - 95040alpha * h11)/(6(43801 + 600(alpha^3) * (21h11 - 10) + 450(alpha^4) * (8h11 - 7) + 720(alpha^2) * (25h11 - 6) + 360(alpha^5) * (h11 - 2) + 8640alpha * h11 - 60(alpha^6))),
                                                    (6974 + 8400(alpha^3) * (21h11 - 10) + 6300(alpha^4) * (8h11 - 7) + 10080(alpha^2) * (25h11 - 6) + 5040(alpha^5) * (h11 - 2) + 120960alpha * h11 - 840(alpha^6))/(2(43801 + 600(alpha^3) * (21h11 - 10) + 450(alpha^4) * (8h11 - 7) + 720(alpha^2) * (25h11 - 6) + 360(alpha^5) * (h11 - 2) + 8640alpha * h11 - 60(alpha^6))),
                                                    (1541 + 1140(alpha^6) - 11400(alpha^3) * (21h11 - 10) - 8550(alpha^4) * (8h11 - 7) - 13680(alpha^2) * (25h11 - 6) - 6840(alpha^5) * (h11 - 2) - 164160alpha * h11)/(2(43801 + 600(alpha^3) * (21h11 - 10) + 450(alpha^4) * (8h11 - 7) + 720(alpha^2) * (25h11 - 6) + 360(alpha^5) * (h11 - 2) + 8640alpha * h11 - 60(alpha^6))),
                                                    (15600(alpha^3) * (21h11 - 10) + 11700(alpha^4) * (8h11 - 7) + 18720(alpha^2) * (25h11 - 6) + 9360(alpha^5) * (h11 - 2) + 224640alpha * h11 - 178774 - 1560(alpha^6))/(6(43801 + 600(alpha^3) * (21h11 - 10) + 450(alpha^4) * (8h11 - 7) + 720(alpha^2) * (25h11 - 6) + 360(alpha^5) * (h11 - 2) + 8640alpha * h11 - 60(alpha^6))),
                                                    0, 32400/(43801 + 600(alpha^3) * (21h11 - 10) + 450(alpha^4) * (8h11 - 7) + 720(alpha^2) * (25h11 - 6) + 360(alpha^5) * (h11 - 2) + 8640alpha * h11 - 60(alpha^6)),
                                                    -6480/(43801 + 600(alpha^3) * (21h11 - 10) + 450(alpha^4) * (8h11 - 7) + 720(alpha^2) * (25h11 - 6) + 360(alpha^5) * (h11 - 2) + 8640alpha * h11 - 60(alpha^6)),
                                                    720/(43801 + 600(alpha^3) * (21h11 - 10) + 450(alpha^4) * (8h11 - 7) + 720(alpha^2) * (25h11 - 6) + 360(alpha^5) * (h11 - 2) + 8640alpha * h11 - 60(alpha^6)))),
    )

    weights = SVector(
                        h11,
                        13649 // 43200 +
                        alpha +
                        alpha^6 / 720 +
                        alpha^2 * (137 // 120 - (15 * h11) / 8) +
                        alpha^3 * (5 // 8 - (17 * h11) / 24) +
                        alpha^4 * (17 // 96 - h11 / 8) - (1 // 120) * alpha^5 * (-3 + h11) - h11 -
                        (137 / 60) * alpha * h11,
                        (12013 - 60 * alpha^6 +
                        43200 * alpha * h11 +
                        72 * alpha^5 * (-14 + 5 * h11) +
                        90 * alpha^4 * (-71 + 56 * h11) +
                        720 * alpha^2 * (-30 + 77 * h11) +
                        120 * alpha^3 * (-154 + 213 * h11)) / 8640,
                        (2711 + 60 * alpha^6 - 21600 * alpha * h11 - 72 * alpha^5 * (-13 + 5 * h11) -
                        90 * alpha^4 * (-59 + 52 * h11) - 360 * alpha^2 * (-30 + 107 * h11) -
                        120 * alpha^3 * (-107 + 177 * h11)) / 4320,
                        5359 // 4320 - alpha^6 / 72 +
                        alpha^4 * (-(49 // 48) + h11) +
                        (10 // 3) * alpha * h11 +
                        (1 // 60) * alpha^5 * (-12 + 5 * h11) +
                        (1 // 6) * alpha^2 * (-10 + 39 * h11) +
                        (1 // 12) * alpha^3 * (-26 + 49 * h11),
                        (7877 + 60 * alpha^6 - 10800 * alpha * h11 - 72 * alpha^5 * (-11 + 5 * h11) -
                        90 * alpha^4 * (-41 + 44 * h11) - 360 * alpha^2 * (-15 + 61 * h11) -
                        120 * alpha^3 * (-61 + 123 * h11)) / 8640,
                        (43801 - 60 * alpha^6 +
                        360 * alpha^5 * (-2 + h11) +
                        8640 * alpha * h11 +
                        450 * alpha^4 * (-7 + 8 * h11) +
                        600 * alpha^3 * (-10 + 21 * h11) +
                        720 * alpha^2 * (-6 + 25 * h11)) / 43200,
                     )

    return coeffs, weights
end


function first_derivative_coefficients(source::SharanBradyLivescu2022, order::Int,
                                       T = Float64, mode = FastMode())
    @unpack alpha_left, alpha_right = source

    if order == 2
        left_boundary, left_weights = boundary_coefficients_weights_order_2(alpha_left, T)
        right_boundary, right_weights = boundary_coefficients_weights_order_2(alpha_right, T, sign = -1)
        upper_coef = SVector(T(1//2))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, mode, 1, order, source)
    elseif order == 4
        left_boundary, left_weights = boundary_coefficients_weights_order_4(alpha_left, T)
        right_boundary, right_weights = boundary_coefficients_weights_order_4(alpha_right, T; sign = -1)
        upper_coef = SVector(T(2//3), T(-1//12))
        central_coef = zero(T)
        lower_coef = -upper_coef
        left_boundary_derivatives = Tuple{}()
        right_boundary_derivatives = left_boundary_derivatives

        DerivativeCoefficients(left_boundary, right_boundary,
                                left_boundary_derivatives, right_boundary_derivatives,
                                lower_coef, central_coef, upper_coef,
                                left_weights, right_weights, mode, 1, order, source)
    elseif order == 6
        left_boundary, left_weights = boundary_coefficients_weights_order_6(alpha_left, T)
        right_boundary, right_weights = boundary_coefficients_weights_order_6(alpha_right, T; sign = -1)
        upper_coef = SVector(T(3//4), T(-3//20), T(1//60))
        central_coef = zero(T)
        lower_coef = -upper_coef
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
