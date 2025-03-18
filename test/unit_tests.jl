using Test
using SummationByPartsOperators

@testset "DerivativeCoefficientRow" begin
    row1 = SummationByPartsOperators.DerivativeCoefficientRow{Rational, 1, 3}([1, 2, 3])
    row2 = SummationByPartsOperators.DerivativeCoefficientRow{Rational, 1, 3}([2, 4, 6])
    row3 = SummationByPartsOperators.DerivativeCoefficientRow{Rational, 2, 3}([2, 4, 6])

    @test -row1 ==
          SummationByPartsOperators.DerivativeCoefficientRow{Rational, 1, 3}([-1, -2, -3])

    @test row1 + row2 ==
          SummationByPartsOperators.DerivativeCoefficientRow{Rational, 1, 3}([3, 6, 9])

    @test row1 + row3 ==
          SummationByPartsOperators.DerivativeCoefficientRow{Rational, 1, 4}([1, 4, 7, 6])

    @test row1 / 2 ==
          SummationByPartsOperators.DerivativeCoefficientRow{Rational, 1, 3}([
                                                                                 1 // 2,
                                                                                 2 // 2,
                                                                                 3 // 2
                                                                             ])

    @test (row1 + -row3) / 2 ==
          SummationByPartsOperators.DerivativeCoefficientRow{Rational, 1, 4}([
                                                                                 1 // 2,
                                                                                 0,
                                                                                 -1 // 2,
                                                                                 -3
                                                                             ])
end
