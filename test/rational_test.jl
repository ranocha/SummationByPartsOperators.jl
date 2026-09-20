# Tests of the order of accuracy and of the SBP (symmetry) properties of the
# operators using exact rational arithmetic. Since all relations checked here
# hold exactly for the exact coefficients, they can be checked with `==`
# instead of `≈`, which makes these tests much stronger than the corresponding
# floating point tests.
#
# `Rational{Int}` overflows quickly for the more complicated operators, so
# `Rational{BigInt}` is used throughout. It is slower but always works.
#
# Operators whose coefficients are only available as (truncated) floating point
# numbers cannot be tested here. These are
# - `DienerDorbandSchnetterTiglio2007` with accuracy orders 6 and 8
# - `MattssonAlmquistCarpenter2014Optimal` (the grid points are truncated decimals)
# - `MattssonAlmquistVanDerWeide2018Minimal`, `MattssonAlmquistVanDerWeide2018Accurate`
# - `MattssonNiemeläWinters2026`
# - `Mattsson2012` with accuracy order 6
# - `LegendreDerivativeOperator`, `FourierDerivativeOperator`, and
#   `function_space_operator`s

using Test
using LinearAlgebra
using SummationByPartsOperators

const RT = Rational{BigInt}

# Domain and number of nodes used for all tests below. The grid spacing
# `Δx = 5 // 64` is not equal to unity to make sure that the scaling of the
# coefficients by powers of `Δx` is tested as well.
const XMIN = RT(-1 // 2)
const XMAX = RT(2)
const NNODES = 33

"""
    derivative_of_monomial(x, k, d)

The `d`-th derivative of the monomial `x^k`.
"""
function derivative_of_monomial(x, k, d)
    k < d && return zero(x)
    factor = one(x)
    for j in 0:(d - 1)
        factor *= (k - j)
    end
    return factor * x^(k - d)
end

"""
    exactness_degrees(A, x, der_order, max_degree)

For every row of the matrix `A` representing an approximation of the
`der_order`-th derivative on the grid `x`, compute the maximal degree
`p ≤ max_degree` such that all monomials `x^k` with `k ≤ p` are differentiated
exactly. A value of `-1` means that not even constants are handled exactly.
"""
function exactness_degrees(A::AbstractMatrix, x::AbstractVector, der_order::Integer,
                           max_degree::Integer)
    n = size(A, 1)
    degrees = fill(-1, n)
    active = trues(n)
    for k in 0:max_degree
        residual = A * (x .^ k)
        for i in 1:n
            active[i] || continue
            if residual[i] == derivative_of_monomial(x[i], k, der_order)
                degrees[i] = k
            else
                active[i] = false
            end
        end
        any(active) || break
    end
    return degrees
end

"""
    test_exactness(D; interior, boundary)

Check that the interior stencils of the nonperiodic operator `D` differentiate
monomials up to degree `interior` (and no more) exactly and that the boundary
closures do so up to degree `boundary` (and no more).
"""
function test_exactness(D; interior, boundary)
    x = collect(grid(D))
    degrees = exactness_degrees(Matrix(D), x, derivative_order(D), interior + 1)
    nleft = length(D.coefficients.left_boundary)
    nright = length(D.coefficients.right_boundary)
    @test all(==(interior), degrees[(nleft + 1):(end - nright)])
    @test minimum(degrees[1:nleft]) == boundary
    @test minimum(degrees[(end - nright + 1):end]) == boundary
    return nothing
end

"""
    test_periodic_exactness(D; degree)

Check that the stencils of the periodic operator `D` differentiate monomials up
to degree `degree` (and no more) exactly. Rows with a stencil wrapping around
the periodic boundary are skipped since monomials are not periodic.
"""
function test_periodic_exactness(D; degree)
    x = collect(grid(D))
    degrees = exactness_degrees(Matrix(D), x, derivative_order(D), degree + 1)
    nlower = length(D.coefficients.lower_coef)
    nupper = length(D.coefficients.upper_coef)
    @test all(==(degree), degrees[(nlower + 1):(end - nupper)])
    return nothing
end

"""
    test_quadrature_exactness(D; degree)

Check that the quadrature rule given by the mass matrix of `D` is exact for
polynomials up to degree `degree` (and no more).
"""
function test_quadrature_exactness(D; degree)
    x = collect(grid(D))
    M = mass_matrix(D)
    xmin, xmax = first(x), last(x)
    for k in 0:degree
        @test sum(M * (x .^ k)) == (xmax^(k + 1) - xmin^(k + 1)) / (k + 1)
        @test integrate(x .^ k, D) == (xmax^(k + 1) - xmin^(k + 1)) / (k + 1)
    end
    k = degree + 1
    @test sum(M * (x .^ k)) != (xmax^(k + 1) - xmin^(k + 1)) / (k + 1)
    return nothing
end

# Unit vectors picking out the boundary nodes
function boundary_vectors(D)
    x = grid(D)
    eL = zeros(eltype(x), length(x))
    eL[begin] = 1
    eR = zeros(eltype(x), length(x))
    eR[end] = 1
    return eL, eR
end

@testset "First-derivative operators" begin
    # accuracy orders of the rational first-derivative operators
    sources_and_orders = ((MattssonNordström2004(), (2, 4, 6, 8)),
                          (MattssonSvärdNordström2004(), (2, 4, 6, 8)),
                          (MattssonSvärdShoeybi2008(), (2, 4, 6, 8)),
                          (Mattsson2014(), (2, 4, 6)),
                          (MattssonAlmquistCarpenter2014Extended(), (2, 4, 6)),
                          (DienerDorbandSchnetterTiglio2007(), (2, 4)),
                          (SharanBradyLivescu2022(RT(1 // 2), RT(2 // 3)), (2, 4, 6)))

    @testset "$(nameof(typeof(source)))" for (source, acc_orders) in sources_and_orders
        @testset "accuracy order $acc_order" for acc_order in acc_orders
            D = derivative_operator(source, 1, acc_order, XMIN, XMAX, NNODES)
            @test eltype(grid(D)) == RT
            @test derivative_order(D) == 1
            @test accuracy_order(D) == acc_order

            # Order of accuracy of the interior stencil and of the boundary
            # closure. The boundary closures of the first-derivative operators
            # of Mattsson (2014) are one order less accurate than usual since
            # they share the norm of the third- and fourth-derivative operators
            # derived there.
            boundary = if source isa Mattsson2014 && acc_order > 2
                acc_order ÷ 2 - 1
            else
                acc_order ÷ 2
            end
            test_exactness(D; interior = acc_order, boundary)

            # SBP property M D + Dᵀ M = B
            M = mass_matrix(D)
            A = Matrix(D)
            @test M * A + A' * M == mass_matrix_boundary(D)

            # the mass matrix is a diagonal, positive definite quadrature rule
            @test M isa Diagonal
            @test all(>(0), diag(M))
            test_quadrature_exactness(D; degree = 2 * (acc_order ÷ 2) - 1)

            # boundary functionals
            x = collect(grid(D))
            @test derivative_left(D, x, Val{0}()) == first(x)
            @test derivative_right(D, x, Val{0}()) == last(x)
            @test integrate_boundary(x, D) == last(x) - first(x)
        end
    end
end

@testset "Second-derivative operators" begin
    sources_and_orders = ((MattssonNordström2004(), (2, 4, 6, 8)),
                          (MattssonSvärdNordström2004(), (2, 4, 6, 8)),
                          (MattssonSvärdShoeybi2008(), (2, 4, 6, 8)),
                          (Mattsson2014(), (2, 4, 6)))

    @testset "$(nameof(typeof(source)))" for (source, acc_orders) in sources_and_orders
        @testset "accuracy order $acc_order" for acc_order in acc_orders
            D = derivative_operator(source, 2, acc_order, XMIN, XMAX, NNODES)
            @test derivative_order(D) == 2
            @test accuracy_order(D) == acc_order

            # As for the first derivative, the boundary closures of
            # Mattsson (2014) are one order less accurate than usual.
            boundary = if source isa Mattsson2014 && acc_order > 2
                acc_order ÷ 2
            else
                acc_order ÷ 2 + 1
            end
            test_exactness(D; interior = acc_order + 1, boundary)

            # SBP property M D₂ = -A + eR dRᵀ - eL dLᵀ with a symmetric
            # (negative semidefinite) matrix A
            M = mass_matrix(D)
            A = Matrix(D)
            eL, eR = boundary_vectors(D)
            dL = derivative_left(D, Val{1}())
            dR = derivative_right(D, Val{1}())
            R = M * A - (eR * dR' - eL * dL')
            @test R == R'
            # the skew-symmetric part, which is also checked in the
            # floating point tests
            @test M * A - A' * M == eR * dR' - eL * dL' - dR * eR' + dL * eL'

            # the boundary derivative functionals are one order more accurate
            # than the boundary closure of the first-derivative operators
            x = collect(grid(D))
            for k in 0:(acc_order ÷ 2 + 1)
                @test derivative_left(D, x .^ k, Val{1}()) ==
                      derivative_of_monomial(first(x), k, 1)
                @test derivative_right(D, x .^ k, Val{1}()) ==
                      derivative_of_monomial(last(x), k, 1)
            end
        end
    end
end

@testset "Third- and fourth-derivative operators (Mattsson2014)" begin
    # The boundary closures of the narrow-stencil third- and fourth-derivative
    # operators are of very low order; some of them are not even pointwise
    # consistent, i.e., they do not differentiate x^der_order exactly.
    # Convergence is still obtained since the corresponding boundary weights of
    # the norm are O(Δx). Hence, the degrees of exactness of the boundary
    # closures are listed explicitly here.
    boundary_degrees = Dict((3, 2) => 1, (3, 4) => 3, (3, 6) => 4,
                            (4, 2) => 2, (4, 4) => 3, (4, 6) => 4)

    @testset "derivative order $der_order" for der_order in (3, 4)
        @testset "accuracy order $acc_order" for acc_order in (2, 4, 6)
            D = derivative_operator(Mattsson2014(), der_order, acc_order, XMIN, XMAX,
                                    NNODES)
            @test derivative_order(D) == der_order
            @test accuracy_order(D) == acc_order

            test_exactness(D; interior = acc_order + der_order - 1,
                           boundary = boundary_degrees[(der_order, acc_order)])

            M = mass_matrix(D)
            A = Matrix(D)
            eL, eR = boundary_vectors(D)
            dL1 = derivative_left(D, Val{1}())
            dR1 = derivative_right(D, Val{1}())
            dL2 = derivative_left(D, Val{2}())
            dR2 = derivative_right(D, Val{2}())

            if der_order == 3
                # Integrating by parts twice yields
                #   2 ∫ u u‴ = 2 [u u″] - [(u′)²]
                # so that the SBP property of D₃ reads
                #   M D₃ + D₃ᵀ M = eR dR2ᵀ + dR2 eRᵀ - eL dL2ᵀ - dL2 eLᵀ
                #                  - (dR1 dR1ᵀ - dL1 dL1ᵀ)
                boundary_terms = eR * dR2' + dR2 * eR' - eL * dL2' - dL2 * eL' -
                                 (dR1 * dR1' - dL1 * dL1')
                if acc_order == 2
                    # The second-order accurate third-derivative operator does
                    # not satisfy the SBP property. The symmetric part of M D₃
                    # contains an additional contribution coupling dL2 (dR2) to
                    # the fourth (fourth to last) node. It would vanish if the
                    # left boundary block had a fourth row
                    #   (1//16, -5//8, 17//16, 0, -1, 1//2)
                    # instead of the interior stencil, see
                    # https://github.com/ranocha/SummationByPartsOperators.jl/issues/210
                    @test_broken M * A + A' * M == boundary_terms
                else
                    @test M * A + A' * M == boundary_terms
                end
            else
                # Integrating by parts twice yields
                #   ∫ u u⁗ = [u u‴] - [u′ u″] + ∫ u″ u″
                # so that M D₄ minus the boundary terms is symmetric (and the
                # remaining part is positive semidefinite)
                dL3 = derivative_left(D, Val{3}())
                dR3 = derivative_right(D, Val{3}())
                R = M * A - ((eR * dR3' - eL * dL3') - (dR1 * dR2' - dL1 * dL2'))
                @test R == R'
            end
        end
    end
end

@testset "Upwind operators" begin
    sources_and_orders = ((Mattsson2017, (2, 3, 4, 5, 6, 7, 8, 9)),
                          (WilliamsDuru2024, (4, 5, 6, 7)))

    @testset "$source_type" for (source_type, acc_orders) in sources_and_orders
        @testset "accuracy order $acc_order" for acc_order in acc_orders
            D = upwind_operators(source_type; derivative_order = 1,
                                 accuracy_order = acc_order,
                                 xmin = XMIN, xmax = XMAX, N = NNODES)

            # The central operator of an upwind pair of odd order is one order
            # more accurate in the interior.
            interior_central = isodd(acc_order) ? acc_order + 1 : acc_order
            test_exactness(D.minus; interior = acc_order, boundary = acc_order ÷ 2)
            test_exactness(D.plus; interior = acc_order, boundary = acc_order ÷ 2)
            test_exactness(D.central; interior = interior_central,
                           boundary = acc_order ÷ 2)

            Am = Matrix(D.minus)
            Ac = Matrix(D.central)
            Ap = Matrix(D.plus)
            @test Ac == (Am + Ap) / 2

            # all three operators share the same norm
            M = mass_matrix(D.minus)
            @test mass_matrix(D.central) == M
            @test mass_matrix(D.plus) == M

            # SBP property M D₊ + D₋ᵀ M = B
            @test M * Ap + Am' * M == mass_matrix_boundary(D.minus)
            # the central operator is a classical SBP operator
            @test M * Ac + Ac' * M == mass_matrix_boundary(D.minus)
            # M (D₊ - D₋) is symmetric (and negative semidefinite)
            @test M * (Ap - Am) == (M * (Ap - Am))'

            test_quadrature_exactness(D.minus; degree = 2 * (acc_order ÷ 2) - 1)
        end
    end
end

@testset "Periodic operators" begin
    @testset "derivative order $der_order" for der_order in 1:4
        @testset "accuracy order $acc_order" for acc_order in der_order:8
            D = periodic_derivative_operator(der_order, acc_order, XMIN, XMAX, NNODES)
            @test derivative_order(D) == der_order
            @test accuracy_order(D) == acc_order

            # The stencil uses `acc_order + 1` nodes. Symmetric stencils of
            # even-order derivatives are exact for one additional degree.
            degree = acc_order + (iseven(der_order) && iseven(acc_order) ? 1 : 0)
            test_periodic_exactness(D; degree)

            A = Matrix(D)
            if iseven(der_order)
                @test A == A'
            elseif iseven(acc_order)
                @test A == -A'
                # SBP property M D + Dᵀ M = 0
                M = mass_matrix(D)
                @test iszero(M * A + A' * M)
            end
        end
    end

    @testset "upwind, accuracy order $acc_order" for acc_order in 1:8
        D = upwind_operators(periodic_derivative_operator; derivative_order = 1,
                             accuracy_order = acc_order,
                             xmin = XMIN, xmax = XMAX, N = NNODES)
        interior_central = isodd(acc_order) ? acc_order + 1 : acc_order
        test_periodic_exactness(D.minus; degree = acc_order)
        test_periodic_exactness(D.plus; degree = acc_order)
        test_periodic_exactness(D.central; degree = interior_central)

        Am = Matrix(D.minus)
        Ac = Matrix(D.central)
        Ap = Matrix(D.plus)
        @test Ac == (Am + Ap) / 2
        M = mass_matrix(D.minus)
        @test iszero(M * Ap + Am' * M)
        @test M * (Ap - Am) == (M * (Ap - Am))'
    end
end

@testset "Variable coefficient operators (Mattsson2012)" begin
    # The coefficients of the sixth-order operators are only given as truncated
    # decimals, so only the second- and fourth-order operators can be checked.
    @testset "accuracy order $acc_order" for acc_order in (2, 4)
        bfunc = x -> 1 + x^2
        D = var_coef_derivative_operator(Mattsson2012(), 2, acc_order, XMIN, XMAX, NNODES,
                                         bfunc)
        @test derivative_order(D) == 2
        @test accuracy_order(D) == acc_order

        A = Matrix(D)
        M = mass_matrix(D)
        x = collect(grid(D))
        b = bfunc.(x)

        # SBP property M D₂(b) = -A + eR b(xmax) dRᵀ - eL b(xmin) dLᵀ with a
        # symmetric (negative semidefinite) matrix A. The boundary derivative
        # functionals are the ones of the constant coefficient operator.
        D2 = derivative_operator(Mattsson2012(), 2, acc_order, XMIN, XMAX, NNODES)
        eL, eR = boundary_vectors(D2)
        dL = derivative_left(D2, Val{1}())
        dR = derivative_right(D2, Val{1}())
        R = M * A - (last(b) * eR * dR' - first(b) * eL * dL')
        @test R == R'

        # For a constant coefficient b ≡ 1, the operator is a wide-stencil
        # approximation of the second derivative.
        D1 = var_coef_derivative_operator(Mattsson2012(), 2, acc_order, XMIN, XMAX, NNODES,
                                          one)
        degrees = exactness_degrees(Matrix(D1), x, 2, acc_order + 2)
        nboundary = 2 * acc_order
        @test all(==(acc_order + 1), degrees[(nboundary + 1):(end - nboundary)])
        @test minimum(degrees[1:nboundary]) == acc_order ÷ 2 + 1
    end
end

@testset "Dissipation operators" begin
    @testset "accuracy order $acc_order" for acc_order in (2, 4, 6, 8)
        D = derivative_operator(MattssonSvärdNordström2004(), 1, acc_order, XMIN, XMAX,
                                NNODES)
        M = mass_matrix(D)
        x = collect(grid(D))
        # dissipation operators are implemented up to order 8
        orders = filter(<=(8), (acc_order, acc_order + 2))
        @testset "dissipation order $order" for order in orders
            Di = dissipation_operator(MattssonSvärdNordström2004(), D; order)
            A = Matrix(Di)
            # M Diss is symmetric (and negative semidefinite)
            @test M * A == (M * A)'
            # The dissipation operator annihilates polynomials of degree
            # less than half its order.
            for k in 0:(order ÷ 2 - 1)
                @test iszero(A * (x .^ k))
            end
            @test !iszero(A * (x .^ (order ÷ 2)))
        end
    end
end

@testset "Coupled operators" begin
    @testset "accuracy order $acc_order" for acc_order in (2, 4)
        D = derivative_operator(MattssonNordström2004(), 1, acc_order, RT(0), RT(1), 9)
        @testset "$(nameof(typeof(mesh)))" for mesh in (UniformMesh1D(RT(0), RT(4), 4),
                                                        UniformPeriodicMesh1D(RT(0), RT(4),
                                                                              4))
            @testset "$coupling" for coupling in (couple_continuously,
                                                  couple_discontinuously)
                Dc = coupling(D, mesh)
                A = Matrix(Dc)
                M = mass_matrix(Dc)
                @test M * A + A' * M == mass_matrix_boundary(Dc)
            end
        end
    end
end

@testset "Rational{Int} arithmetic" begin
    # `Rational{Int}` is sufficient for the simplest operators; it overflows for
    # most of the higher-order ones, which is why `Rational{BigInt}` is used
    # above.
    D = derivative_operator(MattssonNordström2004(), 1, 2, -1 // 2, 2 // 1, 11)
    @test eltype(grid(D)) == Rational{Int}
    M = mass_matrix(D)
    A = Matrix(D)
    @test M * A + A' * M == mass_matrix_boundary(D)
    x = collect(grid(D))
    @test A * x == ones(Rational{Int}, length(x))

    Dp = periodic_derivative_operator(1, 2, -1 // 2, 2 // 1, 12)
    @test eltype(grid(Dp)) == Rational{Int}
    Ap = Matrix(Dp)
    @test Ap == -Ap'
end
