# Tests of the order of accuracy and of the SBP (symmetry) properties of the
# operators using exact rational arithmetic. Since all relations checked here
# hold exactly for the exact coefficients, they can be checked with `==`
# instead of `≈`, which makes these tests much stronger than the corresponding
# floating point tests.
#
# `Rational{Int}` overflows quickly for the more complicated operators, so
# `Rational{BigInt}` is used throughout. It is slower but always works.
#
# The order of accuracy cannot be checked for operators whose coefficients are
# only available as (truncated) floating point numbers. These are
# - `DienerDorbandSchnetterTiglio2007` with accuracy orders 6 and 8
# - `MattssonAlmquistCarpenter2014Optimal` (both the grid points and the
#   coefficients are truncated decimals)
# - `MattssonAlmquistVanDerWeide2018Minimal`, `MattssonAlmquistVanDerWeide2018Accurate`
# - `MattssonNiemeläWinters2026`
# - `Mattsson2012` variable coefficient operators with accuracy order 6
# - `LegendreDerivativeOperator`, `FourierDerivativeOperator`, and
#   `function_space_operator`s
# The SBP property is structural and can still hold exactly for such operators;
# it is checked whenever it does, see the testset
# "Operators with floating point coefficients" and the sixth-order case of the
# testset "Variable coefficient operators (Mattsson2012)" below.

module RationalArithmeticTest

using Test
using LinearAlgebra
using SummationByPartsOperators

const RT = Rational{BigInt}

# Domain and number of nodes used for all tests below. On uniform grids, this
# yields the grid spacing `Δx = 5 // 64`, which is deliberately not equal to
# unity to test the scaling of the coefficients by powers of `Δx` as well.
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
    exactness_degrees(A, x, exact, max_degree)

For every row of the matrix `A` representing a differential operator on the
grid `x`, compute the maximal degree `p ≤ max_degree` such that all monomials
`x^k` with `k ≤ p` are handled exactly. Here, `exact(xᵢ, k)` must return the
exact value of the differential operator applied to `x^k` at the node `xᵢ`.
A value of `-1` means that not even constants are handled exactly.
"""
function exactness_degrees(A::AbstractMatrix, x::AbstractVector, exact,
                           max_degree::Integer)
    n = size(A, 1)
    degrees = fill(-1, n)
    active = trues(n)
    for k in 0:max_degree
        residual = A * (x .^ k)
        for i in 1:n
            active[i] || continue
            if residual[i] == exact(x[i], k)
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
monomials up to degree `interior` (and no more) exactly and that the least
accurate rows of the boundary closures do so up to degree `boundary`
(and no more).
"""
function test_exactness(D; interior, boundary)
    x = collect(grid(D))
    der_order = derivative_order(D)
    exact = (xi, k) -> derivative_of_monomial(xi, k, der_order)
    degrees = exactness_degrees(Matrix(D), x, exact, interior + 1)
    nleft = length(D.coefficients.left_boundary)
    nright = length(D.coefficients.right_boundary)
    # make sure there are interior nodes left to check
    @test nleft + nright < length(x)
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
    der_order = derivative_order(D)
    exact = (xi, k) -> derivative_of_monomial(xi, k, der_order)
    degrees = exactness_degrees(Matrix(D), x, exact, degree + 1)
    nlower = length(D.coefficients.lower_coef)
    nupper = length(D.coefficients.upper_coef)
    # make sure there are nodes left whose stencil does not wrap around
    @test nlower + nupper < length(x)
    @test all(==(degree), degrees[(nlower + 1):(end - nupper)])
    return nothing
end

"""
    test_quadrature_exactness(D; degree)

Check that the quadrature rule given by the mass matrix of `D` is exact for
polynomials up to degree `degree` (and no more) on the domain of `D`. Note
that the boundaries of the domain are not necessarily the first and last node
of the grid, e.g., for operators coupled continuously on a periodic mesh.
"""
function test_quadrature_exactness(D; degree)
    xmin = SummationByPartsOperators.xmin(D)
    xmax = SummationByPartsOperators.xmax(D)
    x = collect(grid(D))
    M = mass_matrix(D)
    for k in 0:degree
        @test sum(M * (x .^ k)) == (xmax^(k + 1) - xmin^(k + 1)) / (k + 1)
        @test integrate(x .^ k, D) == (xmax^(k + 1) - xmin^(k + 1)) / (k + 1)
    end
    k = degree + 1
    @test sum(M * (x .^ k)) != (xmax^(k + 1) - xmin^(k + 1)) / (k + 1)
    @test integrate(x .^ k, D) != (xmax^(k + 1) - xmin^(k + 1)) / (k + 1)
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
            # closure. The boundary closures of the operators of Mattsson
            # (2014) are one order less accurate than usual: they share the
            # norm of the third- and fourth-derivative operators derived
            # there, and the accuracy conditions combined with the SBP
            # property have no solution of the usual order for that norm.
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
            # Note that this does not apply to the boundary derivative
            # functionals checked below.
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

            # the mass matrix is a diagonal, positive definite quadrature rule
            @test M isa Diagonal
            @test all(>(0), diag(M))
            test_quadrature_exactness(D; degree = 2 * (acc_order ÷ 2) - 1)

            # The boundary derivative functionals are exact for polynomials of
            # degree `acc_order ÷ 2 + 1`, i.e., one more than the usual
            # boundary closure of the first-derivative operators.
            x = collect(grid(D))
            for k in 0:(acc_order ÷ 2 + 1)
                @test derivative_left(D, x .^ k, Val{1}()) ==
                      derivative_of_monomial(first(x), k, 1)
                @test derivative_right(D, x .^ k, Val{1}()) ==
                      derivative_of_monomial(last(x), k, 1)
            end
            k = acc_order ÷ 2 + 2
            @test derivative_left(D, x .^ k, Val{1}()) !=
                  derivative_of_monomial(first(x), k, 1)
            @test derivative_right(D, x .^ k, Val{1}()) !=
                  derivative_of_monomial(last(x), k, 1)

            # boundary functionals
            @test derivative_left(D, x, Val{0}()) == first(x)
            @test derivative_right(D, x, Val{0}()) == last(x)
            @test integrate_boundary(x, D) == last(x) - first(x)
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

            # These operators share the norm of the first-derivative operators
            # of Mattsson (2014), which is a diagonal, positive definite
            # quadrature rule.
            @test mass_matrix(D) ==
                  mass_matrix(derivative_operator(Mattsson2014(), 1, acc_order, XMIN, XMAX,
                                                  NNODES))
            @test mass_matrix(D) isa Diagonal
            @test all(>(0), diag(mass_matrix(D)))
            test_quadrature_exactness(D; degree = 2 * (acc_order ÷ 2) - 1)

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
                @test M * A + A' * M == boundary_terms

                # Mattsson (2014) writes
                #   D₃ = M⁻¹ (R + dL1 dL1ᵀ / 2 - dR1 dR1ᵀ / 2
                #             - eL dL2ᵀ + eR dR2ᵀ)
                # with an antisymmetric matrix R, which is equivalent to the
                # SBP property checked above.
                R = M * A - (dL1 * dL1' / 2 - dR1 * dR1' / 2 - eL * dL2' + eR * dR2')
                @test R == -R'
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

            # all three operators share the same norm, which is a diagonal,
            # positive definite quadrature rule
            M = mass_matrix(D.minus)
            @test mass_matrix(D.central) == M
            @test mass_matrix(D.plus) == M
            @test M isa Diagonal
            @test all(>(0), diag(M))

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

    # The coefficients used by `periodic_central_derivative_operator` are
    # evaluated exactly following Beljadid, LeFloch, Mishra, Parés (2017). They
    # agree with the central stencils computed via the algorithm of Fornberg
    # (1998) used by `periodic_derivative_operator`.
    @testset "central stencils, derivative order $der_order" for der_order in 1:3
        @testset "accuracy order $acc_order" for acc_order in (2, 4, 6, 8)
            Dc = periodic_central_derivative_operator(der_order, acc_order, XMIN, XMAX,
                                                      NNODES)
            @test eltype(grid(Dc)) == RT
            @test source_of_coefficients(Dc) isa BeljaddLeFlochMishraParés2017
            @test Matrix(Dc) ==
                  Matrix(periodic_derivative_operator(der_order, acc_order, XMIN, XMAX,
                                                      NNODES))
        end
    end

    # Stencils that are not centered around the node where the derivative is
    # approximated
    @testset "non-central stencils, derivative order $der_order" for der_order in 1:2
        @testset "accuracy order $acc_order" for acc_order in der_order:6
            @testset "left offset $left_offset" for left_offset in (-acc_order):0
                D = periodic_derivative_operator(der_order, acc_order, XMIN, XMAX, NNODES,
                                                 left_offset)
                @test accuracy_order(D) == acc_order
                # The stencil uses `acc_order + 1` nodes. It is symmetric if
                # `left_offset` is `-acc_order / 2`; symmetric stencils of
                # even-order derivatives are exact for one additional degree.
                symmetric = 2 * left_offset == -acc_order
                degree = acc_order + (iseven(der_order) && symmetric ? 1 : 0)
                test_periodic_exactness(D; degree)
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

        # all three operators share the same norm
        M = mass_matrix(D.minus)
        @test mass_matrix(D.central) == M
        @test mass_matrix(D.plus) == M

        @test iszero(M * Ap + Am' * M)
        @test M * (Ap - Am) == (M * (Ap - Am))'
    end
end

@testset "Periodic operators with wide stencils" begin
    # The smooth noise-robust differentiators of Holoborodko (2008) and the
    # (super) Lanczos low-noise differentiators use wider stencils than needed
    # for their order of accuracy to reduce the amplification of noise. Their
    # coefficients are exact rational numbers.
    @testset "Holoborodko2008" begin
        # the implemented stencil widths per derivative and accuracy order
        widths = Dict((1, 2) => (5, 7, 9, 11), (1, 4) => (7, 9, 11),
                      (2, 2) => (5, 7, 9), (2, 4) => (7, 9))
        @testset "derivative order $der_order" for der_order in (1, 2)
            @testset "accuracy order $acc_order" for acc_order in (2, 4)
                @testset "stencil width $stencil_width" for stencil_width in widths[(der_order,
                                                                                     acc_order)]
                    D = periodic_derivative_operator(Holoborodko2008(), der_order,
                                                     acc_order, XMIN, XMAX, NNODES;
                                                     stencil_width)
                    @test eltype(grid(D)) == RT
                    @test derivative_order(D) == der_order
                    @test accuracy_order(D) == acc_order
                    @test length(D.coefficients.lower_coef) == stencil_width ÷ 2

                    # The stencils are symmetric, so those of even-order
                    # derivatives are exact for one additional degree.
                    degree = acc_order + (iseven(der_order) ? 1 : 0)
                    test_periodic_exactness(D; degree)

                    A = Matrix(D)
                    if iseven(der_order)
                        @test A == A'
                    else
                        @test A == -A'
                    end
                end
            end
        end
    end

    @testset "LanczosLowNoise" begin
        # the implemented stencil widths per accuracy order
        widths = Dict(2 => (5, 7, 9, 11), 4 => (7, 9, 11))
        @testset "accuracy order $acc_order" for acc_order in (2, 4)
            @testset "stencil width $stencil_width" for stencil_width in widths[acc_order]
                D = periodic_derivative_operator(LanczosLowNoise(); derivative_order = 1,
                                                 accuracy_order = acc_order, stencil_width,
                                                 xmin = XMIN, xmax = XMAX, N = NNODES)
                @test eltype(grid(D)) == RT
                @test derivative_order(D) == 1
                @test accuracy_order(D) == acc_order
                @test length(D.coefficients.lower_coef) == stencil_width ÷ 2

                test_periodic_exactness(D; degree = acc_order)
                A = Matrix(D)
                @test A == -A'
            end
        end
    end
end

@testset "Periodic dissipation operators" begin
    # The dissipation operator depends on `D` only via its grid; the order of
    # dissipation defaults to the accuracy order of `D`.
    @testset "accuracy order $acc_order" for acc_order in (2, 4, 6, 8)
        D = periodic_derivative_operator(1, acc_order, XMIN, XMAX, NNODES)
        M = mass_matrix(D)
        x = collect(grid(D))
        Δx = x[2] - x[1]
        Di = dissipation_operator(D)
        A = Matrix(Di)

        # The dissipation operator uses undivided differences, i.e., it is the
        # `acc_order`-th derivative operator of accuracy order `acc_order`
        # scaled by Δx^acc_order and by a sign making it negative semidefinite.
        Dref = periodic_derivative_operator(acc_order, acc_order, XMIN, XMAX, NNODES)
        @test A == (-1)^(1 + acc_order ÷ 2) * Δx^acc_order * Matrix(Dref)

        # M Diss is symmetric (and negative semidefinite)
        @test M * A == (M * A)'
    end
end

@testset "Variable coefficient operators (Mattsson2012)" begin
    @testset "accuracy order $acc_order" for acc_order in (2, 4, 6)
        # a variable coefficient that is a polynomial of degree `degree_b`
        bfunc = x -> 1 + x^2
        dbfunc = x -> 2 * x
        degree_b = 2

        D = var_coef_derivative_operator(Mattsson2012(), 2, acc_order, XMIN, XMAX, NNODES,
                                         bfunc)
        @test derivative_order(D) == 2
        @test accuracy_order(D) == acc_order

        A = Matrix(D)
        M = mass_matrix(D)
        x = collect(grid(D))
        b = bfunc.(x)

        # Mattsson (2012) reuses the constant coefficient operators, and thus
        # also the boundary derivative functionals, of Mattsson & Nordström
        # (2004).
        D2 = derivative_operator(Mattsson2012(), 2, acc_order, XMIN, XMAX, NNODES)
        eL, eR = boundary_vectors(D2)
        dL = derivative_left(D2, Val{1}())
        dR = derivative_right(D2, Val{1}())

        # SBP property M D₂(b) = -A + b(xmax) eR dRᵀ - b(xmin) eL dLᵀ with a
        # symmetric (negative semidefinite) matrix A. This is a structural
        # property of the coefficients and holds exactly also for the
        # sixth-order operators.
        R = M * A - (last(b) * eR * dR' - first(b) * eL * dL')
        @test R == R'
        @test mass_matrix(D) == mass_matrix(D2)

        # The coefficients of the sixth-order variable coefficient operators
        # are only given as truncated decimals. Hence, the reduction to the
        # constant coefficient operator and the order of accuracy can only be
        # checked exactly for the second- and fourth-order operators.
        if acc_order != 6
            # For a constant coefficient b ≡ 1, the operator reduces to the
            # constant coefficient second-derivative operator.
            D1 = var_coef_derivative_operator(Mattsson2012(), 2, acc_order, XMIN, XMAX,
                                              NNODES, one)
            @test Matrix(D1) == Matrix(D2)
            @test mass_matrix(D1) == mass_matrix(D2)

            # The operator applied to u is exact whenever b u′ is a polynomial
            # of sufficiently low degree. Hence, the degrees of exactness are
            # reduced by the degree of b compared to the constant coefficient
            # operator.
            function exact(xi, k)
                return dbfunc(xi) * derivative_of_monomial(xi, k, 1) +
                       bfunc(xi) * derivative_of_monomial(xi, k, 2)
            end
            interior = acc_order + 1 - degree_b
            boundary = acc_order ÷ 2 + 1 - degree_b
            degrees = exactness_degrees(A, x, exact, interior + 1)
            cache = D.coefficients.coefficient_cache
            nleft = SummationByPartsOperators.left_length(cache)
            nright = SummationByPartsOperators.right_length(cache)
            @test nleft + nright < length(x)
            @test all(==(interior), degrees[(nleft + 1):(end - nright)])
            @test minimum(degrees[1:nleft]) == boundary
            @test minimum(degrees[(end - nright + 1):end]) == boundary
        end
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
    @testset "accuracy order $acc_order" for acc_order in (2, 4, 6, 8)
        # the operator used on each element of the mesh
        D = derivative_operator(MattssonNordström2004(), 1, acc_order, RT(0), RT(1), 17)
        meshes = (UniformMesh1D(XMIN, XMAX, 4), UniformPeriodicMesh1D(XMIN, XMAX, 4))

        @testset "$(nameof(typeof(mesh)))" for mesh in meshes
            @testset "$coupling" for coupling in (couple_continuously,
                                                  couple_discontinuously)
                Dc = coupling(D, mesh)
                A = Matrix(Dc)
                M = mass_matrix(Dc)

                # SBP property M D + Dᵀ M = B
                @test M * A + A' * M == mass_matrix_boundary(Dc)

                # the mass matrix is a diagonal, positive definite quadrature
                # rule
                @test M isa Diagonal
                @test all(>(0), diag(M))

                # The quadrature rule is exact for polynomials of the same
                # degree as the one of a single element. The continuous
                # coupling on a periodic mesh is the exception: it identifies
                # the first and the last node, so that only constants are
                # integrated exactly.
                degree = if mesh isa UniformPeriodicMesh1D &&
                            coupling === couple_continuously
                    0
                else
                    2 * (acc_order ÷ 2) - 1
                end
                test_quadrature_exactness(Dc; degree)
            end
        end

        # The discontinuous coupling using the upwind numerical fluxes
        # `Val{:minus}()` and `Val{:plus}()` yields an upwind SBP pair sharing
        # the norm of the central coupling.
        @testset "upwind coupling, $(nameof(typeof(mesh)))" for mesh in meshes
            Dminus = couple_discontinuously(D, mesh, Val{:minus}())
            Dcentral = couple_discontinuously(D, mesh, Val{:central}())
            Dplus = couple_discontinuously(D, mesh, Val{:plus}())
            Am = Matrix(Dminus)
            Ac = Matrix(Dcentral)
            Ap = Matrix(Dplus)

            M = mass_matrix(Dminus)
            @test mass_matrix(Dcentral) == M
            @test mass_matrix(Dplus) == M
            @test Ac == (Am + Ap) / 2

            # SBP property M D₊ + D₋ᵀ M = B
            @test M * Ap + Am' * M == mass_matrix_boundary(Dminus)
            # M (D₊ - D₋) is symmetric (and negative semidefinite)
            @test M * (Ap - Am) == (M * (Ap - Am))'
        end
    end
end

@testset "Operators with floating point coefficients" begin
    # The coefficients of these operators are truncated decimals, so their
    # order of accuracy cannot be checked exactly. The SBP property is a
    # structural property of the coefficients as they are stored, though, and
    # still holds exactly for some of them. It does not hold exactly for
    # `DienerDorbandSchnetterTiglio2007` with accuracy orders 6 and 8,
    # `MattssonAlmquistVanDerWeide2018Minimal`,
    # `MattssonAlmquistVanDerWeide2018Accurate`, and
    # `MattssonNiemeläWinters2026`, which are thus not checked here.
    @testset "MattssonAlmquistCarpenter2014Optimal" begin
        @testset "accuracy order $acc_order" for acc_order in (2, 4, 6, 8)
            D = derivative_operator(MattssonAlmquistCarpenter2014Optimal(), 1, acc_order,
                                    XMIN, XMAX, NNODES)
            M = mass_matrix(D)
            A = Matrix(D)
            @test M * A + A' * M == mass_matrix_boundary(D)
            @test M isa Diagonal
            @test all(>(0), diag(M))
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

end # module
