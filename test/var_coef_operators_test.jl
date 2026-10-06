using Test
using LinearAlgebra, SparseArrays
using SummationByPartsOperators

# Sources of coefficients for variable coefficient second-derivative operators
# together with the orders of accuracy they provide.
test_list = ((Mattsson2012(), (2, 4, 6)),
             (StiernströmAlmquistMattsson2023(), (4, 6, 8, 10, 12)))

@testset "Test consistency with constant coefficient operators" begin
    for (source, acc_orders) in test_list, acc_order in acc_orders,
        T in (Float32, Float64)

        xmin = -one(T)
        xmax = 2 * one(T)
        N = 101

        D2 = derivative_operator(source, 2, acc_order, xmin, xmax, N)
        D2var = try
            var_coef_derivative_operator(source, 2, acc_order, xmin, xmax, N, one)
        catch err
            !isa(err, ArgumentError) && throw(err)
            nothing
        end
        D2var === nothing && continue

        for compact in (true, false)
            show(IOContext(devnull, :compact => false), D2var)
            show(IOContext(devnull, :compact => false), D2var.coefficients)
        end

        @test maximum(abs, Matrix(D2) - Matrix(D2var)) < 10000 * eps(T)
    end
end

@testset "Test consistency for vanishing coefficients" begin
    for (source, acc_orders) in test_list, T in (Float32, Float64),
        acc_order in acc_orders

        xmin = zero(T)
        xmax = 5 * one(T)
        N = 51

        D2var = var_coef_derivative_operator(source, 2, acc_order, xmin, xmax, N, zero)
        x = grid(D2var)
        u = rand(T, length(x))
        @test maximum(abs, D2var * u) < 10 * eps(T)
    end
end

@testset "Compare mul! with β=0 and mul! without β" begin
    for (source, acc_orders) in test_list, T in (Float32, Float64),
        acc_order in acc_orders

        xmin = zero(T)
        xmax = 5 * one(T)
        N = 51

        D2var_serial = var_coef_derivative_operator(source, 2, acc_order, xmin, xmax, N,
                                                    one)
        D2var_threads = var_coef_derivative_operator(source, 2, acc_order, xmin, xmax, N,
                                                     one, ThreadedMode())
        D2var_safe = var_coef_derivative_operator(source, 2, acc_order, xmin, xmax, N, one,
                                                  SafeMode())
        D2var_full = Matrix(D2var_serial)
        D2var_sparse = sparse(D2var_serial)

        x = grid(D2var_serial)
        u = x .^ 5
        dest1 = fill(zero(eltype(u)), length(u))
        dest2 = fill(zero(eltype(u)), length(u))

        # The different execution modes and matrix representations sum the same
        # coefficients in a different order. Near the zeros of the result, the
        # round-off errors of the sums dominate, so a relative comparison alone
        # is not meaningful; `atol` bounds the round-off error of a single row.
        rtol = sqrt(eps(T))
        atol = 10 * eps(T) * maximum(sum(abs, D2var_full, dims = 2)) * maximum(abs, u)
        approx_equal(a, b) = all(i -> isapprox(a[i], b[i]; atol, rtol), eachindex(a))

        mul!(dest1, D2var_serial, u, one(T), zero(T))
        mul!(dest2, D2var_serial, u, one(T))
        @test approx_equal(dest1, dest2)
        mul!(dest1, D2var_safe, u, one(T), zero(T))
        mul!(dest2, D2var_safe, u, one(T))
        @test approx_equal(dest1, dest2)
        mul!(dest1, D2var_threads, u, one(T), zero(T))
        mul!(dest2, D2var_threads, u, one(T))
        @test approx_equal(dest1, dest2)
        mul!(dest2, D2var_serial, u)
        @test approx_equal(dest1, dest2)
        mul!(dest2, D2var_full, u)
        @test approx_equal(dest1, dest2)
        mul!(dest2, D2var_sparse, u)
        @test approx_equal(dest1, dest2)
        mul!(dest2, D2var_safe, u)
        @test approx_equal(dest1, dest2)

        # test allocations
        allocs = @allocated mul!(dest1, D2var_serial, u)
        @test iszero(allocs)
        allocs = @allocated mul!(dest1, D2var_serial, u, one(T))
        @test iszero(allocs)
        allocs = @allocated mul!(dest1, D2var_serial, u, one(T), zero(T))
        @test iszero(allocs)
    end
end

@testset "Test interior symmetry" begin
    for (source, acc_orders) in test_list, T in (Float32, Float64),
        acc_order in acc_orders

        xmin = one(T)
        xmax = 2 * one(T)
        N = 51

        D2var = var_coef_derivative_operator(source, 2, acc_order, xmin, xmax, N, inv)
        M = mass_matrix(D2var)
        D2 = sparse(D2var)
        A = (M * D2)[2:(end - 1), 2:(end - 1)]
        # The entries of `A` grow like `1 / Δx`, so the tolerance is scaled by
        # the size of the operator instead of being a fixed absolute value.
        @test maximum(abs, A - A') < 10 * eps(T) * maximum(abs, A)
    end
end

@testset "Setting coefficients manually" begin
    D0 = var_coef_derivative_operator(Mattsson2012(), 2, 2, -1.0, 1.0, 11, zero)
    @test iszero(@inferred Matrix(D0))

    D0.b .= grid(D0) .^ 2
    @test !iszero(@inferred Matrix(D0))

    D2 = var_coef_derivative_operator(Mattsson2012(), 2, 2, -1.0, 1.0, 11, abs2)
    @test (@inferred Matrix(D2)) ≈ (@inferred Matrix(D0))
end

# The variable coefficient second-derivative operators satisfy the SBP property
#   H * D2(b) = -M(b) + B̄ * S,
# where `M(b)` is symmetric and positive semidefinite, `S` approximates the
# first derivative at the boundaries, and `B̄ = diag(-b[1], 0, ..., 0, b[end])`.
# This is the basis of energy stability/conservation proofs using these
# operators, see also
# https://github.com/ranocha/SummationByPartsOperators.jl/issues/344
@testset "SBP property of M(b)" begin
    for (source, acc_orders) in test_list, acc_order in acc_orders
        T = Float64
        xmin = zero(T)
        xmax = one(T)
        N = 40

        # `S` is the same as for the constant coefficient operators
        D2 = derivative_operator(source, 2, acc_order, xmin, xmax, N)
        S_left = derivative_left(D2, Val(1))
        S_right = derivative_right(D2, Val(1))

        D2var = var_coef_derivative_operator(source, 2, acc_order, xmin, xmax, N, one)
        H = Matrix(mass_matrix(D2var))

        for b in (one.(grid(D2var)), 1 .+ grid(D2var) .^ 2,
                  1 .+ sinpi.(3 .* grid(D2var)) ./ 2)
            D2var.b .= b
            M = -H * Matrix(D2var)
            M[begin, :] .-= b[begin] * S_left
            M[end, :] .+= b[end] * S_right

            # `M(b)` is symmetric
            tol = 1000 * eps(T) * max(one(T), maximum(abs, M))
            @test maximum(abs, M - M') < tol

            # `M(b)` is positive semidefinite with constants in its kernel
            λ = eigvals(Symmetric((M + M') / 2))
            @test minimum(λ) > -tol
            @test abs(λ[1]) < tol
            @test λ[2] > 1 // 10
        end
    end
end

# The boundary-optimized operators of Stiernström, Almquist, Mattsson (2023)
# are fully compatible with the boundary-optimized first-derivative operators
# of Mattsson, Almquist, van der Weide (2018), i.e., they share their
# non-uniform grid, their norm matrix `H`, and their first derivative `D₁`, and
# satisfy `D₂(b) = D₁ * Diagonal(b) * D₁ - H⁻¹ * R(b)`.
@testset "StiernströmAlmquistMattsson2023" begin
    source = StiernströmAlmquistMattsson2023()
    reference = MattssonAlmquistVanDerWeide2018Accurate()

    @testset "not implemented" begin
        @test_throws ArgumentError derivative_operator(source, 2, 5, 0.0, 1.0, 40)
        @test_throws ArgumentError var_coef_derivative_operator(source, 2, 5, 0.0, 1.0, 40,
                                                                one)
        @test_throws ArgumentError var_coef_derivative_operator(source, 3, 4, 0.0, 1.0, 40,
                                                                one)
    end

    @testset "accuracy order $acc_order" for acc_order in (4, 6, 8, 10, 12)
        p = acc_order ÷ 2
        T = Float64
        xmin = -one(T)
        xmax = 2 * one(T)
        N = 40

        # The operators share grid, norm, and first derivative with the
        # boundary-optimized operators they are fully compatible with.
        D1 = derivative_operator(source, 1, acc_order, xmin, xmax, N)
        D1_reference = derivative_operator(reference, 1, acc_order, xmin, xmax, N)
        @test grid(D1) == grid(D1_reference)
        @test Matrix(D1) == Matrix(D1_reference)
        @test mass_matrix(D1) == mass_matrix(D1_reference)

        D2 = var_coef_derivative_operator(source, 2, acc_order, xmin, xmax, N, one)
        @test derivative_order(D2) == 2
        @test accuracy_order(D2) == acc_order
        @test source_of_coefficients(D2) === source
        @test grid(D2) == grid(D1_reference)
        @test mass_matrix(D2) == mass_matrix(D1_reference)

        # The constant coefficient operator is `D₂(b)` for `b ≡ 1`.
        D2_constant = derivative_operator(source, 2, acc_order, xmin, xmax, N)
        @test Matrix(D2_constant) == Matrix(D2)
        @test derivative_order(D2_constant) == 2
        @test accuracy_order(D2_constant) == acc_order
        @test issymmetric(D2_constant) == false
        @test isperiodic(D2_constant) == false
        @test mass_matrix(D2_constant) == mass_matrix(D2)

        # SBP property of the constant coefficient operator
        eL = zeros(T, N)
        eL[begin] = 1
        eR = zeros(T, N)
        eR[end] = 1
        dL = derivative_left(D2_constant, Val{1}())
        dR = derivative_right(D2_constant, Val{1}())
        H = Matrix(mass_matrix(D2_constant))
        A_constant = Matrix(D2_constant)
        @test H * A_constant - A_constant' * H ≈
              eR * dR' - eL * dL' - dR * eR' + dL * eL'

        # The boundary derivatives are exact for linear functions.
        @test derivative_left(D2_constant, grid(D2_constant), Val{1}()) ≈ one(T)
        @test derivative_right(D2_constant, grid(D2_constant), Val{1}()) ≈ one(T)

        # scaling by the mass matrix is invertible
        u = sinpi.(grid(D2_constant))
        v = copy(u)
        scale_by_mass_matrix!(v, D2_constant)
        @test v ≈ H * u
        scale_by_inverse_mass_matrix!(v, D2_constant)
        @test v ≈ u

        cache = D2.coefficients.coefficient_cache
        nleft = SummationByPartsOperators.left_length(cache)
        nright = SummationByPartsOperators.right_length(cache)
        @test nleft == 3 * p
        @test nright == 3 * p

        # Not enough nodes to keep the two boundary closures apart
        @test_throws DimensionMismatch var_coef_derivative_operator(source, 2, acc_order,
                                                                    xmin, xmax, 6 * p - 1,
                                                                    one)

        x = grid(D2)
        b = 1 .+ sinpi.(3 .* x) ./ 2
        D2.b .= b
        A = Matrix(D2)

        # The operator is a narrow-stencil operator: the interior stencils have
        # the minimal width `2p + 1`. In the boundary closure, the first `2p`
        # rows use `3p` points and the following `p` rows use one point more
        # each, cf. the reference implementation.
        for i in (nleft + 1):(N - nright)
            @test all(j -> iszero(A[i, j]), filter(j -> abs(i - j) > p, 1:N))
        end
        for i in 1:nleft
            jmax = i <= 2 * p ? 3 * p : p + i
            @test all(j -> iszero(A[i, j]), (jmax + 1):N)
            @test !iszero(A[i, jmax])
        end
        @test SummationByPartsOperators.lower_bandwidth(D2) == 3 * p - 1
        @test SummationByPartsOperators.upper_bandwidth(D2) == 3 * p - 1

        # The grid and the operators are symmetric with respect to the center of
        # the domain, i.e., reversing the variable coefficients mirrors `D₂(b)`.
        D2_mirrored = var_coef_derivative_operator(source, 2, acc_order, xmin, xmax, N,
                                                   one)
        D2_mirrored.b .= reverse(b)
        @test maximum(abs, Matrix(D2_mirrored) - reverse(A)) <
              10 * eps(T) * maximum(abs, A)

        # Order of accuracy. For `b ≡ 1`, the interior stencils are exact for
        # polynomials up to degree `2p + 1` and the boundary closures up to
        # degree `p`, i.e., the boundary closures are `p - 1`-th order accurate.
        # The tolerances account for the coefficients of `D₁` and of the grid
        # being available only as truncated decimals.
        D2.b .= one.(x)
        A = Matrix(D2)
        scale = maximum(abs, A)
        residual = function (k)
            u = x .^ k
            e = A * u - (k < 2 ? zero(x) : (k * (k - 1)) .* x .^ (k - 2))
            return (maximum(abs, e[1:nleft]) / (scale * maximum(abs, u)),
                    maximum(abs, e[(nleft + 1):(end - nright)]) / (scale * maximum(abs, u)))
        end
        for k in 0:(acc_order + 1)
            boundary, interior = residual(k)
            @test interior < 1e-15
            if k <= p
                @test boundary < 1e-12
            end
        end
        # The boundary closures are not exact for one degree more.
        @test first(residual(p + 1)) > 1000 * first(residual(p))

        # For a variable coefficient `b` of degree two, the degrees of exactness
        # are reduced accordingly.
        D2.b .= 1 .+ x .^ 2
        A = Matrix(D2)
        scale = maximum(abs, A)
        residual = function (k)
            u = x .^ k
            du = k < 1 ? zero(x) : k .* x .^ (k - 1)
            d2u = k < 2 ? zero(x) : (k * (k - 1)) .* x .^ (k - 2)
            e = A * u - (2 .* x .* du + (1 .+ x .^ 2) .* d2u)
            return (maximum(abs, e[1:nleft]) / (scale * maximum(abs, u)),
                    maximum(abs, e[(nleft + 1):(end - nright)]) / (scale * maximum(abs, u)))
        end
        for k in 0:(acc_order - 1)
            boundary, interior = residual(k)
            @test interior < 1e-15
            if k <= p - 1
                @test boundary < 1e-12
            end
        end
        # The separation between the exact and the inexact degrees is smaller
        # here than for `b ≡ 1` since the coefficients of `D₁` and the grid
        # points are only available as decimals truncated to 14 digits.
        @test first(residual(p)) > 20 * first(residual(p - 1))
        @test last(residual(acc_order + 1)) > 1000 * last(residual(acc_order - 1))
    end

    # Convergence of `D₂(b) u` to `(b u′)′` in the discrete `H` norm. The
    # expected rate is `p - 1/2`, resulting from the `p - 1`-th order accurate
    # boundary closures and the `1/2` gained by the `H` norm weighting a fixed
    # number of boundary nodes with `Δx`.
    @testset "convergence order $acc_order" for acc_order in (4, 6, 8, 10)
        p = acc_order ÷ 2
        bfunc(x) = 1 + x^2 / 2 + sin(x) / 4
        dbfunc(x) = x + cos(x) / 4
        ufunc(x) = exp(sin(3x))
        dufunc(x) = 3cos(3x) * ufunc(x)
        d2ufunc(x) = (9cos(3x)^2 - 9sin(3x)) * ufunc(x)

        errors = map((50, 100)) do N
            D = var_coef_derivative_operator(StiernströmAlmquistMattsson2023(), 2,
                                             acc_order, -1.0, 1.0, N, bfunc)
            x = grid(D)
            e = D * ufunc.(x) - (dbfunc.(x) .* dufunc.(x) + bfunc.(x) .* d2ufunc.(x))
            return sqrt(sum(diag(mass_matrix(D)) .* e .^ 2))
        end
        @test log2(errors[1] / errors[2]) > p - 1.5
    end
end

# `integrate` uses the quadrature rule given by the diagonal norm of the
# operator. For an operator with interior order of accuracy `acc_order`, this
# quadrature is exact for polynomials of degree `acc_order - 1`.
@testset "integrate" begin
    for (source, acc_orders) in test_list, T in (Float32, Float64),
        acc_order in acc_orders

        xmin = -one(T)
        xmax = 2 * one(T)
        N = 41

        D = var_coef_derivative_operator(source, 2, acc_order, xmin, xmax, N, one)
        x = grid(D)
        M = mass_matrix(D)

        u = sinpi.(x)
        @test integrate(u, D) == integrate(identity, u, D)
        @test integrate(u, D) ≈ sum(diag(M) .* u)
        @test integrate(abs2, u, D) ≈ sum(diag(M) .* u .^ 2)

        for k in 0:(acc_order - 1)
            @test integrate(x .^ k, D) ≈ (xmax^(k + 1) - xmin^(k + 1)) / (k + 1)
        end

        @test_throws DimensionMismatch integrate(u[(begin + 1):end], D)
    end
end
