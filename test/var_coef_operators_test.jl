using Test
using LinearAlgebra, SparseArrays
using SummationByPartsOperators

test_list = (Mattsson2012(),)

@testset "Test consistency with constant coefficient operators" begin
    for source in test_list, acc_order in (2, 4, 6), T in (Float32, Float64)
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
    for source in test_list, T in (Float32, Float64), acc_order in (2, 4, 6)
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
    for T in (Float32, Float64), acc_order in (2, 4, 6)
        xmin = zero(T)
        xmax = 5 * one(T)
        N = 51
        source = Mattsson2012()

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

        mul!(dest1, D2var_serial, u, one(T), zero(T))
        mul!(dest2, D2var_serial, u, one(T))
        @test all(i -> dest1[i] ≈ dest2[i], eachindex(u))
        mul!(dest1, D2var_safe, u, one(T), zero(T))
        mul!(dest2, D2var_safe, u, one(T))
        @test all(i -> dest1[i] ≈ dest2[i], eachindex(u))
        mul!(dest1, D2var_threads, u, one(T), zero(T))
        mul!(dest2, D2var_threads, u, one(T))
        @test all(i -> dest1[i] ≈ dest2[i], eachindex(u))
        mul!(dest2, D2var_serial, u)
        @test all(i -> dest1[i] ≈ dest2[i], eachindex(u))
        mul!(dest2, D2var_full, u)
        @test all(i -> dest1[i] ≈ dest2[i], eachindex(u))
        mul!(dest2, D2var_safe, u)
        @test all(i -> dest1[i] ≈ dest2[i], eachindex(u))

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
    for source in test_list, T in (Float32, Float64), acc_order in (2, 4, 6)
        xmin = one(T)
        xmax = 2 * one(T)
        N = 51

        D2var = var_coef_derivative_operator(source, 2, acc_order, xmin, xmax, N, inv)
        M = mass_matrix(D2var)
        D2 = sparse(D2var)
        A = (M * D2)[2:(end - 1), 2:(end - 1)]
        @test maximum(abs, A - A') < 200 * eps(T)
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

# The operators of Mattsson (2012) satisfy the SBP property
#   H * D2(b) = -M(b) + B̄ * S,
# where `M(b)` is symmetric and positive semidefinite, `S` approximates the
# first derivative at the boundaries, and `B̄ = diag(-b[1], 0, ..., 0, b[end])`.
# This is the basis of energy stability/conservation proofs using these
# operators, see also
# https://github.com/ranocha/SummationByPartsOperators.jl/issues/344
@testset "SBP property of M(b)" begin
    for source in test_list, acc_order in (2, 4, 6)
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
            @test maximum(abs, M - M') < 1000 * eps(T)

            # `M(b)` is positive semidefinite with constants in its kernel
            λ = eigvals(Symmetric((M + M') / 2))
            @test minimum(λ) > -1000 * eps(T)
            @test abs(λ[1]) < 1000 * eps(T)
            @test λ[2] > 1 // 10
        end
    end
end
