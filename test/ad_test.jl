module TestAD

using ForwardDiff
using StructArrays
using SummationByPartsOperators

using LinearAlgebra: Diagonal, I, mul!
using Test

@testset "Jacobian" begin
    @testset "periodic_derivative_operator" begin
        D = periodic_derivative_operator(derivative_order = 1, accuracy_order = 2,
                                         xmin = 0.0, xmax = 1.0, N = 8)
        u = rand(size(D, 2))
        f = let D = D
            f(u) = D * u
        end
        J = ForwardDiff.jacobian(f, u)
        @test J ≈ @inferred Matrix(D)
    end

    @testset "derivative_operator" begin
        D = derivative_operator(MattssonNordström2004(),
                                derivative_order = 1, accuracy_order = 2,
                                xmin = 0.0, xmax = 1.0, N = 8)
        u = rand(size(D, 2))
        f = let D = D
            f(u) = D * u
        end
        J = ForwardDiff.jacobian(f, u)
        @test J ≈ @inferred Matrix(D)
    end
end

@testset "mul! with vectors of Duals" begin
    # `FastMode` used to be much slower than `SafeMode` for `Dual`s since
    # `@turbo` falls back to a scalar loop using `Base.FastMath` operations for
    # element types it cannot handle, see
    # https://github.com/ranocha/SummationByPartsOperators.jl/issues/421
    function make_dual(value, partials)
        ForwardDiff.Dual{:ad_test}(value, ForwardDiff.Partials(Tuple(partials)))
    end

    for mode in (FastMode(), SafeMode(), ThreadedMode()),
        accuracy_order in (2, 4), N in (17, 30), num_partials in (1, 2, 5)

        for D in (derivative_operator(MattssonNordström2004();
                                      derivative_order = 1,
                                      accuracy_order = accuracy_order,
                                      xmin = 0.0, xmax = 1.0, N = N, mode = mode),
                  periodic_derivative_operator(derivative_order = 1,
                                               accuracy_order = accuracy_order,
                                               xmin = 0.0, xmax = 1.0, N = N,
                                               mode = mode))
            A = Matrix(D)
            values = randn(N)
            partials = randn(N, num_partials)
            u = [make_dual(values[i], view(partials, i, :)) for i in 1:N]

            du = similar(u)
            mul!(du, D, u)
            @test ForwardDiff.value.(du) ≈ A * values
            for j in 1:num_partials
                @test ForwardDiff.partials.(du, j) ≈ A * partials[:, j]
            end

            mul!(du, D, u, 2.0)
            @test ForwardDiff.value.(du) ≈ 2 * A * values
            for j in 1:num_partials
                @test ForwardDiff.partials.(du, j) ≈ 2 * A * partials[:, j]
            end

            dest = [make_dual(values[i], view(partials, i, :)) for i in 1:N]
            mul!(dest, D, u, 2.0, 3.0)
            @test ForwardDiff.value.(dest) ≈ 2 * A * values + 3 * values
            for j in 1:num_partials
                @test ForwardDiff.partials.(dest, j) ≈
                      2 * A * partials[:, j] + 3 * partials[:, j]
            end

            # Scaling by a `Dual` must not be vectorized component-wise
            α = make_dual(2.0, fill(0.5, num_partials))
            mul!(du, D, u, α)
            @test ForwardDiff.value.(du) ≈ 2 * A * values
            for j in 1:num_partials
                @test ForwardDiff.partials.(du, j) ≈
                      2 * A * partials[:, j] + 0.5 * A * values
            end

            # Views must work as well
            u_long = [make_dual(0.0, zeros(num_partials)) for _ in 1:(N + 4)]
            u_view = view(u_long, 3:(N + 2))
            copyto!(u_view, u)
            du_view = view(similar(u_long), 3:(N + 2))
            mul!(du_view, D, u_view)
            @test ForwardDiff.value.(du_view) ≈ A * values
        end
    end

    # Nested `Dual`s
    D = derivative_operator(MattssonNordström2004(); derivative_order = 1,
                            accuracy_order = 4, xmin = 0.0, xmax = 1.0, N = 20)
    A = Matrix(D)
    values = randn(20)
    u = map(x -> ForwardDiff.Dual{:outer}(ForwardDiff.Dual{:inner}(x, 2 * x),
                                          ForwardDiff.Dual{:inner}(3 * x, 4 * x)),
            values)
    du = similar(u)
    mul!(du, D, u)
    @test ForwardDiff.value.(ForwardDiff.value.(du)) ≈ A * values
    @test ForwardDiff.partials.(ForwardDiff.value.(du), 1) ≈ A * (2 * values)
    @test ForwardDiff.value.(ForwardDiff.partials.(du, 1)) ≈ A * (3 * values)
    @test ForwardDiff.partials.(ForwardDiff.partials.(du, 1), 1) ≈ A * (4 * values)

    # Scaling factors that are no native numbers must not end up inside the
    # vectorized loops over the components
    u = map(x -> ForwardDiff.Dual{:ad_test}(x, 2 * x), values)
    du = similar(u)
    mul!(du, D, u, big(2.0))
    @test ForwardDiff.value.(du) ≈ 2 * (A * values)
    @test ForwardDiff.partials.(du, 1) ≈ 2 * (A * (2 * values))
end

@testset "Jacobian-vector product" begin
    function StructDual(x::AbstractVector{T}, w::AbstractVector{T}) where {T}
        @assert length(x) == length(w)
        # This was the original suggestion. However, it is currently not stable
        # under broadcasting. Thus, we use a slightly different version.
        # partials = StructArray{ForwardDiff.Partials{1, T}}(
        #     (StructArray{Tuple{T}}(
        #         (w,)
        #     ),)
        # )
        partials = reinterpret(reshape, ForwardDiff.Partials{1, T}, w)
        duals = StructArray{ForwardDiff.Dual{Nothing, T, 1}}((x, partials))
        return duals
    end

    function ForwardDiff.value(dx::StructArray{D}) where {D <: ForwardDiff.Dual}
        return dx.value
    end

    function ForwardDiff.partials(dx::StructArray{<:ForwardDiff.Dual{Tag, T, 1}},
                                  i) where {Tag, T}
        # This was the original suggestion. We need to update it (see above).
        # return getproperty(dx.partials.values, i)
        @assert i == 1
        return reinterpret(reshape, T, dx.partials)
    end

    @testset "fourier_derivative_operator" begin
        D = fourier_derivative_operator(xmin = 0.0, xmax = 1.0, N = 8)

        u = randn(size(D, 2))
        v = randn(size(D, 2))
        u_v = StructDual(u, v)
        f_df = @inferred(D*u_v)
        @test ForwardDiff.value(f_df) ≈ @inferred(D*u)
        @test ForwardDiff.partials(f_df, 1) ≈ @inferred(D*v)

        f = let D = D
            f(u) = u .* (D * (u .^ 2))
        end
        f_df = f(u_v)
        J = Diagonal(D * u .^ 2) + 2 .* u .* Matrix(D) * Diagonal(u)
        @test ForwardDiff.value(f_df) ≈ f(u)
        @test ForwardDiff.partials(f_df, 1) ≈ J * v
    end

    @testset "FourierPolynomialDerivativeOperator" begin
        D = fourier_derivative_operator(xmin = 0.0, xmax = 1.0, N = 8)
        D = I - D^2

        u = randn(size(D, 2))
        v = randn(size(D, 2))
        u_v = StructDual(u, v)
        f_df = @inferred(D*u_v)
        @test ForwardDiff.value(f_df) ≈ @inferred(D*u)
        @test ForwardDiff.partials(f_df, 1) ≈ @inferred(D*v)

        f = let D = D
            f(u) = u .* (D * (u .^ 2))
        end
        f_df = f(u_v)
        J = Diagonal(D * u .^ 2) + 2 .* u .* Matrix(D) * Diagonal(u)
        @test ForwardDiff.value(f_df) ≈ f(u)
        @test ForwardDiff.partials(f_df, 1) ≈ J * v
    end

    @testset "FourierRationalDerivativeOperator" begin
        D = fourier_derivative_operator(xmin = 0.0, xmax = 1.0, N = 8)
        D = inv(I - D^2)

        u = randn(size(D, 2))
        v = randn(size(D, 2))
        u_v = StructDual(u, v)
        f_df = @inferred(D*u_v)
        @test ForwardDiff.value(f_df) ≈ @inferred(D*u)
        @test ForwardDiff.partials(f_df, 1) ≈ @inferred(D*v)

        f = let D = D
            f(u) = u .* (D * (u .^ 2))
        end
        f_df = f(u_v)
        J = Diagonal(D * u .^ 2) + 2 .* u .* Matrix(D) * Diagonal(u)
        @test ForwardDiff.value(f_df) ≈ f(u)
        @test ForwardDiff.partials(f_df, 1) ≈ J * v
    end

    @testset "PeriodicRationalDerivativeOperator, 1" begin
        D = periodic_derivative_operator(derivative_order = 1, accuracy_order = 4,
                                         xmin = 0.0, xmax = 2.0, N = 20)
        D = I - D^2

        u = randn(size(D, 2))
        v = randn(size(D, 2))
        u_v = StructDual(u, v)
        f_df = @inferred(D*u_v)
        @test ForwardDiff.value(f_df) ≈ @inferred(D*u)
        @test ForwardDiff.partials(f_df, 1) ≈ @inferred(D*v)

        f = let D = D
            f(u) = u .* (D * (u .^ 2))
        end
        f_df = f(u_v)
        J = Diagonal(D * u .^ 2) + 2 .* u .* Matrix(D) * Diagonal(u)
        @test ForwardDiff.value(f_df) ≈ f(u)
        @test ForwardDiff.partials(f_df, 1) ≈ J * v
    end

    @testset "PeriodicRationalDerivativeOperator, 2" begin
        D = periodic_derivative_operator(derivative_order = 1, accuracy_order = 4,
                                         xmin = 0.0, xmax = 2.0, N = 20)
        D = inv(I - D^2)

        u = randn(size(D, 2))
        v = randn(size(D, 2))
        u_v = StructDual(u, v)
        f_df = @inferred(D*u_v)
        @test ForwardDiff.value(f_df) ≈ @inferred(D*u)
        @test ForwardDiff.partials(f_df, 1) ≈ @inferred(D*v)

        f = let D = D
            f(u) = u .* (D * (u .^ 2))
        end
        f_df = f(u_v)
        J = Diagonal(D * u .^ 2) + 2 .* u .* Matrix(D) * Diagonal(u)
        @test ForwardDiff.value(f_df) ≈ f(u)
        @test ForwardDiff.partials(f_df, 1) ≈ J * v
    end
end

end # module
