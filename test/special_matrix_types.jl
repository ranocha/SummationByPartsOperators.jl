module SpecialMatrixTypesTest

using Test
using LinearAlgebra
using SparseArrays
using StaticArrays
using StructArrays
using SummationByPartsOperators

struct Vec4{T} <: FieldVector{4, T}
    x1::T
    x2::T
    x3::T
    x4::T
end

function mul_aos!(du, D, u, args...)
    for i in 1:size(du, 1)
        mul!(view(du, i, :), D, view(u, i, :), args...)
    end
end

for T in (Float32, Float64)
    D = derivative_operator(MattssonNordström2004(), derivative_order = 1,
                            accuracy_order = 4,
                            xmin = zero(T), xmax = one(T), N = 51)
    D_threaded = derivative_operator(MattssonNordström2004(), derivative_order = 1,
                                     accuracy_order = 4,
                                     xmin = zero(T), xmax = one(T), N = 51,
                                     mode = ThreadedMode())
    D_safe = derivative_operator(MattssonNordström2004(), derivative_order = 1,
                                 accuracy_order = 4,
                                 xmin = zero(T), xmax = one(T), N = 51, mode = SafeMode())
    D_sparse = sparse(D)

    # 3-arg mul!
    u_aos_plain = randn(T, 4, size(D, 1))
    du_aos_plain = similar(u_aos_plain)
    mul_aos!(du_aos_plain, D, u_aos_plain)

    u_aos_r = reinterpret(reshape, Vec4{T}, u_aos_plain)
    du_aos_r = similar(u_aos_r)
    mul!(du_aos_r, D, u_aos_r)
    @test reinterpret(reshape, T, du_aos_r) ≈ du_aos_plain
    u_aos_r_sparse = D_sparse * u_aos_r
    @test D * u_aos_r ≈ u_aos_r_sparse
    @test u_aos_r_sparse ≈ D_threaded * u_aos_r
    @test u_aos_r_sparse ≈ D_safe * u_aos_r

    u_aos = Array(u_aos_r)
    du_aos = similar(u_aos)
    mul!(du_aos, D, u_aos)
    @test du_aos ≈ du_aos_r
    @test D * u_aos ≈ D_sparse * u_aos

    u_soa = StructArray(u_aos)
    du_soa = similar(u_soa)
    mul!(du_soa, D, u_soa)
    @test du_soa ≈ du_aos
    @test D * u_soa ≈ D_sparse * u_soa

    # 4-arg mul!
    α = 2 * one(T)
    fill!(du_aos_plain, zero(eltype(du_aos_plain)))
    mul_aos!(du_aos_plain, D, u_aos_plain, α)

    fill!(du_aos_r, zero(eltype(du_aos_r)))
    mul!(du_aos_r, D, u_aos_r, α)
    @test reinterpret(reshape, T, du_aos_r) ≈ du_aos_plain

    fill!(du_aos, zero(eltype(du_aos)))
    mul!(du_aos, D, u_aos, α)
    @test du_aos ≈ du_aos_r

    fill!(du_aos, zero(eltype(du_aos)))
    mul!(du_aos, D_threaded, u_aos, α)
    @test du_aos ≈ du_aos_r

    fill!(du_aos, zero(eltype(du_aos)))
    mul!(du_aos, D_safe, u_aos, α)
    @test du_aos ≈ du_aos_r

    fill!(du_soa, zero(eltype(du_soa)))
    mul!(du_soa, D, u_soa, α)
    @test du_soa ≈ du_aos

    fill!(du_soa, zero(eltype(du_soa)))
    mul!(du_soa, D_threaded, u_soa, α)
    @test du_soa ≈ du_aos

    fill!(du_soa, zero(eltype(du_soa)))
    mul!(du_soa, D_safe, u_soa, α)
    @test du_soa ≈ du_aos

    # 5-arg mul!
    α = 2 * one(T)
    β = 3 * one(T)
    fill!(du_aos_plain, zero(eltype(du_aos_plain)))
    mul_aos!(du_aos_plain, D, u_aos_plain, α, β)

    fill!(du_aos_r, zero(eltype(du_aos_r)))
    mul!(du_aos_r, D, u_aos_r, α, β)
    @test reinterpret(reshape, T, du_aos_r) ≈ du_aos_plain

    fill!(du_aos, zero(eltype(du_aos)))
    mul!(du_aos, D, u_aos, α, β)
    @test du_aos ≈ du_aos_r

    fill!(du_aos, zero(eltype(du_aos)))
    mul!(du_aos, D_threaded, u_aos, α, β)
    @test du_aos ≈ du_aos_r

    fill!(du_aos, zero(eltype(du_aos)))
    mul!(du_aos, D_safe, u_aos, α, β)
    @test du_aos ≈ du_aos_r

    fill!(du_soa, zero(eltype(du_soa)))
    mul!(du_soa, D, u_soa, α, β)
    @test du_soa ≈ du_aos

    fill!(du_soa, zero(eltype(du_soa)))
    mul!(du_soa, D_threaded, u_soa, α, β)
    @test du_soa ≈ du_aos

    fill!(du_soa, zero(eltype(du_soa)))
    mul!(du_soa, D_safe, u_soa, α, β)
    @test du_soa ≈ du_aos

    # `SVector` with only one element
    u_scalar = randn(T, size(D, 1))
    du_scalar = similar(u_scalar)
    u_vector = reinterpret(SVector{1, T}, u_scalar)
    du_vector = similar(u_vector)

    mul!(du_scalar, D, u_scalar)
    mul!(du_vector, D, u_vector)
    @test du_scalar ≈ reinterpret(T, du_vector)

    α = 2 * one(T)
    mul!(du_scalar, D, u_scalar, α)
    mul!(du_vector, D, u_vector, α)
    @test du_scalar ≈ reinterpret(T, du_vector)

    β = 3 * one(T)
    mul!(du_scalar, D, u_scalar, α, β)
    mul!(du_vector, D, u_vector, α, β)
    @test du_scalar ≈ reinterpret(T, du_vector)
end

# Vectors of `Complex` numbers are reinterpreted as matrices of their real and
# imaginary parts, cf.
# https://github.com/ranocha/SummationByPartsOperators.jl/issues/421
@testset "Complex element types" begin
    for T in (ComplexF32, ComplexF64), accuracy_order in (2, 4), N in (17, 30),
        mode in (FastMode(), SafeMode(), ThreadedMode())

        for D in (derivative_operator(MattssonNordström2004();
                                      derivative_order = 1,
                                      accuracy_order = accuracy_order,
                                      xmin = 0.0, xmax = 1.0, N = N, mode = mode),
                  periodic_derivative_operator(derivative_order = 1,
                                               accuracy_order = accuracy_order,
                                               xmin = 0.0, xmax = 1.0, N = N,
                                               mode = mode))
            A = Matrix(D)
            u = rand(T, N)
            du = similar(u)

            mul!(du, D, u)
            @test du ≈ A * u

            mul!(du, D, u, 2)
            @test du ≈ 2 * A * u

            dest = copy(u)
            mul!(dest, D, u, 2, 3)
            @test dest ≈ 2 * A * u + 3 * u

            # Complex scaling factors mix the real and imaginary parts and must
            # therefore not be applied component-wise
            α = 2 + 3im
            mul!(du, D, u, α)
            @test du ≈ α * (A * u)

            β = 4 - 5im
            dest = copy(u)
            mul!(dest, D, u, α, β)
            @test dest ≈ α * (A * u) + β * u
        end
    end
end

# Element types that cannot be reinterpreted as arrays of native numbers must be
# handled by plain loops over the values, cf.
# https://github.com/ranocha/SummationByPartsOperators.jl/issues/421
@testset "Element types without a native memory layout" begin
    N = 20
    # the coefficients of the operators are `Float64`s
    rtol = sqrt(eps(Float64))
    for accuracy_order in (2, 4), mode in (FastMode(), SafeMode(), ThreadedMode())
        for D in (derivative_operator(MattssonNordström2004();
                                      derivative_order = 1,
                                      accuracy_order = accuracy_order,
                                      xmin = 0.0, xmax = 1.0, N = N, mode = mode),
                  periodic_derivative_operator(derivative_order = 1,
                                               accuracy_order = accuracy_order,
                                               xmin = 0.0, xmax = 1.0, N = N,
                                               mode = mode))
            A = Matrix(D)
            # `BigFloat`s are no `isbits` types, so neither vectors of `BigFloat`s
            # nor vectors of `Complex{BigFloat}`s or `SVector`s of `BigFloat`s
            # can be reinterpreted. The same holds for vectors of the mutable
            # `MVector`s.
            for u in (big.(randn(N)),
                      big.(randn(N)) .+ im .* big.(randn(N)),
                      [SVector(big(randn()), big(randn())) for _ in 1:N],
                      [MVector(randn(), randn()) for _ in 1:N])
                du = similar(u)
                mul!(du, D, u)
                @test isapprox(du, A * u; rtol = rtol)

                mul!(du, D, u, 2)
                @test isapprox(du, 2 * (A * u); rtol = rtol)

                dest = copy(u)
                mul!(dest, D, u, 2, 3)
                @test isapprox(dest, 2 * (A * u) + 3 * u; rtol = rtol)
            end
        end
    end
end

# Scaling factors that are no native numbers must not end up inside the
# vectorized loops, cf.
# https://github.com/ranocha/SummationByPartsOperators.jl/issues/421
@testset "Scaling factors without a native number type" begin
    N = 20
    rtol = sqrt(eps(Float64))
    α = big(2.0)
    β = big(3.0)
    for T in (Float64, ComplexF64, SVector{2, Float64}),
        mode in (FastMode(), SafeMode(), ThreadedMode())

        for D in (derivative_operator(MattssonNordström2004();
                                      derivative_order = 1, accuracy_order = 4,
                                      xmin = 0.0, xmax = 1.0, N = N, mode = mode),
                  periodic_derivative_operator(derivative_order = 1, accuracy_order = 4,
                                               xmin = 0.0, xmax = 1.0, N = N,
                                               mode = mode))
            A = Matrix(D)
            u = [rand(T) for _ in 1:N]

            du = similar(u)
            mul!(du, D, u, α)
            @test isapprox(du, α * (A * u); rtol = rtol)

            dest = copy(u)
            mul!(dest, D, u, α, β)
            @test isapprox(dest, α * (A * u) + β * u; rtol = rtol)
        end
    end
end

# User-defined scalar types that are neither `Real` nor `Complex` need not
# support `real`, cf.
# https://github.com/ranocha/SummationByPartsOperators.jl/issues/421
struct Quaternion{T} <: Number
    s::T
    i::T
    j::T
    k::T
end
function Base.zero(::Type{Quaternion{T}}) where {T}
    Quaternion(zero(T), zero(T), zero(T), zero(T))
end
Base.one(::Type{Quaternion{T}}) where {T} = Quaternion(one(T), zero(T), zero(T), zero(T))
function Base.:+(a::Quaternion, b::Quaternion)
    Quaternion(a.s + b.s, a.i + b.i, a.j + b.j, a.k + b.k)
end
Base.:*(x::Real, a::Quaternion) = Quaternion(x * a.s, x * a.i, x * a.j, x * a.k)
Base.:*(a::Quaternion, x::Real) = x * a
# `mul!(dest, D, u)` scales by `one(Quaternion)`, which requires the full product
function Base.:*(a::Quaternion, b::Quaternion)
    Quaternion(a.s * b.s - a.i * b.i - a.j * b.j - a.k * b.k,
               a.s * b.i + a.i * b.s + a.j * b.k - a.k * b.j,
               a.s * b.j - a.i * b.k + a.j * b.s + a.k * b.i,
               a.s * b.k + a.i * b.j - a.j * b.i + a.k * b.s)
end
Base.muladd(x::Real, a::Quaternion, b::Quaternion) = x * a + b

@testset "User-defined scalar element types" begin
    N = 20
    for accuracy_order in (2, 4), mode in (FastMode(), SafeMode(), ThreadedMode())
        for D in (derivative_operator(MattssonNordström2004();
                                      derivative_order = 1,
                                      accuracy_order = accuracy_order,
                                      xmin = 0.0, xmax = 1.0, N = N, mode = mode),
                  periodic_derivative_operator(derivative_order = 1,
                                               accuracy_order = accuracy_order,
                                               xmin = 0.0, xmax = 1.0, N = N,
                                               mode = mode))
            # the components are transformed independently of each other
            components = ntuple(_ -> randn(N), 4)
            u = map(Quaternion, components...)
            names = (:s, :i, :j, :k)

            du = similar(u)
            mul!(du, D, u)
            for (component, name) in zip(components, names)
                @test getproperty.(du, name) ≈ D * component
            end

            for (component, name) in zip(components, names)
                @test getproperty.(D * u, name) ≈ D * component
            end

            mul!(du, D, u, 2)
            for (component, name) in zip(components, names)
                @test getproperty.(du, name) ≈ 2 * (D * component)
            end

            dest = copy(u)
            mul!(dest, D, u, 2, 3)
            for (component, name) in zip(components, names)
                @test getproperty.(dest, name) ≈ 2 * (D * component) + 3 * component
            end
        end
    end
end

end # module
