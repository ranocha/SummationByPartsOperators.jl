module TestAD

using ForwardDiff
using StructArrays
using SummationByPartsOperators

using LinearAlgebra: Diagonal, I
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

  function ForwardDiff.partials(dx::StructArray{<: ForwardDiff.Dual{Tag, T, 1}}, i) where {Tag, T}
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
    f_df = @inferred(D * u_v)
    @test ForwardDiff.value(f_df) ≈ @inferred(D * u)
    @test ForwardDiff.partials(f_df, 1) ≈ @inferred(D * v)

    f = let D = D
      f(u) = u .* (D * (u.^2))
    end
    f_df = f(u_v)
    J = Diagonal(D * u.^2) + 2 .* u .* Matrix(D) * Diagonal(u)
    @test ForwardDiff.value(f_df) ≈ f(u)
    @test ForwardDiff.partials(f_df, 1) ≈ J * v
  end

  @testset "FourierPolynomialDerivativeOperator" begin
    D = fourier_derivative_operator(xmin = 0.0, xmax = 1.0, N = 8)
    D = I - D^2

    u = randn(size(D, 2))
    v = randn(size(D, 2))
    u_v = StructDual(u, v)
    f_df = @inferred(D * u_v)
    @test ForwardDiff.value(f_df) ≈ @inferred(D * u)
    @test ForwardDiff.partials(f_df, 1) ≈ @inferred(D * v)

    f = let D = D
      f(u) = u .* (D * (u.^2))
    end
    f_df = f(u_v)
    J = Diagonal(D * u.^2) + 2 .* u .* Matrix(D) * Diagonal(u)
    @test ForwardDiff.value(f_df) ≈ f(u)
    @test ForwardDiff.partials(f_df, 1) ≈ J * v
  end

  @testset "FourierRationalDerivativeOperator" begin
    D = fourier_derivative_operator(xmin = 0.0, xmax = 1.0, N = 8)
    D = inv(I - D^2)

    u = randn(size(D, 2))
    v = randn(size(D, 2))
    u_v = StructDual(u, v)
    f_df = @inferred(D * u_v)
    @test ForwardDiff.value(f_df) ≈ @inferred(D * u)
    @test ForwardDiff.partials(f_df, 1) ≈ @inferred(D * v)

    f = let D = D
      f(u) = u .* (D * (u.^2))
    end
    f_df = f(u_v)
    J = Diagonal(D * u.^2) + 2 .* u .* Matrix(D) * Diagonal(u)
    @test ForwardDiff.value(f_df) ≈ f(u)
    @test ForwardDiff.partials(f_df, 1) ≈ J * v
  end
end

end # module
