module SpecialMatrixTypesTest

using Test
using LinearAlgebra
using StaticArrays
using StructArrays
using SummationByPartsOperators

struct Vec4{T} <: FieldVector{4,T}
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
  D = derivative_operator(MattssonNordström2004(), derivative_order=1, accuracy_order=4,
                          xmin=zero(T), xmax=one(T), N=51)
  D_sparse = sparse(D)

  # 3-arg mul!
  u_aos_plain = randn(T, 4, size(D, 1)); du_aos_plain = similar(u_aos_plain)
  mul_aos!(du_aos_plain, D, u_aos_plain)

  u_aos_r = reinterpret(reshape, Vec4{T}, u_aos_plain); du_aos_r = similar(u_aos_r)
  mul!(du_aos_r, D, u_aos_r)
  @test reinterpret(reshape, T, du_aos_r) ≈ du_aos_plain
  @test D * u_aos_r ≈ D_sparse * u_aos_r

  u_aos = Array(u_aos_r); du_aos = similar(u_aos)
  mul!(du_aos, D, u_aos)
  @test du_aos ≈ du_aos_r
  @test D * u_aos ≈ D_sparse * u_aos

  u_soa = StructArray(u_aos); du_soa = similar(u_soa)
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

  fill!(du_soa, zero(eltype(du_soa)))
  mul!(du_soa, D, u_soa, α)
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

  fill!(du_soa, zero(eltype(du_soa)))
  mul!(du_soa, D, u_soa, α, β)
  @test du_soa ≈ du_aos
end

end # module
