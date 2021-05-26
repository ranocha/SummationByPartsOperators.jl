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

function mul_aos!(du, D, u)
  for i in 1:size(du, 1)
    mul!(view(du, i, :), D, view(u, i, :))
  end
end

for T in (Float32, Float64)
  D = derivative_operator(MattssonNordström2004(), derivative_order=1, accuracy_order=4,
                          xmin=zero(T), xmax=one(T), N=51)
  x = grid(D)

  u_aos_plain = randn(T, 4, length(x)); du_aos_plain = similar(u_aos_plain)
  mul_aos!(du_aos_plain, D, u_aos_plain)

  u_aos_r = reinterpret(reshape, Vec4{T}, u_aos_plain); du_aos_r = similar(u_aos_r)
  mul!(du_aos_r, D, u_aos_r)
  @test reinterpret(reshape, T, du_aos_r) ≈ du_aos_plain

  u_aos = Array(u_aos_r); du_aos = similar(u_aos)
  mul!(du_aos, D, u_aos)
  @test du_aos ≈ du_aos_r

  u_soa = StructArray(u_aos); du_soa = similar(u_soa)
  mul!(du_soa, D, u_soa)
  @test du_soa ≈ du_aos
end

end # module
