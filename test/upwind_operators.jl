using Test
using LinearAlgebra
using SummationByPartsOperators

# check construction of interior part of upwind operators
@testset "Check interior parts" begin
  N = 21
  xmin = 0//1
  xmax = (N + 1) // 1
  interior = 10:N-10

  acc_order = 2
  Dp_bounded = derivative_operator(Mattsson2017(:plus   ), 1, acc_order, xmin, xmax, N)
  Dm_bounded = derivative_operator(Mattsson2017(:minus  ), 1, acc_order, xmin, xmax, N)
  Dc_bounded = derivative_operator(Mattsson2017(:central), 1, acc_order, xmin, xmax, N)
  println(devnull, Dp_bounded)
  M = mass_matrix(Dp_bounded)
  @test M == mass_matrix(Dm_bounded)
  @test M == mass_matrix(Dc_bounded)
  Dp_periodic = periodic_derivative_operator(1, acc_order, xmin, xmax, N, -(acc_order - 1) ÷ 2)
  Dm_periodic = periodic_derivative_operator(1, acc_order, xmin, xmax, N, -acc_order + (acc_order - 1) ÷ 2)
  Dp = Matrix(Dp_bounded)
  Dm = Matrix(Dm_bounded)
  @test Dp[interior,interior] ≈ Matrix(Dp_periodic)[interior,interior]
  @test Dm[interior,interior] ≈ Matrix(Dm_periodic)[interior,interior]
  res = M * Dp + Dm' * M
  res[1,1] += 1
  res[end,end] -= 1
  @test norm(res) < 10N * eps()
end
