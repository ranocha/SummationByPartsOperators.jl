using Test
using LinearAlgebra
using SummationByPartsOperators

D_test_list = vcat(
  [(Mattsson2017, acc_order) for acc_order in 2:9],
  [(WilliamsDuru2024, acc_order) for acc_order in 4:7]
)

# check construction of interior part of upwind operators
@testset "Check interior parts" begin
  N = 21
  xmin = 0.0
  xmax = Float64(N + 1)

  # bounded derivative operators
  for (source, acc_order) in D_test_list
    Dp_bounded = derivative_operator(source(:plus   ), 1, acc_order, xmin, xmax, N)
    Dm_bounded = derivative_operator(source(:minus  ), 1, acc_order, xmin, xmax, N)
    Dc_bounded = derivative_operator(source(:central), 1, acc_order, xmin, xmax, N)
    for compact in (true, false)
      show(IOContext(devnull, :compact=>compact), Dp_bounded)
      show(IOContext(devnull, :compact=>compact), Dm_bounded)
      show(IOContext(devnull, :compact=>compact), Dc_bounded)
      summary(IOContext(devnull, :compact=>compact), Dp_bounded)
      summary(IOContext(devnull, :compact=>compact), Dm_bounded)
      summary(IOContext(devnull, :compact=>compact), Dc_bounded)
    end
    M = mass_matrix(Dp_bounded)
    @test M == mass_matrix(Dm_bounded)
    @test M == mass_matrix(Dc_bounded)
    Dp = Matrix(Dp_bounded)
    Dm = Matrix(Dm_bounded)
    Dc = Matrix(Dc_bounded)
    @test Dc ≈ (Dm + Dp) / 2

    # upwind summation by parts
    res = M * Dp + Dm' * M
    res[1,1] += 1
    res[end,end] -= 1
    @test norm(res) < N * eps()

    # derivative accuracy order
    for D in (Dp_bounded, Dm_bounded, Dc_bounded)
      x = grid(D)
      interior = length(D.coefficients.left_boundary)+1:N-length(D.coefficients.right_boundary)-1
      @test norm(D * x.^0) < N * eps()
      for k in 1:acc_order÷2
        @test D * x.^k ≈ k .* x.^(k-1)
      end
      for k in acc_order÷2+1:acc_order
        @test (D * x.^k)[interior] ≈ (k .* x.^(k-1))[interior]
      end
    end

    # dissipation matrix
    diss = M * (Dp - Dm)
    @test diss ≈ diss'
    @test maximum(eigvals(Symmetric(diss))) < N * eps()
  end

  # periodic derivative operators
  for acc_order in 2:9
    Dp_bounded = derivative_operator(Mattsson2017(:plus   ), 1, acc_order, xmin, xmax, N)
    Dm_bounded = derivative_operator(Mattsson2017(:minus  ), 1, acc_order, xmin, xmax, N)
    Dc_bounded = derivative_operator(Mattsson2017(:central), 1, acc_order, xmin, xmax, N)

    Dp_periodic = periodic_derivative_operator(1, acc_order, xmin, xmax, N-1, -(acc_order - 1) ÷ 2)
    Dm_periodic = periodic_derivative_operator(1, acc_order, xmin, xmax, N-1, -acc_order + (acc_order - 1) ÷ 2)

    interior = length(Dc_bounded.coefficients.left_boundary)+1:N-length(Dc_bounded.coefficients.right_boundary)-1
    @test Matrix(Dp_bounded)[interior,interior] ≈ Matrix(Dp_periodic)[interior,interior]
    @test Matrix(Dm_bounded)[interior,interior] ≈ Matrix(Dm_periodic)[interior,interior]

    D_periodic = upwind_operators(periodic_derivative_operator;
                                  accuracy_order = acc_order,
                                  xmin, xmax, N = N - 1)
    @test D_periodic.minus == Dm_periodic
    @test D_periodic.plus == Dp_periodic
    @test Matrix(D_periodic.central) ≈ (Matrix(Dm_periodic) + Matrix(Dp_periodic)) / 2

    # mass matrix scaling
    x1 = grid(D_periodic)
    M = @inferred mass_matrix(D_periodic)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D_periodic)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin + 1):(end - 1)]), D_periodic)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D_periodic)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(@view(v[(begin + 1):(end - 1)]), D_periodic)
    @test v ≈ u
  end
end

@testset "Rational Tests" begin
  N = 21
  xmin = big(0 // 1)
  xmax = big((N + 1) // 1)

  # bounded derivative operators
  for (source, acc_order) in D_test_list
    Dp_bounded = derivative_operator(source(:plus   ), 1, acc_order, xmin, xmax, N)
    Dm_bounded = derivative_operator(source(:minus  ), 1, acc_order, xmin, xmax, N)
    Dc_bounded = derivative_operator(source(:central), 1, acc_order, xmin, xmax, N)

    M = mass_matrix(Dp_bounded)
    Dp = Matrix(Dp_bounded)
    Dm = Matrix(Dm_bounded)
    Dc = Matrix(Dc_bounded)
    @test Dc == (Dm + Dp) / 2

    # upwind summation by parts
    res = M * Dp + Dm' * M
    res[1,1] += 1
    res[end,end] -= 1
    @test all(isequal(0), res)

    # derivative accuracy order
    for D in (Dp_bounded, Dm_bounded, Dc_bounded)
      x = grid(D)
      interior = length(D.coefficients.left_boundary)+1:N-length(D.coefficients.right_boundary)-1
      @test all(isequal(0), D * x.^0)
      for k in 1:acc_order÷2
        @test D * x.^k == k .* x.^(k-1)
      end
      for k in acc_order÷2+1:acc_order
        @test (D * x.^k)[interior] == (k .* x.^(k-1))[interior]
      end
    end
  end
end

@testset "UpwindOperators" begin
  N = 21
  xmin = 0.
  xmax = Float64(N + 1)
  acc_order = 2

  D = upwind_operators(Mattsson2017, derivative_order=1, accuracy_order=acc_order,
                       xmin=xmin, xmax=xmax, N=N)
  Dp = derivative_operator(Mattsson2017(:plus   ), 1, acc_order, xmin, xmax, N)
  Dm = derivative_operator(Mattsson2017(:minus  ), 1, acc_order, xmin, xmax, N)
  Dc = derivative_operator(Mattsson2017(:central), 1, acc_order, xmin, xmax, N)
  for compact in (true, false)
    show(IOContext(devnull, :compact=>compact), D)
    show(IOContext(devnull, :compact=>compact), Dp)
    show(IOContext(devnull, :compact=>compact), Dm)
    show(IOContext(devnull, :compact=>compact), Dc)
    summary(IOContext(devnull, :compact=>compact), D)
    summary(IOContext(devnull, :compact=>compact), Dp)
    summary(IOContext(devnull, :compact=>compact), Dm)
    summary(IOContext(devnull, :compact=>compact), Dc)
  end
  @test D.minus == Dm
  @test D.plus == Dp
  @test D.central == Dc
  @test derivative_order(D) == 1

  @test grid(D) == grid(Dp)
  @test SummationByPartsOperators.xmin(D) == SummationByPartsOperators.xmin(Dp)
  @test SummationByPartsOperators.xmax(D) == SummationByPartsOperators.xmax(Dp)

  @test mass_matrix(D) == mass_matrix(Dp)
  @test left_boundary_weight(D) == left_boundary_weight(Dp)
  @test right_boundary_weight(D) == right_boundary_weight(Dp)

  x = grid(D)
  @test integrate(x, D) == integrate(x, Dp)

  M = mass_matrix(D)
  @test M * Matrix(Dp) + Matrix(Dm)' * M ≈ mass_matrix_boundary(D)

  @test_throws ArgumentError UpwindOperators(
    derivative_operator(Mattsson2017(:minus  ), 1, acc_order, xmin, xmax, N),
    derivative_operator(Mattsson2017(:central), 1, acc_order, xmin, xmax, N+1),
    derivative_operator(Mattsson2017(:plus   ), 1, acc_order, xmin, xmax, N)
  )

  # mass matrix scaling
  x1 = grid(D)
  M = @inferred mass_matrix(D)
  u = sinpi.(x1)
  v = copy(u)
  scale_by_mass_matrix!(v, D)
  @test v ≈ M * u
  scale_by_inverse_mass_matrix!(v, D)
  @test v ≈ u
end

@testset "Empty lower/upper coefficients" begin
  D = upwind_operators(periodic_derivative_operator, accuracy_order = 2,
                       xmin = 0.0, xmax = 1.0, N = 10)
  x = grid(D)
  u = @. sinpi(2 * x)

  du = zero(u)
  @test_nowarn mul!(du, D.minus, u)
  @test du ≈ D.minus * u
  @test_nowarn mul!(du, D.minus, u, 2.0)
  @test du ≈ 2 * D.minus * u
  @test_nowarn mul!(du, D.minus, u, 2.0, 3.0)
  @test du ≈ 8 * D.minus * u # 5 = 2 + 3 * 2

  du = zero(u)
  @test_nowarn mul!(du, D.plus, u)
  @test du ≈ D.plus * u
  @test_nowarn mul!(du, D.plus, u, 2.0)
  @test du ≈ 2 * D.plus * u
  @test_nowarn mul!(du, D.plus, u, 2.0, 3.0)
  @test du ≈ 8 * D.plus * u # 5 = 2 + 3 * 2
end
