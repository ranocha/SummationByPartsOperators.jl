using Test
using LinearAlgebra
using SummationByPartsOperators

# check construction of interior part of upwind operators
@testset "Check against some upwind operators" begin
  N = 14
  xmin_construction = 0.5
  xmax_construction = 1.0
  xmin_application = -0.25
  xmax_application = 0.75

  for acc_order in 2:6
    Dm_bounded = derivative_operator(Mattsson2017(:minus  ), 1, acc_order,
                                     xmin_construction, xmax_construction, N)
    Dp_bounded = derivative_operator(Mattsson2017(:plus   ), 1, acc_order,
                                     xmin_construction, xmax_construction, N)
    Dc_bounded = derivative_operator(Mattsson2017(:central), 1, acc_order,
                                     xmin_construction, xmax_construction, N)

    nodes = collect(grid(Dc_bounded))
    weights = diag(mass_matrix(Dc_bounded))
    Dm = MatrixDerivativeOperator(xmin_application, xmax_application,
                                  nodes, weights, Matrix(Dm_bounded), acc_order)
    Dp = MatrixDerivativeOperator(xmin_application, xmax_application,
                                  nodes, weights, Matrix(Dp_bounded), acc_order)
    Dc = MatrixDerivativeOperator(xmin_application, xmax_application,
                                  nodes, weights, Matrix(Dc_bounded), acc_order)
    D = UpwindOperators(Dm, Dc, Dp)

    for compact in (true, false)
      show(IOContext(devnull, :compact=>compact), D)
      show(IOContext(devnull, :compact=>compact), Dm)
      show(IOContext(devnull, :compact=>compact), Dp)
      show(IOContext(devnull, :compact=>compact), Dc)
      summary(IOContext(devnull, :compact=>compact), D)
      summary(IOContext(devnull, :compact=>compact), Dm)
      summary(IOContext(devnull, :compact=>compact), Dp)
      summary(IOContext(devnull, :compact=>compact), Dc)
    end
    M = mass_matrix(D)
    @test M == mass_matrix(Dm)
    @test M == mass_matrix(Dp)
    @test M == mass_matrix(Dc)

    x = grid(D)
    @test x == grid(Dm)
    @test x == grid(Dp)
    @test x == grid(Dc)

    D_reference = upwind_operators(Mattsson2017;
                                   derivative_order = 1, accuracy_order = acc_order,
                                   xmin = xmin_application, xmax = xmax_application,
                                   N = N)
    @test x ≈ grid(D_reference)
    @test M ≈ mass_matrix(D_reference)

    u = sinpi.(x)
    @test D.minus * u ≈ D_reference.minus * u
    @test D.plus * u ≈ D_reference.plus * u
    @test D.central * u ≈ D_reference.central * u
  end
end
