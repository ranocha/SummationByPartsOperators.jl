using Test
using SparseArrays
using LinearAlgebra
using SummationByPartsOperators

# check construction of interior part of upwind operators
@testset "Check against some upwind operators (dense)" begin
    N = 14
    xmin_construction = 0.5
    xmax_construction = 1.0
    xmin_application = -0.25
    xmax_application = 0.75

    for acc_order in 2:6
        Dm_bounded = derivative_operator(Mattsson2017(:minus),
                                         1,
                                         acc_order,
                                         xmin_construction,
                                         xmax_construction,
                                         N)
        Dp_bounded = derivative_operator(Mattsson2017(:plus),
                                         1,
                                         acc_order,
                                         xmin_construction,
                                         xmax_construction,
                                         N)
        Dc_bounded = derivative_operator(Mattsson2017(:central),
                                         1,
                                         acc_order,
                                         xmin_construction,
                                         xmax_construction,
                                         N)

        nodes = collect(grid(Dc_bounded))
        weights = diag(mass_matrix(Dc_bounded))
        Dm = MatrixDerivativeOperator(xmin_application,
                                      xmax_application,
                                      nodes,
                                      weights,
                                      Matrix(Dm_bounded),
                                      acc_order,
                                      source_of_coefficients(Dm_bounded))
        Dp = MatrixDerivativeOperator(xmin_application,
                                      xmax_application,
                                      nodes,
                                      weights,
                                      Matrix(Dp_bounded),
                                      acc_order,
                                      source_of_coefficients(Dp_bounded))
        Dc = MatrixDerivativeOperator(xmin_application,
                                      xmax_application,
                                      nodes,
                                      weights,
                                      Matrix(Dc_bounded),
                                      acc_order,
                                      source_of_coefficients(Dc_bounded))
        D = UpwindOperators(Dm, Dc, Dp)

        for compact in (true, false)
            show(IOContext(devnull, :compact => compact), D)
            show(IOContext(devnull, :compact => compact), Dm)
            show(IOContext(devnull, :compact => compact), Dp)
            show(IOContext(devnull, :compact => compact), Dc)
            summary(IOContext(devnull, :compact => compact), D)
            summary(IOContext(devnull, :compact => compact), Dm)
            summary(IOContext(devnull, :compact => compact), Dp)
            summary(IOContext(devnull, :compact => compact), Dc)
        end
        M = mass_matrix(D)
        @test M == mass_matrix(Dm)
        @test M == mass_matrix(Dp)
        @test M == mass_matrix(Dc)

        x = grid(D)
        @test x == grid(Dm)
        @test x == grid(Dp)
        @test x == grid(Dc)

        @test issymmetric(Dm) == false
        @test issymmetric(Dp) == false
        @test issymmetric(Dc) == false

        D_reference = upwind_operators(Mattsson2017;
                                       derivative_order = 1,
                                       accuracy_order = acc_order,
                                       xmin = xmin_application,
                                       xmax = xmax_application,
                                       N = N,)
        @test x ≈ grid(D_reference)
        @test M ≈ mass_matrix(D_reference)

        u = sinpi.(x)
        @test D.minus * u ≈ D_reference.minus * u
        @test D.plus * u ≈ D_reference.plus * u
        @test D.central * u ≈ D_reference.central * u

        du = copy(u)
        du_reference = copy(u)
        α = 1.23
        mul!(du, D.minus, u, α)
        mul!(du_reference, D_reference.minus, u, α)
        @test du ≈ du_reference
        mul!(du, D.plus, u, α)
        mul!(du_reference, D_reference.plus, u, α)
        @test du ≈ du_reference
        mul!(du, D.central, u, α)
        mul!(du_reference, D_reference.central, u, α)
        @test du ≈ du_reference

        du = copy(u)
        du_reference = copy(u)
        β = 4.56
        mul!(du, D.minus, u, α, β)
        mul!(du_reference, D_reference.minus, u, α, β)
        @test du ≈ du_reference
        mul!(du, D.plus, u, α, β)
        mul!(du_reference, D_reference.plus, u, α, β)
        @test du ≈ du_reference
        mul!(du, D.central, u, α, β)
        mul!(du_reference, D_reference.central, u, α, β)
        @test du ≈ du_reference

        @test integrate(abs2, u, Dm) ≈ integrate(abs2, u, D_reference.minus)
        @test integrate(abs2, u, Dp) ≈ integrate(abs2, u, D_reference.plus)
        @test integrate(abs2, u, Dc) ≈ integrate(abs2, u, D_reference.central)

        u = sinpi.(x)
        u_reference = copy(u)
        SummationByPartsOperators.scale_by_mass_matrix!(u, Dm)
        @test_throws DimensionMismatch scale_by_mass_matrix!(@view(u[(begin + 1):(end - 1)]),
                                                             Dm)
        SummationByPartsOperators.scale_by_mass_matrix!(u_reference, D_reference.minus)
        @test_throws DimensionMismatch scale_by_mass_matrix!(@view(u_reference[(begin + 1):(end - 1)]),
                                                             D_reference.minus)
        @test u ≈ u_reference
        SummationByPartsOperators.scale_by_inverse_mass_matrix!(u, Dm)
        @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(@view(u[(begin + 1):(end - 1)]),
                                                                     Dm)
        SummationByPartsOperators.scale_by_inverse_mass_matrix!(u_reference,
                                                                D_reference.minus)
        @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(@view(u_reference[(begin + 1):(end - 1)]),
                                                                     D_reference.minus)
        @test u ≈ u_reference

        @test SummationByPartsOperators.get_weight(Dm, 1) == left_boundary_weight(D)
        @test SummationByPartsOperators.get_weight(Dm, N) == right_boundary_weight(D)

        @test SummationByPartsOperators.lower_bandwidth(Dm) == size(Dm, 1) - 1
        @test SummationByPartsOperators.upper_bandwidth(Dm) == size(Dm, 1) - 1
    end
end

@testset "Check against some upwind operators (sparse)" begin
    N = 14
    xmin_construction = 0.5
    xmax_construction = 1.0
    xmin_application = -0.25
    xmax_application = 0.75

    for acc_order in 2:6
        Dm_bounded = derivative_operator(Mattsson2017(:minus),
                                         1,
                                         acc_order,
                                         xmin_construction,
                                         xmax_construction,
                                         N)
        Dp_bounded = derivative_operator(Mattsson2017(:plus),
                                         1,
                                         acc_order,
                                         xmin_construction,
                                         xmax_construction,
                                         N)
        Dc_bounded = derivative_operator(Mattsson2017(:central),
                                         1,
                                         acc_order,
                                         xmin_construction,
                                         xmax_construction,
                                         N)

        nodes = collect(grid(Dc_bounded))
        weights = diag(mass_matrix(Dc_bounded))
        Dm = MatrixDerivativeOperator(xmin_application,
                                      xmax_application,
                                      nodes,
                                      weights,
                                      sparse(Dm_bounded),
                                      acc_order,
                                      source_of_coefficients(Dm_bounded))
        Dp = MatrixDerivativeOperator(xmin_application,
                                      xmax_application,
                                      nodes,
                                      weights,
                                      sparse(Dp_bounded),
                                      acc_order,
                                      source_of_coefficients(Dp_bounded))
        Dc = MatrixDerivativeOperator(xmin_application,
                                      xmax_application,
                                      nodes,
                                      weights,
                                      sparse(Dc_bounded),
                                      acc_order,
                                      source_of_coefficients(Dc_bounded))
        D = UpwindOperators(Dm, Dc, Dp)

        for compact in (true, false)
            show(IOContext(devnull, :compact => compact), D)
            show(IOContext(devnull, :compact => compact), Dm)
            show(IOContext(devnull, :compact => compact), Dp)
            show(IOContext(devnull, :compact => compact), Dc)
            summary(IOContext(devnull, :compact => compact), D)
            summary(IOContext(devnull, :compact => compact), Dm)
            summary(IOContext(devnull, :compact => compact), Dp)
            summary(IOContext(devnull, :compact => compact), Dc)
        end
        M = mass_matrix(D)
        @test M == mass_matrix(Dm)
        @test M == mass_matrix(Dp)
        @test M == mass_matrix(Dc)

        x = grid(D)
        @test x == grid(Dm)
        @test x == grid(Dp)
        @test x == grid(Dc)

        @test issymmetric(Dm) == false
        @test issymmetric(Dp) == false
        @test issymmetric(Dc) == false

        D_reference = upwind_operators(Mattsson2017;
                                       derivative_order = 1,
                                       accuracy_order = acc_order,
                                       xmin = xmin_application,
                                       xmax = xmax_application,
                                       N = N,)
        @test x ≈ grid(D_reference)
        @test M ≈ mass_matrix(D_reference)

        u = sinpi.(x)
        @test D.minus * u ≈ D_reference.minus * u
        @test D.plus * u ≈ D_reference.plus * u
        @test D.central * u ≈ D_reference.central * u

        du = copy(u)
        du_reference = copy(u)
        α = 1.23
        mul!(du, D.minus, u, α)
        mul!(du_reference, D_reference.minus, u, α)
        @test du ≈ du_reference
        mul!(du, D.plus, u, α)
        mul!(du_reference, D_reference.plus, u, α)
        @test du ≈ du_reference
        mul!(du, D.central, u, α)
        mul!(du_reference, D_reference.central, u, α)
        @test du ≈ du_reference

        du = copy(u)
        du_reference = copy(u)
        β = 4.56
        mul!(du, D.minus, u, α, β)
        mul!(du_reference, D_reference.minus, u, α, β)
        @test du ≈ du_reference
        mul!(du, D.plus, u, α, β)
        mul!(du_reference, D_reference.plus, u, α, β)
        @test du ≈ du_reference
        mul!(du, D.central, u, α, β)
        mul!(du_reference, D_reference.central, u, α, β)
        @test du ≈ du_reference

        @test integrate(abs2, u, Dm) ≈ integrate(abs2, u, D_reference.minus)
        @test integrate(abs2, u, Dp) ≈ integrate(abs2, u, D_reference.plus)
        @test integrate(abs2, u, Dc) ≈ integrate(abs2, u, D_reference.central)

        u = sinpi.(x)
        u_reference = copy(u)
        SummationByPartsOperators.scale_by_mass_matrix!(u, Dm)
        @test_throws DimensionMismatch scale_by_mass_matrix!(@view(u[(begin + 1):(end - 1)]),
                                                             Dm)
        SummationByPartsOperators.scale_by_mass_matrix!(u_reference, D_reference.minus)
        @test_throws DimensionMismatch scale_by_mass_matrix!(@view(u_reference[(begin + 1):(end - 1)]),
                                                             D_reference.minus)
        @test u ≈ u_reference
        SummationByPartsOperators.scale_by_inverse_mass_matrix!(u, Dm)
        @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(@view(u[(begin + 1):(end - 1)]),
                                                                     Dm)
        SummationByPartsOperators.scale_by_inverse_mass_matrix!(u_reference,
                                                                D_reference.minus)
        @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(@view(u_reference[(begin + 1):(end - 1)]),
                                                                     D_reference.minus)
        @test u ≈ u_reference

        @test SummationByPartsOperators.get_weight(Dm, 1) == left_boundary_weight(D)
        @test SummationByPartsOperators.get_weight(Dm, N) == right_boundary_weight(D)

        @test SummationByPartsOperators.lower_bandwidth(Dm) == size(Dm, 1) - 1
        @test SummationByPartsOperators.upper_bandwidth(Dm) == size(Dm, 1) - 1
    end
end
