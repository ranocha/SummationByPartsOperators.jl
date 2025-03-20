using Test
using SummationByPartsOperators
using LinearAlgebra
using StaticArrays
using SparseArrays

@testset "Check against MattssonNordström2004() in 1D" begin
    N = 14
    xmin_construction = 0.5
    xmax_construction = 1.0

    for acc_order in [2, 4, 6]
        D = derivative_operator(MattssonNordström2004(), 1, acc_order, xmin_construction,
                                xmax_construction, N)

        nodes = SVector.(grid(D))
        boundary_indices = [1, N]
        normals = [SVector(-1.0), SVector(1.0)]
        weights = diag(mass_matrix(D))
        weights_boundary = [1.0, 1.0]
        Ds_dense = (Matrix(D),)
        D_multi_dense = MultidimensionalMatrixDerivativeOperator(nodes, boundary_indices,
                                                                 normals, weights,
                                                                 weights_boundary, Ds_dense,
                                                                 acc_order,
                                                                 source_of_coefficients(D))
        Ds_sparse = (sparse(D),)
        D_multi_sparse = MultidimensionalMatrixDerivativeOperator(nodes, boundary_indices,
                                                                  normals, weights,
                                                                  weights_boundary,
                                                                  Ds_sparse, acc_order,
                                                                  source_of_coefficients(D))
        @test D_multi_sparse[1] isa SparseMatrixCSC

        for D_multi in (D_multi_dense, D_multi_sparse)
            for compact in (true, false)
                show(IOContext(devnull, :compact => compact), D_multi)
                summary(IOContext(devnull, :compact => compact), D_multi)
            end

            @test ndims(D_multi) == 1
            @test derivative_order(D_multi) == 1
            @test accuracy_order(D_multi) == acc_order
            @test source_of_coefficients(D_multi) == source_of_coefficients(D)
            @test eltype(D_multi) == eltype(D)

            @test grid(D_multi) == SVector.(grid(D))
            M = mass_matrix(D_multi)
            @test M == mass_matrix(D)
            D_x = D_multi[1]
            @test Matrix(D_x) == Matrix(D)
            @test sparse(D_x) == sparse(D)
            B_x = mass_matrix_boundary(D_multi, 1)
            @test B_x == mass_matrix_boundary(D)
            @test M * D_x + D_x' * M ≈ B_x

            x = grid(D)
            u = sinpi.(x)
            u_reference = copy(u)

            @test integrate(abs2, u, D_multi) ≈ integrate(abs2, u, D)
            @test integrate_boundary(abs2, u, D_multi, 1) ≈ integrate_boundary(abs2, u, D)
            @test integrate_boundary(u, D_multi, 1) ≈ integrate_boundary(u, D)

            SummationByPartsOperators.scale_by_mass_matrix!(u, D_multi)
            @test_throws DimensionMismatch scale_by_mass_matrix!(@view(u[(begin + 1):(end - 1)]),
                                                                 D_multi)
            SummationByPartsOperators.scale_by_mass_matrix!(u_reference, D)
            @test_throws DimensionMismatch scale_by_mass_matrix!(@view(u_reference[(begin + 1):(end - 1)]),
                                                                 D)
            @test u ≈ u_reference
            SummationByPartsOperators.scale_by_inverse_mass_matrix!(u, D_multi)
            @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(@view(u[(begin + 1):(end - 1)]),
                                                                         D_multi)
            SummationByPartsOperators.scale_by_inverse_mass_matrix!(u_reference, D)
            @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(@view(u_reference[(begin + 1):(end - 1)]),
                                                                         D)
            @test u ≈ u_reference

            @test SummationByPartsOperators.weights_boundary(D_multi) == [1.0, 1.0]
            @test SummationByPartsOperators.get_weight(D_multi, 1) ==
                  left_boundary_weight(D_multi) == left_boundary_weight(D)
            @test SummationByPartsOperators.get_weight(D_multi, N) ==
                  right_boundary_weight(D_multi) == right_boundary_weight(D)

            @test SummationByPartsOperators.lower_bandwidth(D_multi) ==
                  size(D_multi[1], 1) - 1
            @test SummationByPartsOperators.upper_bandwidth(D_multi) ==
                  size(D_multi[1], 1) - 1
        end
    end
end

@testset "2D tensor product operators" begin
    N_x = 14
    N_y = 18
    xmin_construction = 0.5
    xmax_construction = 1.0
    ymin_construction = -1.0
    ymax_construction = -0.5

    for acc_order in [2, 4, 6], T in (Float32, Float64)
        D_1 = derivative_operator(MattssonNordström2004(), 1, acc_order,
                                  T(xmin_construction), T(xmax_construction), N_x)
        D_2 = derivative_operator(MattssonAlmquistCarpenter2014Extended(), 1, acc_order,
                                  T(ymin_construction), T(ymax_construction), N_y)
        D_t = tensor_product_operator_2D(D_1, D_2)

        for compact in (true, false)
            show(IOContext(devnull, :compact => compact), D_t)
            summary(IOContext(devnull, :compact => compact), D_t)
        end

        @test ndims(D_t) == 2
        @test derivative_order(D_t) == 1
        @test accuracy_order(D_t) == acc_order
        @test source_of_coefficients(D_t) isa SourceOfCoefficientsCombination{2}
        for compact in (true, false)
            show(IOContext(devnull, :compact => compact), source_of_coefficients(D_t))
            summary(IOContext(devnull, :compact => compact), source_of_coefficients(D_t))
        end
        @test eltype(D_t) == eltype(D_1) == eltype(D_2) == T
        @test real(D_t) == real(D_1) == real(D_2) == T
        @test normals(D_t) == D_t.normals
        @test boundary_indices(D_t) == D_t.boundary_indices

        nodes = grid(D_t)
        @test length(nodes) == N_x * N_y
        u = sin.(first.(nodes) .^ 2) .* cospi.(last.(nodes))
        u_reference = copy(u)

        SummationByPartsOperators.scale_by_mass_matrix!(u, D_t)
        SummationByPartsOperators.scale_by_inverse_mass_matrix!(u, D_t)
        @test u ≈ u_reference
        @test length(restrict_interior(u, D_t)) == (N_x - 2) * (N_y - 2)
        @test length(restrict_boundary(u, D_t)) == 2 * (N_x + N_y)

        # derivative operators
        M = mass_matrix(D_t)
        D_x = D_t[1]
        @test D_x isa SparseMatrixCSC

        D_y = D_t[2]
        @test D_y isa SparseMatrixCSC

        M_1D_1 = mass_matrix(D_1)
        M_1D_2 = mass_matrix(D_2)
        Q_1D_1 = M_1D_1 * Matrix(D_1)
        Q_x = kron(Q_1D_1, M_1D_2)
        Q_1D_2 = M_1D_2 * Matrix(D_2)
        Q_y = kron(M_1D_1, Q_1D_2)
        @test Q_x ≈ M * D_x
        @test Q_y ≈ M * D_y

        B_x = mass_matrix_boundary(D_t, 1)
        B_y = mass_matrix_boundary(D_t, 2)
        @test M * D_x + D_x' * M ≈ B_x
        @test M * D_y + D_y' * M ≈ B_y
        B_1D_1 = mass_matrix_boundary(D_1)
        B_1D_2 = mass_matrix_boundary(D_2)
        @test B_x ≈ Diagonal(kron(B_1D_1, M_1D_2))
        @test B_y ≈ Diagonal(kron(M_1D_1, B_1D_2))

        # integrals
        @test integrate(abs2, u, D_t) ≈ sum(M * abs2.(u))
        @test integrate_boundary(abs2, u, D_t, 1) ≈ sum(B_x * abs2.(u))
        @test integrate_boundary(abs2, u, D_t, 2) ≈ sum(B_y * abs2.(u))

        # accuracy test
        x = grid(D_t)
        x0 = ones(N_x * N_y)
        x1 = first.(x)
        y1 = last.(x)
        res = D_t[1] * x0
        @test all(i -> res[i] < 1000 * eps(T), eachindex(res))
        res = D_t[2] * x0
        @test all(i -> res[i] < 1000 * eps(T), eachindex(res))
        res = D_t[1] * x1
        @test all(i -> isapprox(res[i], x0[i], atol = 1000 * eps(T)), eachindex(res))
        res = D_t[1] * y1
        @test all(i -> res[i] < 1000 * eps(T), eachindex(res))
        res = D_t[2] * x1
        @test all(i -> res[i] < 1000 * eps(T), eachindex(res))
        res = D_t[2] * y1
        @test all(i -> isapprox(res[i], x0[i], atol = 1000 * eps(T)), eachindex(res))
        if acc_order >= 4
            res = D_t[1] * x1 .^ 2
            @test all(i -> isapprox(res[i], 2 * x1[i], atol = 1000 * eps(T)),
                      eachindex(res))
            res = D_t[2] * y1 .^ 2
            @test all(i -> isapprox(res[i], 2 * y1[i], atol = 1000 * eps(T)),
                      eachindex(res))
        end
        if acc_order >= 6
            res = D_t[1] * x1 .^ 3
            @test all(i -> isapprox(res[i], 3 * x1[i]^2, atol = 1000 * eps(T)),
                      eachindex(res))
            res = D_t[2] * y1 .^ 3
            @test all(i -> isapprox(res[i], 3 * y1[i]^2, atol = 1000 * eps(T)),
                      eachindex(res))
        end
    end
end
