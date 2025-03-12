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
        D = derivative_operator(MattssonNordström2004(), 1, acc_order, xmin_construction, xmax_construction, N)

        nodes = SVector.(grid(D))
        on_boundary = fill(false, N)
        on_boundary[1] = true
        on_boundary[end] = true
        normals = [SVector(-1.0), SVector(1.0)]
        weights = diag(mass_matrix(D))
        weights_boundary = [1.0, 1.0]
        Ds_dense = (Matrix(D),)
        D_multi_dense = MultidimensionalMatrixDerivativeOperator(nodes, on_boundary, normals, weights, weights_boundary, Ds_dense, acc_order, source_of_coefficients(D))
        Ds_sparse = (sparse(D),)
        D_multi_sparse = MultidimensionalMatrixDerivativeOperator(nodes, on_boundary, normals, weights, weights_boundary, Ds_sparse, acc_order, source_of_coefficients(D))
        @test D_multi_sparse[1] isa SparseMatrixCSC

        for D_multi in (D_multi_dense, D_multi_sparse)
            for compact in (true, false)
                show(IOContext(devnull, :compact=>compact), D_multi)
                summary(IOContext(devnull, :compact=>compact), D_multi)
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

            SummationByPartsOperators.scale_by_mass_matrix!(u, D_multi)
            @test_throws DimensionMismatch scale_by_mass_matrix!(@view(u[(begin + 1):(end - 1)]), D_multi)
            SummationByPartsOperators.scale_by_mass_matrix!(u_reference, D)
            @test_throws DimensionMismatch scale_by_mass_matrix!(@view(u_reference[(begin + 1):(end - 1)]), D)
            @test u ≈ u_reference
            SummationByPartsOperators.scale_by_inverse_mass_matrix!(u, D_multi)
            @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(@view(u[(begin + 1):(end - 1)]), D_multi)
            SummationByPartsOperators.scale_by_inverse_mass_matrix!(u_reference, D)
            @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(@view(u_reference[(begin + 1):(end - 1)]), D)
            @test u ≈ u_reference

            @test SummationByPartsOperators.weights_boundary(D_multi) == [1.0, zeros(eltype(D_multi), N - 2)..., 1.0]
            @test SummationByPartsOperators.get_weight(D_multi, 1) == left_boundary_weight(D_multi) == left_boundary_weight(D)
            @test SummationByPartsOperators.get_weight(D_multi, N) == right_boundary_weight(D_multi) == right_boundary_weight(D)

            @test SummationByPartsOperators.lower_bandwidth(D_multi) == size(D_multi[1], 1) - 1
            @test SummationByPartsOperators.upper_bandwidth(D_multi) == size(D_multi[1], 1) - 1
        end
    end
end

@testset "2D tensor product operators" begin
    N = 14
    xmin_construction = 0.5
    xmax_construction = 1.0

    for acc_order in [2, 4, 6]
        D = derivative_operator(MattssonNordström2004(), 1, acc_order, xmin_construction, xmax_construction, N)
        D_t = tensor_product_operator_2D(D)

        for compact in (true, false)
            show(IOContext(devnull, :compact=>compact), D_t)
            summary(IOContext(devnull, :compact=>compact), D_t)
        end

        @test ndims(D_t) == 2
        @test derivative_order(D_t) == 1
        @test accuracy_order(D_t) == acc_order
        @test source_of_coefficients(D_t) == source_of_coefficients(D)
        @test eltype(D_t) == eltype(D)

        @test length(grid(D_t)) == N^2
        M = mass_matrix(D_t)
        D_x = D_t[1]
        @test D_x isa SparseMatrixCSC
        B_x = mass_matrix_boundary(D_t, 1)
        @test M * D_x + D_x' * M ≈ B_x
        D_y = D_t[2]
        @test D_y isa SparseMatrixCSC
        B_y = mass_matrix_boundary(D_t, 2)
        @test M * D_y + D_y' * M ≈ B_y

        M_1D = mass_matrix(D)
        Q_1D = M_1D * Matrix(D)
        Q_x = kron(Q_1D, M_1D)
        Q_y = kron(M_1D, Q_1D)
        @test Q_x ≈ M * D_x
        @test Q_y ≈ M * D_y

        B_1D = mass_matrix_boundary(D)
        @test B_x ≈ Diagonal(kron(B_1D, M_1D))
        @test B_y ≈ Diagonal(kron(M_1D, B_1D))

        # accuracy test
        T = eltype(D)
        x = grid(D_t)
        x0 = ones(N^2)
        x1 = first.(x)
        y1 = last.(x)
        res = D_t[1] * x0; @test all(i->res[i] < 1000 * eps(T), eachindex(res))
        res = D_t[2] * x0; @test all(i->res[i] < 1000 * eps(T), eachindex(res))
        res = D_t[1] * x1; @test all(i->res[i] ≈ x0[i], eachindex(res))
        res = D_t[1] * y1; @test all(i->res[i] < 1000 * eps(T), eachindex(res))
        res = D_t[2] * x1; @test all(i->res[i] < 1000 * eps(T), eachindex(res))
        res = D_t[2] * y1; @test all(i->res[i] ≈ x0[i], eachindex(res))
        if acc_order >= 4
            res = D_t[1] * x1.^2; @test all(i->res[i] ≈ 2 * x1[i], eachindex(res))
            res = D_t[2] * y1.^2; @test all(i->res[i] ≈ 2 * y1[i], eachindex(res))
        end
        if acc_order >= 6
            res = D_t[1] * x1.^3; @test all(i->res[i] ≈ 3 * x1[i]^2, eachindex(res))
            res = D_t[2] * y1.^3; @test all(i->res[i] ≈ 3 * y1[i]^2, eachindex(res))
        end
    end
end
