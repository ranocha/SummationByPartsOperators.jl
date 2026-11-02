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

@testset "Boundary mass matrix with nodes appearing several times" begin
    # The corners of a tensor product operator appear twice in `boundary_indices`, once for
    # each adjacent face. The boundary mass matrix must not depend on the order of these
    # entries, so we build a `MultidimensionalMatrixDerivativeOperator` from the same data
    # with the boundary entries in reverse order, which also swaps the two entries of each
    # corner.
    N_x = 10
    N_y = 12
    for acc_order in [2, 4], T in (Float32, Float64)
        D_1 = derivative_operator(MattssonNordström2004(), 1, acc_order,
                                  T(0.0), T(1.0), N_x)
        D_2 = derivative_operator(MattssonNordström2004(), 1, acc_order,
                                  T(-1.0), T(2.0), N_y)
        D_t = tensor_product_operator_2D(D_1, D_2)
        @test length(unique(boundary_indices(D_t))) < length(boundary_indices(D_t))
        # The wrapped operator is concretely typed, so forwarding to it is type stable
        @test isconcretetype(fieldtype(typeof(D_t), :D_multi))
        @inferred grid(D_t)
        @inferred mass_matrix_boundary(D_t, 1)

        perm = reverse(eachindex(D_t.boundary_indices))
        D_multi = MultidimensionalMatrixDerivativeOperator(grid(D_t),
                                                           D_t.boundary_indices[perm],
                                                           D_t.normals[perm],
                                                           D_t.weights,
                                                           D_t.weights_boundary[perm],
                                                           D_t.Ds,
                                                           accuracy_order(D_t),
                                                           source_of_coefficients(D_t))
        M = mass_matrix(D_multi)
        # boundary integrals clearly nonzero in both directions, so that `≈` is meaningful
        u = exp.(first.(grid(D_t)) .+ 2 .* last.(grid(D_t)))
        for dim in 1:2
            B = mass_matrix_boundary(D_multi, dim)
            @test B == mass_matrix_boundary(D_t, dim)
            @test M * D_multi[dim] + D_multi[dim]' * M ≈ B
            @test integrate_boundary(abs2, u, D_multi, dim) ≈ sum(B * abs2.(u))
        end
    end
end

@testset "Boundary mass matrix with several nonzero contributions at corners" begin
    # We map a tensor product operator on the unit square affinely to a parallelogram,
    # x = J ξ. Then, the normals of both faces adjacent to a corner have nonzero components
    # in both directions, i.e., the boundary mass matrix receives two nonzero contributions
    # at each corner in each direction. By the chain rule ∇_x = J^{-T} ∇_ξ and Nanson's
    # formula n ds = det(J) J^{-T} n̂ dŝ, the mapped operator satisfies the SBP property.
    N_x = 9
    N_y = 10
    for acc_order in [2, 4], T in (Float32, Float64)
        D_1 = derivative_operator(MattssonNordström2004(), 1, acc_order,
                                  zero(T), one(T), N_x)
        D_2 = derivative_operator(MattssonNordström2004(), 1, acc_order,
                                  zero(T), one(T), N_y)
        D_t = tensor_product_operator_2D(D_1, D_2)

        J = SMatrix{2, 2, T}(2, 1, 1, 3)
        det_J = det(J)
        J_inv_T = inv(J)'
        nodes = Ref(J) .* grid(D_t)
        weights = det_J .* D_t.weights
        # boundary weight times normal
        scaled_normals = Ref(det_J * J_inv_T) .* D_t.normals
        normals = normalize.(scaled_normals)
        weights_boundary = D_t.weights_boundary .* norm.(scaled_normals)
        Ds = ntuple(i -> J_inv_T[i, 1] * D_t[1] + J_inv_T[i, 2] * D_t[2], 2)
        D = MultidimensionalMatrixDerivativeOperator(nodes, D_t.boundary_indices, normals,
                                                     weights, weights_boundary, Ds,
                                                     accuracy_order(D_t),
                                                     source_of_coefficients(D_t))

        # positions of the entries of the corners in `boundary_indices`
        corner_entries = findall(k -> count(==(k), D.boundary_indices) > 1,
                                 D.boundary_indices)
        @test length(corner_entries) == 8
        @test length(unique(D.boundary_indices[corner_entries])) == 4
        for j in corner_entries, dim in 1:2
            @test abs(weights_boundary[j] * normals[j][dim]) > sqrt(eps(T))
        end

        M = mass_matrix(D)
        x = first.(grid(D))
        y = last.(grid(D))
        @test D[1] * x ≈ ones(T, length(x))
        @test D[2] * y ≈ ones(T, length(y))
        @test maximum(abs, D[1] * y) < 1000 * eps(T)
        @test maximum(abs, D[2] * x) < 1000 * eps(T)
        u = exp.(x .- y)
        for dim in 1:2
            B = mass_matrix_boundary(D, dim)
            @test M * D[dim] + D[dim]' * M ≈ B
            @test B ≈
                  det_J * (J_inv_T[dim, 1] * mass_matrix_boundary(D_t, 1) +
                   J_inv_T[dim, 2] * mass_matrix_boundary(D_t, 2))
            @test integrate_boundary(abs2, u, D, dim) ≈ sum(B * abs2.(u))
        end
    end
end

@testset "2D tensor product operators with different types of 1D operators" begin
    D_fd = derivative_operator(MattssonNordström2004(), 1, 4, -1.0, 1.0, 12)
    D_lgl = legendre_derivative_operator(0.0, 2.0, 6)
    for (D_1, D_2) in ((D_fd, D_lgl), (D_lgl, D_fd))
        D_t = tensor_product_operator_2D(D_1, D_2)
        source = source_of_coefficients(D_t)
        @test source isa SourceOfCoefficientsCombination{2}
        @test source.sources ==
              (source_of_coefficients(D_1), source_of_coefficients(D_2))
        for compact in (true, false)
            show(IOContext(devnull, :compact => compact), source)
            summary(IOContext(devnull, :compact => compact), source)
        end
        @test accuracy_order(D_t) == min(accuracy_order(D_1), accuracy_order(D_2))

        M = mass_matrix(D_t)
        M_1D_1 = mass_matrix(D_1)
        M_1D_2 = mass_matrix(D_2)
        @test M ≈ kron(M_1D_1, M_1D_2)
        @test M * D_t[1] ≈ kron(M_1D_1 * Matrix(D_1), M_1D_2)
        @test M * D_t[2] ≈ kron(M_1D_1, M_1D_2 * Matrix(D_2))
        B_x = mass_matrix_boundary(D_t, 1)
        B_y = mass_matrix_boundary(D_t, 2)
        @test M * D_t[1] + D_t[1]' * M ≈ B_x
        @test M * D_t[2] + D_t[2]' * M ≈ B_y
        @test B_x ≈ Diagonal(kron(mass_matrix_boundary(D_1), M_1D_2))
        @test B_y ≈ Diagonal(kron(M_1D_1, mass_matrix_boundary(D_2)))
    end
end
