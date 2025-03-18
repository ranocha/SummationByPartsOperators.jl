using Test
using LinearAlgebra, SparseArrays
using SummationByPartsOperators

# Accuracy tests with Float32.
let T = Float32
    xmin = one(T)
    xmax = 2 * one(T)
    N = 99
    x1 = compute_coefficients(identity, periodic_derivative_operator(1, 2, xmin, xmax, N))

    x0 = fill(one(eltype(x1)), length(x1))
    x2 = x1 .* x1
    x3 = x2 .* x1
    x4 = x2 .* x2
    x5 = x2 .* x3
    x6 = x3 .* x3
    x7 = x4 .* x3

    res = fill(zero(eltype(x0)), length(x0))

    # first derivative operators
    derivative_order = 1
    accuracy_order = 2
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    for compact in (true, false)
        show(IOContext(devnull, :compact => false), D)
    end
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == false
    @test SummationByPartsOperators.xmin(D) ≈ xmin
    @test SummationByPartsOperators.xmax(D) ≈ xmax
    M = mass_matrix(D)
    @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
    mul!(res, D, x0)
    @test all(i -> abs(res[i]) < eps(T), accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x1)
    @test all(i -> res[i] ≈ x0[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x2)
    @test all(i -> res[i] ≈ 2 * x1[i], accuracy_order:length(res)-accuracy_order)
    # mass matrix scaling
    M = @inferred mass_matrix(D)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin+1):(end-1)]), D)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(
        @view(v[(begin+1):(end-1)]),
        D,
    )
    @test v ≈ u

    accuracy_order = 4
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    for compact in (true, false)
        show(IOContext(devnull, :compact => false), D)
    end
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == false
    @test SummationByPartsOperators.xmin(D) ≈ xmin
    @test SummationByPartsOperators.xmax(D) ≈ xmax
    M = mass_matrix(D)
    @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
    mul!(res, D, x0)
    @test all(i -> abs(res[i]) < 50 * eps(T), accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x1)
    @test all(i -> res[i] ≈ x0[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x2)
    @test all(i -> res[i] ≈ 2 * x1[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x3)
    @test all(i -> res[i] ≈ 3 * x2[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x4)
    @test all(i -> res[i] ≈ 4 * x3[i], accuracy_order:length(res)-accuracy_order)
    # mass matrix scaling
    M = @inferred mass_matrix(D)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin+1):(end-1)]), D)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(
        @view(v[(begin+1):(end-1)]),
        D,
    )
    @test v ≈ u

    accuracy_order = 6
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    for compact in (true, false)
        show(IOContext(devnull, :compact => false), D)
    end
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == false
    @test SummationByPartsOperators.xmin(D) ≈ xmin
    @test SummationByPartsOperators.xmax(D) ≈ xmax
    M = mass_matrix(D)
    @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
    mul!(res, D, x0)
    @test all(i -> abs(res[i]) < 50 * eps(T), accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x1)
    @test all(i -> res[i] ≈ x0[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x2)
    @test all(i -> res[i] ≈ 2 * x1[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x3)
    @test all(i -> res[i] ≈ 3 * x2[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x4)
    @test all(i -> res[i] ≈ 4 * x3[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x5)
    @test all(i -> res[i] ≈ 5 * x4[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x6)
    @test all(i -> res[i] ≈ 6 * x5[i], accuracy_order:length(res)-accuracy_order)
    # mass matrix scaling
    M = @inferred mass_matrix(D)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin+1):(end-1)]), D)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(
        @view(v[(begin+1):(end-1)]),
        D,
    )
    @test v ≈ u

    # second derivative operators
    derivative_order = 2
    accuracy_order = 2
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    for compact in (true, false)
        show(IOContext(devnull, :compact => false), D)
    end
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == true
    @test SummationByPartsOperators.xmin(D) ≈ xmin
    @test SummationByPartsOperators.xmax(D) ≈ xmax
    M = mass_matrix(D)
    @test M * Matrix(D) - Matrix(D)' * M ≈ zeros(T, N, N)
    mul!(res, D, x0)
    @test all(i -> abs(res[i]) < eps(T), accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x1)
    @test all(i -> abs(res[i]) < 100N * eps(T), accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x2)
    @test all(
        i -> isapprox(res[i], 2 * x0[i], atol = 620N * eps(T)),
        accuracy_order:length(res)-accuracy_order,
    )
    mul!(res, D, x3)
    @test all(
        i -> isapprox(res[i], 6 * x1[i], atol = 2000N * eps(T)),
        accuracy_order:length(res)-accuracy_order,
    )
    mul!(res, D, x4)
    @test any(i -> !(res[i] ≈ 12 * x2[i]), accuracy_order:length(res)-accuracy_order)
    # mass matrix scaling
    M = @inferred mass_matrix(D)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin+1):(end-1)]), D)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(
        @view(v[(begin+1):(end-1)]),
        D,
    )
    @test v ≈ u

    accuracy_order = 4
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    for compact in (true, false)
        show(IOContext(devnull, :compact => false), D)
    end
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == true
    M = mass_matrix(D)
    @test M * Matrix(D) - Matrix(D)' * M ≈ zeros(T, N, N)
    mul!(res, D, x0)
    @test all(i -> abs(res[i]) < 50N * eps(T), accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x1)
    @test all(i -> abs(res[i]) < 1000N * eps(T), accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x2)
    @test all(
        i -> isapprox(res[i], 2 * x0[i], atol = 5000N * eps(T)),
        accuracy_order:length(res)-accuracy_order,
    )
    mul!(res, D, x3)
    @test all(
        i -> isapprox(res[i], 6 * x1[i], atol = 5000N * eps(T)),
        accuracy_order:length(res)-accuracy_order,
    )
    mul!(res, D, x4)
    @test all(
        i -> isapprox(res[i], 12 * x2[i], atol = 8500N * eps(T)),
        accuracy_order:length(res)-accuracy_order,
    )
    mul!(res, D, x5)
    @test any(i -> !(res[i] ≈ 30 * x3[i]), accuracy_order:length(res)-accuracy_order)
    # mass matrix scaling
    M = @inferred mass_matrix(D)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin+1):(end-1)]), D)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(
        @view(v[(begin+1):(end-1)]),
        D,
    )
    @test v ≈ u

    accuracy_order = 6
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    for compact in (true, false)
        show(IOContext(devnull, :compact => false), D)
    end
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == true
    M = mass_matrix(D)
    @test M * Matrix(D) - Matrix(D)' * M ≈ zeros(T, N, N)
    mul!(res, D, x0)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 50N * eps(T)
    mul!(res, D, x1)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 1000N * eps(T)
    mul!(res, D, x2)
    @test norm((res-2 .* x0)[accuracy_order:end-accuracy_order], Inf) < 5000N * eps(T)
    mul!(res, D, x3)
    @test norm((res-6 .* x1)[accuracy_order:end-accuracy_order], Inf) < 5000N * eps(T)
    mul!(res, D, x4)
    @test norm((res-12 .* x2)[accuracy_order:end-accuracy_order], Inf) < 8000N * eps(T)
    mul!(res, D, x5)
    @test norm((res-20 .* x3)[accuracy_order:end-accuracy_order], Inf) < 50000N * eps(T)
    mul!(res, D, x6)
    @test norm((res-30 .* x4)[accuracy_order:end-accuracy_order], Inf) < 52000N * eps(T)
    mul!(res, D, x7)
    @test norm((res-42 .* x5)[accuracy_order:end-accuracy_order], Inf) < 92000N * eps(T)
    # mass matrix scaling
    M = @inferred mass_matrix(D)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin+1):(end-1)]), D)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(
        @view(v[(begin+1):(end-1)]),
        D,
    )
    @test v ≈ u

    tmp = D * x7
    @test typeof(tmp) == typeof(res)
    @test norm(tmp - res, Inf) < eps(T)

    # third derivative operators
    derivative_order = 3
    accuracy_order = 2
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    for compact in (true, false)
        show(IOContext(devnull, :compact => false), D)
    end
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == true # because this operator is zero!
    M = mass_matrix(D)
    @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
    mul!(res, D, x0)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < eps(T)
    mul!(res, D, x1)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 1400 * eps(T)
    mul!(res, D, x2)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 7000 * eps(T)
    mul!(res, D, x3)
    @test norm((res-6 .* x0)[accuracy_order:end-accuracy_order], Inf) > 7000 * eps(T)
    # mass matrix scaling
    M = @inferred mass_matrix(D)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin+1):(end-1)]), D)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(
        @view(v[(begin+1):(end-1)]),
        D,
    )
    @test v ≈ u

    accuracy_order = 4
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    for compact in (true, false)
        show(IOContext(devnull, :compact => false), D)
    end
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == false
    M = mass_matrix(D)
    @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
    mul!(res, D, x0)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 400 * eps(T)
    mul!(res, D, x1)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 10000N * eps(T)
    mul!(res, D, x2)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 50000N * eps(T)
    mul!(res, D, x3)
    @test norm((res-6 .* x0)[accuracy_order:end-accuracy_order], Inf) < 500000N * eps(T)
    mul!(res, D, x4)
    @test norm((res-24 .* x1)[accuracy_order:end-accuracy_order], Inf) < 600000N * eps(T)
    mul!(res, D, x5)
    @test norm((res-60 .* x2)[accuracy_order:end-accuracy_order], Inf) > 600000 * eps(T)
    # mass matrix scaling
    M = @inferred mass_matrix(D)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin+1):(end-1)]), D)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(
        @view(v[(begin+1):(end-1)]),
        D,
    )
    @test v ≈ u

    accuracy_order = 6
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    for compact in (true, false)
        show(IOContext(devnull, :compact => false), D)
    end
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == false
    M = mass_matrix(D)
    @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
    mul!(res, D, x0)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 300 * eps(T)
    mul!(res, D, x1)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 80000N * eps(T)
    mul!(res, D, x2)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 300000N * eps(T)
    mul!(res, D, x3)
    @test norm((res-6 .* x0)[accuracy_order:end-accuracy_order], Inf) < 700000N * eps(T)
    mul!(res, D, x4)
    @test norm((res-24 .* x1)[accuracy_order:end-accuracy_order], Inf) < 1000000N * eps(T)
    mul!(res, D, x5)
    @test norm((res-60 .* x2)[accuracy_order:end-accuracy_order], Inf) < 4000000N * eps(T)
    mul!(res, D, x6)
    @test norm((res-120 .* x3)[accuracy_order:end-accuracy_order], Inf) < 8000000N * eps(T)
    mul!(res, D, x7)
    @test norm((res-240 .* x4)[accuracy_order:end-accuracy_order], Inf) > 220000N * eps(T)

    tmp = D * x7
    @test typeof(tmp) == typeof(res)
    @test norm(tmp - res, Inf) < eps(T)

    # higher derivative operators
    accuracy_order = 2
    @test_throws ArgumentError periodic_central_derivative_operator(
        4,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    @test_throws ArgumentError periodic_central_derivative_operator(
        5,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    @test_throws ArgumentError periodic_central_derivative_operator(
        6,
        accuracy_order,
        xmin,
        xmax,
        N,
    )

    @test integrate(x7, D) ≈ sum(mass_matrix(D) * x7)
    @test integrate(u -> u^2, x7, D) ≈ dot(x7, mass_matrix(D), x7)
    @test integrate_boundary(x7, D) ≈ 0.0
    @test integrate_boundary(u -> u^2, x7, D) ≈ 0.0
end

# Accuracy tests with Float64.
let T = Float64
    xmin = -one(T)
    xmax = 2 * one(T)
    N = 100
    x1 = compute_coefficients(identity, periodic_derivative_operator(1, 2, xmin, xmax, N))

    x0 = fill(one(eltype(x1)), length(x1))
    x2 = x1 .* x1
    x3 = x2 .* x1
    x4 = x2 .* x2
    x5 = x2 .* x3
    x6 = x3 .* x3
    x7 = x4 .* x3
    x8 = x4 .* x4
    x9 = x4 .* x5

    res = fill(zero(eltype(x0)), length(x0))

    # first derivative operators
    derivative_order = 1
    accuracy_order = 2
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    println(devnull, D)
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == false
    M = mass_matrix(D)
    @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
    mul!(res, D, x0)
    @test all(i -> abs(res[i]) < eps(T), accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x1)
    @test all(i -> res[i] ≈ x0[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x2)
    @test all(i -> res[i] ≈ 2 * x1[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x3)
    @test any(i -> !(res[i] ≈ 3 * x2[i]), accuracy_order:length(res)-accuracy_order)
    # mass matrix scaling
    M = @inferred mass_matrix(D)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin+1):(end-1)]), D)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(
        @view(v[(begin+1):(end-1)]),
        D,
    )
    @test v ≈ u

    accuracy_order = 4
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    println(devnull, D)
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == false
    M = mass_matrix(D)
    @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
    mul!(res, D, x0)
    @test all(i -> abs(res[i]) < 10 * eps(T), accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x1)
    @test all(i -> res[i] ≈ x0[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x2)
    @test all(i -> res[i] ≈ 2 * x1[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x3)
    @test all(i -> res[i] ≈ 3 * x2[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x4)
    @test all(i -> res[i] ≈ 4 * x3[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x5)
    @test any(i -> !(res[i] ≈ 5 * x4[i]), accuracy_order:length(res)-accuracy_order)
    # mass matrix scaling
    M = @inferred mass_matrix(D)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin+1):(end-1)]), D)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(
        @view(v[(begin+1):(end-1)]),
        D,
    )
    @test v ≈ u

    accuracy_order = 6
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    println(devnull, D)
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == false
    M = mass_matrix(D)
    @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
    mul!(res, D, x0)
    @test all(i -> abs(res[i]) < 10 * eps(T), accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x1)
    @test all(i -> res[i] ≈ x0[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x2)
    @test all(i -> res[i] ≈ 2 * x1[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x3)
    @test all(
        i -> isapprox(res[i], 3 * x2[i], atol = 100N * eps(T)),
        accuracy_order:length(res)-accuracy_order,
    )
    mul!(res, D, x4)
    @test all(i -> res[i] ≈ 4 * x3[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x5)
    @test all(
        i -> isapprox(res[i], 5 * x4[i], atol = 100N * eps(T)),
        accuracy_order:length(res)-accuracy_order,
    )
    mul!(res, D, x6)
    @test all(i -> res[i] ≈ 6 * x5[i], accuracy_order:length(res)-accuracy_order)
    mul!(res, D, x7)
    @test any(i -> !(res[i] ≈ 7 * x6[i]), accuracy_order:length(res)-accuracy_order)

    tmp = D * x7
    @test typeof(tmp) == typeof(res)
    @test norm(tmp - res, Inf) < eps(T)

    # second derivative operators
    derivative_order = 2
    accuracy_order = 2
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    println(devnull, D)
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == true
    M = mass_matrix(D)
    @test M * Matrix(D) - Matrix(D)' * M ≈ zeros(T, N, N)
    mul!(res, D, x0)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < eps(T)
    mul!(res, D, x1)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 1400 * eps(T)
    mul!(res, D, x2)
    @test norm((res-2 .* x0)[accuracy_order:end-accuracy_order], Inf) < 7000 * eps(T)
    mul!(res, D, x3)
    @test norm((res-6 .* x1)[accuracy_order:end-accuracy_order], Inf) > 7000 * eps(T)
    # mass matrix scaling
    M = @inferred mass_matrix(D)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin+1):(end-1)]), D)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(
        @view(v[(begin+1):(end-1)]),
        D,
    )
    @test v ≈ u

    accuracy_order = 4
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    println(devnull, D)
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == true
    M = mass_matrix(D)
    @test M * Matrix(D) - Matrix(D)' * M ≈ zeros(T, N, N)
    mul!(res, D, x0)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 400 * eps(T)
    mul!(res, D, x1)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 5000 * eps(T)
    mul!(res, D, x2)
    @test norm((res-2 .* x0)[accuracy_order:end-accuracy_order], Inf) < 11000 * eps(T)
    mul!(res, D, x3)
    @test norm((res-6 .* x1)[accuracy_order:end-accuracy_order], Inf) < 23000 * eps(T)
    mul!(res, D, x4)
    @test norm((res-12 .* x2)[accuracy_order:end-accuracy_order], Inf) < 40000 * eps(T)
    mul!(res, D, x5)
    @test norm((res-20 .* x3)[accuracy_order:end-accuracy_order], Inf) > 50000 * eps(T)
    # mass matrix scaling
    M = @inferred mass_matrix(D)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin+1):(end-1)]), D)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(
        @view(v[(begin+1):(end-1)]),
        D,
    )
    @test v ≈ u

    accuracy_order = 6
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    println(devnull, D)
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == true
    M = mass_matrix(D)
    @test M * Matrix(D) - Matrix(D)' * M ≈ zeros(T, N, N)
    mul!(res, D, x0)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 1000 * eps(T)
    mul!(res, D, x1)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 5000 * eps(T)
    mul!(res, D, x2)
    @test norm((res-2 .* x0)[accuracy_order:end-accuracy_order], Inf) < 10020 * eps(T)
    mul!(res, D, x3)
    @test norm((res-6 .* x1)[accuracy_order:end-accuracy_order], Inf) < 30000 * eps(T)
    mul!(res, D, x4)
    @test norm((res-12 .* x2)[accuracy_order:end-accuracy_order], Inf) < 40000 * eps(T)
    mul!(res, D, x5)
    @test norm((res-20 .* x3)[accuracy_order:end-accuracy_order], Inf) < 100000 * eps(T)
    mul!(res, D, x6)
    @test norm((res-30 .* x4)[accuracy_order:end-accuracy_order], Inf) < 220000 * eps(T)
    mul!(res, D, x7)
    @test norm((res-42 .* x5)[accuracy_order:end-accuracy_order], Inf) > 220000 * eps(T)

    tmp = D * x7
    @test typeof(tmp) == typeof(res)
    @test norm(tmp - res, Inf) < eps(T)

    # third derivative operators
    derivative_order = 3
    accuracy_order = 2
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    println(devnull, D)
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == true # because this operator is zero!
    M = mass_matrix(D)
    @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
    mul!(res, D, x0)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < eps(T)
    mul!(res, D, x1)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 1400 * eps(T)
    mul!(res, D, x2)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 7000 * eps(T)
    mul!(res, D, x3)
    @test norm((res-6 .* x0)[accuracy_order:end-accuracy_order], Inf) > 7000 * eps(T)
    # mass matrix scaling
    M = @inferred mass_matrix(D)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin+1):(end-1)]), D)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(
        @view(v[(begin+1):(end-1)]),
        D,
    )
    @test v ≈ u

    accuracy_order = 4
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    println(devnull, D)
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == false
    M = mass_matrix(D)
    @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
    mul!(res, D, x0)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 400 * eps(T)
    mul!(res, D, x1)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 20000 * eps(T)
    mul!(res, D, x2)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 200000 * eps(T)
    mul!(res, D, x3)
    @test norm((res-6 .* x0)[accuracy_order:end-accuracy_order], Inf) < 300000 * eps(T)
    mul!(res, D, x4)
    @test norm((res-24 .* x1)[accuracy_order:end-accuracy_order], Inf) < 600000 * eps(T)
    mul!(res, D, x5)
    @test norm((res-60 .* x2)[accuracy_order:end-accuracy_order], Inf) > 600000 * eps(T)
    # mass matrix scaling
    M = @inferred mass_matrix(D)
    u = sinpi.(x1)
    v = copy(u)
    scale_by_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin+1):(end-1)]), D)
    @test v ≈ M * u
    scale_by_inverse_mass_matrix!(v, D)
    @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(
        @view(v[(begin+1):(end-1)]),
        D,
    )
    @test v ≈ u

    accuracy_order = 6
    D = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    println(devnull, D)
    @test SummationByPartsOperators.derivative_order(D) == derivative_order
    @test SummationByPartsOperators.accuracy_order(D) == accuracy_order
    @test issymmetric(D) == false
    M = mass_matrix(D)
    @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
    mul!(res, D, x0)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 300 * eps(T)
    mul!(res, D, x1)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 70000 * eps(T)
    mul!(res, D, x2)
    @test norm(res[accuracy_order:end-accuracy_order], Inf) < 200000 * eps(T)
    mul!(res, D, x3)
    @test norm((res-6 .* x0)[accuracy_order:end-accuracy_order], Inf) < 600000 * eps(T)
    mul!(res, D, x4)
    @test norm((res-24 .* x1)[accuracy_order:end-accuracy_order], Inf) < 1000000 * eps(T)
    mul!(res, D, x5)
    @test norm((res-60 .* x2)[accuracy_order:end-accuracy_order], Inf) < 4000000 * eps(T)
    mul!(res, D, x6)
    @test norm((res-120 .* x3)[accuracy_order:end-accuracy_order], Inf) < 8000000 * eps(T)
    mul!(res, D, x7)
    @test norm((res-240 .* x4)[accuracy_order:end-accuracy_order], Inf) > 220000 * eps(T)

    tmp = D * x7
    @test typeof(tmp) == typeof(res)
    @test norm(tmp - res, Inf) < eps(T)

    # higher derivative operators
    accuracy_order = 2
    @test_throws ArgumentError periodic_central_derivative_operator(
        4,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    @test_throws ArgumentError periodic_central_derivative_operator(
        5,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    @test_throws ArgumentError periodic_central_derivative_operator(
        6,
        accuracy_order,
        xmin,
        xmax,
        N,
    )

    @test integrate(x7, D) ≈ sum(mass_matrix(D) * x7)
    @test integrate(u -> u^2, x7, D) ≈ dot(x7, mass_matrix(D), x7)
    @test integrate_boundary(x7, D) ≈ 0.0
    @test integrate_boundary(u -> u^2, x7, D) ≈ 0.0
end

# Compare Fornberg algorithm with exact representation of central derivative coefficients.
for T in (Float32, Float64), accuracy_order = 2:2:12, derivative_order = 1:3
    xmin = zero(T)
    xmax = one(T)
    N = 21
    D_central = periodic_central_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
    )
    D_general =
        periodic_derivative_operator(derivative_order, accuracy_order, xmin, xmax, N)
    @test norm(D_central.coefficients.lower_coef - D_general.coefficients.lower_coef) <
          3 * eps(T)
    @test abs(D_central.coefficients.central_coef - D_general.coefficients.central_coef) <
          3 * eps(T)
    @test norm(D_central.coefficients.upper_coef - D_general.coefficients.upper_coef) <
          3 * eps(T)
end

# Compare serial and threaded with full and sparse matrix vector products implementation.
for T in (Float32, Float64), accuracy_order = 1:10, derivative_order = 1:3
    xmin = zero(T)
    xmax = 5 * one(T)
    N = 50
    D_serial = periodic_derivative_operator(derivative_order, accuracy_order, xmin, xmax, N)
    D_threads = periodic_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
        ThreadedMode(),
    )

    u = compute_coefficients(x -> x^5, D_serial)
    dest_serial = fill(zero(eltype(u)), length(u))
    dest_threads = fill(zero(eltype(u)), length(u))
    dest_full = fill(zero(eltype(u)), length(u))
    dest_sparse = fill(zero(eltype(u)), length(u))

    D_full = Matrix(D_serial)
    D_sparse = sparse(D_serial)
    mul!(dest_serial, D_serial, u)
    mul!(dest_threads, D_threads, u)
    mul!(dest_full, D_full, u)
    mul!(dest_sparse, D_sparse, u)
    @test all(i -> dest_serial[i] ≈ dest_threads[i], eachindex(u))
    @test all(
        i -> isapprox(dest_serial[i], dest_full[i], rtol = 5 * sqrt(eps(T))),
        eachindex(u),
    )
    @test all(
        i -> isapprox(dest_serial[i], dest_sparse[i], rtol = 5 * sqrt(eps(T))),
        eachindex(u),
    )
end

# Compare mul! with β=0 and mul! without β.
for T in (Float32, Float64), accuracy_order = 1:10, derivative_order = 1:3
    xmin = zero(T)
    xmax = 5 * one(T)
    N = 51
    D_serial = periodic_derivative_operator(derivative_order, accuracy_order, xmin, xmax, N)
    D_threads = periodic_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
        ThreadedMode(),
    )
    D_safe = periodic_derivative_operator(
        derivative_order,
        accuracy_order,
        xmin,
        xmax,
        N,
        SafeMode(),
    )

    u = compute_coefficients(x -> x^5, D_serial)
    dest1 = fill(zero(eltype(u)), length(u))
    dest2 = fill(zero(eltype(u)), length(u))

    mul!(dest1, D_serial, u, one(T), zero(T))
    mul!(dest2, D_serial, u, one(T))
    @test all(i -> dest1[i] ≈ dest2[i], eachindex(u))
    mul!(dest1, D_threads, u, one(T), zero(T))
    mul!(dest2, D_threads, u, one(T))
    @test all(i -> dest1[i] ≈ dest2[i], eachindex(u))
    @test dest1 ≈ D_safe * u
end

# Accuracy tests for periodic signals
let T = Float32
    xmin = T(-1)
    xmax = T(1)
    N = 50
    res = Array{T}(undef, N)

    # order 2
    accuracy_order = 2
    D = periodic_derivative_operator(1, accuracy_order, xmin, xmax, N)

    factor = 1
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 1.0e-2
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 2
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 0.8
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 5
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 1.2
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    # order 4
    accuracy_order = 4
    D = periodic_derivative_operator(1, accuracy_order, xmin, xmax, N)

    factor = 1
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 5.0e-5
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 2
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 1.0e-3
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 5
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 0.09
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 8
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 0.9
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    # order 6
    accuracy_order = 6
    D = periodic_derivative_operator(1, accuracy_order, xmin, xmax, N)

    factor = 1
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 3.0e-6
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 2
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 3.0e-5
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 5
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 0.008
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 8
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 0.3
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol
end
let T = Float64
    xmin = T(-1)
    xmax = T(1)
    N = 50
    res = Array{T}(undef, N)

    # order 2
    accuracy_order = 2
    D = periodic_derivative_operator(1, accuracy_order, xmin, xmax, N)

    factor = 1
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 0.01
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 2
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 0.09
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 5
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 1.2
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    # order 4
    accuracy_order = 4
    D = periodic_derivative_operator(1, accuracy_order, xmin, xmax, N)

    factor = 1
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 3.0e-5
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 2
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 1.0e-3
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 5
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 0.1
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 8
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 0.9
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    # order 6
    accuracy_order = 6
    D = periodic_derivative_operator(1, accuracy_order, xmin, xmax, N)

    factor = 1
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 1.0e-7
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 2
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 3.0e-5
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 5
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 0.01
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol

    factor = 8
    ufunc = x -> sinpi(factor * x)
    dufunc = x -> factor * T(π) * cospi(factor * x)
    tol = 0.3
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    res = D * u
    @test maximum(abs, res - du) < tol
    xplot, duplot = evaluate_coefficients(res, D)
    @test maximum(abs, duplot - dufunc.(xplot)) < tol
end


# Robust methods of Holoborodko
let T = Float32
    xmin = one(T)
    xmax = 2 * one(T)
    N = 100
    x1 = compute_coefficients(identity, periodic_derivative_operator(1, 2, xmin, xmax, N))

    x0 = fill(one(eltype(x1)), length(x1))
    x2 = x1 .* x1
    x3 = x2 .* x1
    x4 = x2 .* x2
    x5 = x2 .* x3
    x6 = x3 .* x3
    x7 = x4 .* x3

    res = fill(zero(eltype(x0)), length(x0))

    # first derivative operators
    der_order = 1
    acc_order = 2
    for stencil_width = 5:2:11
        D = periodic_derivative_operator(
            Holoborodko2008(),
            der_order,
            acc_order,
            xmin,
            xmax,
            N,
            stencil_width = stencil_width,
        )
        println(devnull, D)
        @test derivative_order(D) == der_order
        @test accuracy_order(D) == acc_order
        @test issymmetric(D) == false
        M = mass_matrix(D)
        @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
        mul!(res, D, x0)
        @test all(i -> abs(res[i]) < eps(T), stencil_width:length(res)-stencil_width)
        mul!(res, D, x1)
        @test all(i -> res[i] ≈ x0[i], stencil_width:length(res)-stencil_width)
        mul!(res, D, x2)
        @test all(i -> res[i] ≈ 2 * x1[i], stencil_width:length(res)-stencil_width)
    end
    @test_throws ArgumentError D = periodic_derivative_operator(
        Holoborodko2008(),
        der_order,
        acc_order,
        xmin,
        xmax,
        N,
        stencil_width = 3,
    )
    @test_throws ArgumentError D = periodic_derivative_operator(
        Holoborodko2008(),
        der_order,
        acc_order,
        xmin,
        xmax,
        N,
        stencil_width = 13,
    )

    acc_order = 4
    for stencil_width = 7:2:11
        D = periodic_derivative_operator(
            Holoborodko2008(),
            der_order,
            acc_order,
            xmin,
            xmax,
            N,
            stencil_width = stencil_width,
        )
        println(devnull, D)
        @test derivative_order(D) == der_order
        @test accuracy_order(D) == acc_order
        @test issymmetric(D) == false
        M = mass_matrix(D)
        @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
        mul!(res, D, x0)
        @test all(i -> abs(res[i]) < 50 * eps(T), stencil_width:length(res)-stencil_width)
        mul!(res, D, x1)
        @test all(i -> res[i] ≈ x0[i], stencil_width:length(res)-stencil_width)
        mul!(res, D, x2)
        @test all(i -> res[i] ≈ 2 * x1[i], stencil_width:length(res)-stencil_width)
        mul!(res, D, x3)
        @test all(i -> res[i] ≈ 3 * x2[i], stencil_width:length(res)-stencil_width)
        mul!(res, D, x4)
        @test all(i -> res[i] ≈ 4 * x3[i], stencil_width:length(res)-stencil_width)
    end
    @test_throws ArgumentError D = periodic_derivative_operator(
        Holoborodko2008(),
        der_order,
        acc_order,
        xmin,
        xmax,
        N,
        stencil_width = 5,
    )
    @test_throws ArgumentError D = periodic_derivative_operator(
        Holoborodko2008(),
        der_order,
        acc_order,
        xmin,
        xmax,
        N,
        stencil_width = 13,
    )

    acc_order = 6
    @test_throws ArgumentError D = periodic_derivative_operator(
        Holoborodko2008(),
        der_order,
        acc_order,
        xmin,
        xmax,
        N,
        stencil_width = 9,
    )
    @test_throws ArgumentError D = periodic_derivative_operator(
        Holoborodko2008(),
        der_order,
        acc_order,
        xmin,
        xmax,
        N,
        stencil_width = 13,
    )

    # second derivative operators
    der_order = 2
    acc_order = 2
    for stencil_width = 5:2:9
        D = periodic_derivative_operator(
            Holoborodko2008(),
            der_order,
            acc_order,
            xmin,
            xmax,
            N,
            stencil_width = stencil_width,
        )
        println(devnull, D)
        @test derivative_order(D) == der_order
        @test accuracy_order(D) == acc_order
        @test issymmetric(D) == true
        M = mass_matrix(D)
        @test M * Matrix(D) - Matrix(D)' * M ≈ zeros(T, N, N)
        mul!(res, D, x0)
        @test all(i -> abs(res[i]) < eps(T), stencil_width:length(res)-stencil_width)
        mul!(res, D, x1)
        @test all(i -> abs(res[i]) < 100N * eps(T), stencil_width:length(res)-stencil_width)
        mul!(res, D, x2)
        @test all(
            i -> isapprox(res[i], 2 * x0[i], atol = 600N * eps(T)),
            stencil_width:length(res)-stencil_width,
        )
        mul!(res, D, x3)
        @test all(
            i -> isapprox(res[i], 6 * x1[i], atol = 2000N * eps(T)),
            stencil_width:length(res)-stencil_width,
        )
    end
    @test_throws ArgumentError D = periodic_derivative_operator(
        Holoborodko2008(),
        der_order,
        acc_order,
        xmin,
        xmax,
        N,
        stencil_width = 3,
    )
    @test_throws ArgumentError D = periodic_derivative_operator(
        Holoborodko2008(),
        der_order,
        acc_order,
        xmin,
        xmax,
        N,
        stencil_width = 11,
    )

    acc_order = 4
    for stencil_width = 7:2:9
        D = periodic_derivative_operator(
            Holoborodko2008(),
            der_order,
            acc_order,
            xmin,
            xmax,
            N,
            stencil_width = stencil_width,
        )
        println(devnull, D)
        @test derivative_order(D) == der_order
        @test accuracy_order(D) == acc_order
        @test issymmetric(D) == true
        M = mass_matrix(D)
        @test M * Matrix(D) - Matrix(D)' * M ≈ zeros(T, N, N)
        mul!(res, D, x0)
        @test all(i -> abs(res[i]) < 50N * eps(T), stencil_width:length(res)-stencil_width)
        mul!(res, D, x1)
        @test all(
            i -> abs(res[i]) < 1000N * eps(T),
            stencil_width:length(res)-stencil_width,
        )
        mul!(res, D, x2)
        @test all(
            i -> isapprox(res[i], 2 * x0[i], atol = 5000N * eps(T)),
            stencil_width:length(res)-stencil_width,
        )
        mul!(res, D, x3)
        @test all(
            i -> isapprox(res[i], 6 * x1[i], atol = 5000N * eps(T)),
            stencil_width:length(res)-stencil_width,
        )
        mul!(res, D, x4)
        @test all(
            i -> isapprox(res[i], 12 * x2[i], atol = 7000N * eps(T)),
            stencil_width:length(res)-stencil_width,
        )
        mul!(res, D, x5)
        @test any(i -> !(res[i] ≈ 30 * x3[i]), stencil_width:length(res)-stencil_width)
    end
    @test_throws ArgumentError D = periodic_derivative_operator(
        Holoborodko2008(),
        der_order,
        acc_order,
        xmin,
        xmax,
        N,
        stencil_width = 5,
    )
    @test_throws ArgumentError D = periodic_derivative_operator(
        Holoborodko2008(),
        der_order,
        acc_order,
        xmin,
        xmax,
        N,
        stencil_width = 11,
    )

    acc_order = 6
    @test_throws ArgumentError D = periodic_derivative_operator(
        Holoborodko2008(),
        der_order,
        acc_order,
        xmin,
        xmax,
        N,
        stencil_width = 9,
    )
    @test_throws ArgumentError D = periodic_derivative_operator(
        Holoborodko2008(),
        der_order,
        acc_order,
        xmin,
        xmax,
        N,
        stencil_width = 13,
    )

    der_order = 3
    acc_order = 2
    @test_throws ArgumentError D = periodic_derivative_operator(
        Holoborodko2008(),
        der_order,
        acc_order,
        xmin,
        xmax,
        N,
        stencil_width = 5,
    )
    @test_throws ArgumentError D = periodic_derivative_operator(
        Holoborodko2008(),
        der_order,
        acc_order,
        xmin,
        xmax,
        N,
        stencil_width = 13,
    )
end

# Robust method LanczosLowNoise
@testset "LanczosLowNoise" begin
    @testset "T = $T" for T in (Float32, Float64)
        xmin = one(T)
        xmax = 2 * one(T)
        N = 100
        x1 = compute_coefficients(
            identity,
            periodic_derivative_operator(1, 2, xmin, xmax, N),
        )

        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        x5 = x2 .* x3
        x6 = x3 .* x3
        x7 = x4 .* x3

        res = fill(zero(eltype(x0)), length(x0))

        # first derivative operators
        der_order = 1
        acc_order = 2
        @testset "second-order operator, stencil width $stencil_width" for stencil_width =
                                                                           5:2:11
            D = periodic_derivative_operator(
                LanczosLowNoise();
                derivative_order = der_order,
                accuracy_order = acc_order,
                xmin = xmin,
                xmax = xmax,
                N = N,
                stencil_width = stencil_width,
            )
            println(devnull, D)
            @test derivative_order(D) == der_order
            @test accuracy_order(D) == acc_order
            @test issymmetric(D) == false
            M = mass_matrix(D)
            @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
            mul!(res, D, x0)
            @test all(
                i -> abs(res[i]) < 20 * eps(T),
                stencil_width:length(res)-stencil_width,
            )
            mul!(res, D, x1)
            @test all(i -> res[i] ≈ x0[i], stencil_width:length(res)-stencil_width)
            mul!(res, D, x2)
            @test all(i -> res[i] ≈ 2 * x1[i], stencil_width:length(res)-stencil_width)
        end
        @test_throws ArgumentError D = periodic_derivative_operator(
            LanczosLowNoise();
            derivative_order = der_order,
            accuracy_order = acc_order,
            xmin = xmin,
            xmax = xmax,
            N = N,
            stencil_width = 3,
        )
        @test_throws ArgumentError D = periodic_derivative_operator(
            LanczosLowNoise();
            derivative_order = der_order,
            accuracy_order = acc_order,
            xmin = xmin,
            xmax = xmax,
            N = N,
            stencil_width = 13,
        )

        acc_order = 4
        @testset "fourth-order operator, stencil width $stencil_width" for stencil_width =
                                                                           7:2:11
            D = periodic_derivative_operator(
                LanczosLowNoise();
                derivative_order = der_order,
                accuracy_order = acc_order,
                xmin = xmin,
                xmax = xmax,
                N = N,
                stencil_width = stencil_width,
            )
            println(devnull, D)
            @test derivative_order(D) == der_order
            @test accuracy_order(D) == acc_order
            @test issymmetric(D) == false
            M = mass_matrix(D)
            @test M * Matrix(D) + Matrix(D)' * M ≈ mass_matrix_boundary(D)
            mul!(res, D, x0)
            @test all(
                i -> abs(res[i]) < 50 * eps(T),
                stencil_width:length(res)-stencil_width,
            )
            mul!(res, D, x1)
            @test all(i -> res[i] ≈ x0[i], stencil_width:length(res)-stencil_width)
            mul!(res, D, x2)
            @test all(i -> res[i] ≈ 2 * x1[i], stencil_width:length(res)-stencil_width)
            mul!(res, D, x3)
            @test all(i -> res[i] ≈ 3 * x2[i], stencil_width:length(res)-stencil_width)
            mul!(res, D, x4)
            @test all(i -> res[i] ≈ 4 * x3[i], stencil_width:length(res)-stencil_width)
        end
        @test_throws ArgumentError D = periodic_derivative_operator(
            LanczosLowNoise();
            derivative_order = der_order,
            accuracy_order = acc_order,
            xmin = xmin,
            xmax = xmax,
            N = N,
            stencil_width = 5,
        )
        @test_throws ArgumentError D = periodic_derivative_operator(
            LanczosLowNoise();
            derivative_order = der_order,
            accuracy_order = acc_order,
            xmin = xmin,
            xmax = xmax,
            N = N,
            stencil_width = 13,
        )
    end

    @testset "T = $T" for T in (Rational{Int},)

        xmin = one(T)
        xmax = 2 * one(T)
        N = 100
        x1 = compute_coefficients(
            identity,
            periodic_derivative_operator(1, 2, xmin, xmax, N),
        )

        x0 = fill(one(eltype(x1)), length(x1))
        x2 = x1 .* x1
        x3 = x2 .* x1
        x4 = x2 .* x2
        x5 = x2 .* x3
        x6 = x3 .* x3
        x7 = x4 .* x3

        res = fill(zero(eltype(x0)), length(x0))

        # first derivative operators
        der_order = 1
        acc_order = 2
        @testset "second-order operator, stencil width $stencil_width" for stencil_width =
                                                                           5:2:11
            D = periodic_derivative_operator(
                LanczosLowNoise();
                derivative_order = der_order,
                accuracy_order = acc_order,
                xmin = xmin,
                xmax = xmax,
                N = N,
                stencil_width = stencil_width,
            )
            println(devnull, D)
            @test derivative_order(D) == der_order
            @test accuracy_order(D) == acc_order
            @test issymmetric(D) == false
            M = mass_matrix(D)
            @test M * Matrix(D) + Matrix(D)' * M == mass_matrix_boundary(D)
            mul!(res, D, x0)
            @test all(i -> iszero(res[i]), stencil_width:length(res)-stencil_width)
            mul!(res, D, x1)
            @test all(i -> res[i] == x0[i], stencil_width:length(res)-stencil_width)
            mul!(res, D, x2)
            @test all(i -> res[i] == 2 * x1[i], stencil_width:length(res)-stencil_width)
        end

        acc_order = 4
        @testset "fourth-order operator, stencil width $stencil_width" for stencil_width =
                                                                           7:2:11
            D = periodic_derivative_operator(
                LanczosLowNoise();
                derivative_order = der_order,
                accuracy_order = acc_order,
                xmin = xmin,
                xmax = xmax,
                N = N,
                stencil_width = stencil_width,
            )
            println(devnull, D)
            @test derivative_order(D) == der_order
            @test accuracy_order(D) == acc_order
            @test issymmetric(D) == false
            M = mass_matrix(D)
            @test M * Matrix(D) + Matrix(D)' * M == mass_matrix_boundary(D)
            mul!(res, D, x0)
            @test all(i -> iszero(res[i]), stencil_width:length(res)-stencil_width)
            mul!(res, D, x1)
            @test all(i -> res[i] == x0[i], stencil_width:length(res)-stencil_width)
            mul!(res, D, x2)
            @test all(i -> res[i] == 2 * x1[i], stencil_width:length(res)-stencil_width)
            mul!(res, D, x3)
            @test all(i -> res[i] == 3 * x2[i], stencil_width:length(res)-stencil_width)
            mul!(res, D, x4)
            @test all(i -> res[i] == 4 * x3[i], stencil_width:length(res)-stencil_width)
        end
    end
end # LanczosLowNoise

# Check rational operators
for T in (Float32, Float64)
    xmin = -one(T)
    xmax = one(T)
    for N in (8, 9), der_order in (1, 2), acc_order in (2, 3, 4)
        D = periodic_derivative_operator(der_order, acc_order, xmin, xmax, N)
        x = grid(D)
        u = @. sinpi(x) - cospi(x)^2 + exp(sinpi(x))
        println(devnull, D)

        @test Matrix(D^2) ≈ Matrix(D)^2

        @test Matrix(D^2) ≈ Matrix(D * D) ≈ Matrix((I * D) * D) ≈ Matrix(D * (D * I))
        @test Matrix(I - D^2) ≈ I - Matrix(D^2)
        @test Matrix(I + D) ≈ I + Matrix(D)
        @test Matrix(D - I) ≈ Matrix(D) - I

        rat = (I + 2D + 5 * D^2) * (2I * D - D^3 * 5I) * (D * 2 - D^2 * 5)
        @test rat.num_coef == (0.0, 0.0, 4.0, -2.0, -10.0, -45.0, 0.0, 125.0)
        @test rat.den_coef == (1.0,)
        println(devnull, rat)

        @test (I + one(T) / 2 * D) * u ≈ (u + D * u ./ 2)

        v = (I - D^2) * u
        @test inv(I - D^2) * v ≈ u

        v = (I - D^2) \ u
        @test D * v ≈ (D / (I - D^2)) * u

        rat = (I - D^2) / (I + D^4)
        println(devnull, rat)
        v = rat * u
        @test (I - D^2) \ (v + D^4 * v) ≈ u

        rat1 = (I - D^2) / (I + D^4)
        rat2 = (I + D^4) / (I - D^2)
        rat3 = (I - D^4) / (I - D^2)
        @test (rat2 + rat3) * u ≈ 2 * ((I - D^2) \ u)
        @test (rat1 * rat2) * u ≈ u

        @test integrate(u, D) ≈ sum(mass_matrix(D) * u)
        @test integrate(u -> u^2, u, D) ≈ dot(u, mass_matrix(D), u)
        @test integrate_boundary(u, D) ≈ 0.0
        @test integrate_boundary(u -> u^2, u, D) ≈ 0.0
    end

    for N in (8, 9), acc_order in (2, 3, 4)
        D1a = periodic_derivative_operator(1, acc_order, xmin, xmax, N)
        D2a = D1a^2
        D1b = fourier_derivative_operator(xmin, xmax, N)
        D2b = D1b^2
        D2c = periodic_derivative_operator(2, acc_order, xmin, xmax, N)

        @test Matrix(D1a // (I - D2a)) ≈ Matrix(D1a) / (I - Matrix(D2a))
        @test Matrix(D1b // (I - D2a)) ≈ Matrix(D1b) / (I - Matrix(D2a))
        @test Matrix(D1a // (I - D2b)) ≈ Matrix(D1a) / (I - Matrix(D2b))
        @test Matrix(D1b // (I - D2b)) ≈ Matrix(D1b) / (I - Matrix(D2b))
        @test Matrix(D1b // (I - D2b)) ≈ Matrix(D1b / (I - D2b))
        @test Matrix(D1a // (I - D2c)) ≈ Matrix(D1a) / (I - Matrix(D2c))
        @test Matrix(D1b // (I - D2c)) ≈ Matrix(D1b) / (I - Matrix(D2c))

        @test Matrix((I + D1a) // (I - D2a)) ≈ (I + Matrix(D1a)) / (I - Matrix(D2a))
        @test Matrix((I + D1b) // (I - D2a)) ≈ (I + Matrix(D1b)) / (I - Matrix(D2a))
        @test Matrix((I + D1a) // (I - D2b)) ≈ (I + Matrix(D1a)) / (I - Matrix(D2b))
        @test Matrix((I + D1b) // (I - D2b)) ≈ (I + Matrix(D1b)) / (I - Matrix(D2b))
        @test Matrix((I + D1b) // (I - D2b)) ≈ Matrix((I + D1b) / (I - D2b))
        @test Matrix((I + D1a) // (I - D2c)) ≈ (I + Matrix(D1a)) / (I - Matrix(D2c))
        @test Matrix((I + D1b) // (I - D2c)) ≈ (I + Matrix(D1b)) / (I - Matrix(D2c))
    end
end


# check construction of upwind operators
let N = 4
    @test_throws ArgumentError periodic_derivative_operator(
        1,
        2,
        0 // 1,
        (N - 1) // 1,
        N,
        1,
    )
    D = periodic_derivative_operator(1, 2, 0 // 1, N, N, 0)
    @test Matrix(D) ≈ [
        -3//2 2//1 -1//2 0//1
        0//1 -3//2 2//1 -1//2
        -1//2 0//1 -3//2 2//1
        2//1 -1//2 0//1 -3//2
    ]
    D = periodic_derivative_operator(1, 2, 0 // 1, N, N, -2)
    @test Matrix(D) ≈ [
        3//2 0//1 1//2 -2//1
        -2//1 3//2 0//1 1//2
        1//2 -2//1 3//2 0//1
        0//1 1//2 -2//1 3//2
    ]
end
