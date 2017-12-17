using Base.Test
using SummationByPartsOperators

let T = Float32
    xmin = -one(T)
    xmax = 2*one(T)
    N = 101
    x1 = linspace(xmin, xmax, N)

    x0 = ones(x1)
    x2 = x1 .* x1
    x3 = x2 .* x1
    x4 = x2 .* x2
    x5 = x2 .* x3
    x6 = x3 .* x3
    x7 = x4 .* x3

    res = zeros(x0)

    # first derivative operators
    accuracy_order = 2
    D = periodic_central_derivative_operator(1, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == false
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < eps(T)
    A_mul_B!(res, D, x1); @test norm((res - x0)[accuracy_order:end-accuracy_order], Inf) < 20*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x1)[accuracy_order:end-accuracy_order], Inf) < 40*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 3 .* x2)[accuracy_order:end-accuracy_order], Inf) > 140*eps(T)

    accuracy_order = 4
    D = periodic_central_derivative_operator(1, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == false
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 10*eps(T)
    A_mul_B!(res, D, x1); @test norm((res - x0)[accuracy_order:end-accuracy_order], Inf) < 40*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x1)[accuracy_order:end-accuracy_order], Inf) < 60*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 3 .* x2)[accuracy_order:end-accuracy_order], Inf) < 140*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 4 .* x3)[accuracy_order:end-accuracy_order], Inf) < 300*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 5 .* x4)[accuracy_order:end-accuracy_order], Inf) > 500*eps(T)

    accuracy_order = 6
    D = periodic_central_derivative_operator(1, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == false
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 10*eps(T)
    A_mul_B!(res, D, x1); @test norm((res - x0)[accuracy_order:end-accuracy_order], Inf) < 40*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x1)[accuracy_order:end-accuracy_order], Inf) < 100*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 3 .* x2)[accuracy_order:end-accuracy_order], Inf) < 200*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 4 .* x3)[accuracy_order:end-accuracy_order], Inf) < 360*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 5 .* x4)[accuracy_order:end-accuracy_order], Inf) < 780*eps(T)
    A_mul_B!(res, D, x6); @test norm((res - 6 .* x5)[accuracy_order:end-accuracy_order], Inf) < 2200*eps(T)
    A_mul_B!(res, D, x7); @test norm((res - 7 .* x6)[accuracy_order:end-accuracy_order], Inf) > 2200*eps(T)

    # second derivative operators
    accuracy_order = 2
    D = periodic_central_derivative_operator(2, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == true
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < eps(T)
    A_mul_B!(res, D, x1); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 1400*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x0)[accuracy_order:end-accuracy_order], Inf) < 7000*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 6 .* x1)[accuracy_order:end-accuracy_order], Inf) > 7000*eps(T)

    accuracy_order = 4
    D = periodic_central_derivative_operator(2, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == true
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 400*eps(T)
    A_mul_B!(res, D, x1); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 6000*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x0)[accuracy_order:end-accuracy_order], Inf) < 11000*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 6 .* x1)[accuracy_order:end-accuracy_order], Inf) < 23000*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 12 .* x2)[accuracy_order:end-accuracy_order], Inf) < 40000*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 20 .* x3)[accuracy_order:end-accuracy_order], Inf) > 50000*eps(T)

    accuracy_order = 6
    D = periodic_central_derivative_operator(2, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == true
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 600*eps(T)
    A_mul_B!(res, D, x1); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 5000*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x0)[accuracy_order:end-accuracy_order], Inf) < 10020*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 6 .* x1)[accuracy_order:end-accuracy_order], Inf) < 30000*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 12 .* x2)[accuracy_order:end-accuracy_order], Inf) < 40000*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 20 .* x3)[accuracy_order:end-accuracy_order], Inf) < 100000*eps(T)
    A_mul_B!(res, D, x6); @test norm((res - 30 .* x4)[accuracy_order:end-accuracy_order], Inf) < 220000*eps(T)
    A_mul_B!(res, D, x7); @test norm((res - 42 .* x5)[accuracy_order:end-accuracy_order], Inf) < 220000*eps(T)

    tmp = D * x7
    @test typeof(tmp) == typeof(res)
    @test norm(tmp - res, Inf) < eps(T)

    # third derivative operators
    accuracy_order = 2
    D = periodic_central_derivative_operator(3, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == true # because this operator is zero!
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < eps(T)
    A_mul_B!(res, D, x1); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 1400*eps(T)
    A_mul_B!(res, D, x2); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 7000*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 6 .* x0)[accuracy_order:end-accuracy_order], Inf) > 7000*eps(T)

    accuracy_order = 4
    D = periodic_central_derivative_operator(3, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == false
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 400*eps(T)
    A_mul_B!(res, D, x1); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 20000*eps(T)
    A_mul_B!(res, D, x2); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 200000*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 6 .* x0)[accuracy_order:end-accuracy_order], Inf) < 300000*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 24 .* x1)[accuracy_order:end-accuracy_order], Inf) < 600000*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 60 .* x2)[accuracy_order:end-accuracy_order], Inf) > 600000*eps(T)

    accuracy_order = 6
    D = periodic_central_derivative_operator(3, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == false
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 300*eps(T)
    A_mul_B!(res, D, x1); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 80000*eps(T)
    A_mul_B!(res, D, x2); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 300000*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 6 .* x0)[accuracy_order:end-accuracy_order], Inf) < 700000*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 24 .* x1)[accuracy_order:end-accuracy_order], Inf) < 1000000*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 60 .* x2)[accuracy_order:end-accuracy_order], Inf) < 4000000*eps(T)
    A_mul_B!(res, D, x6); @test norm((res - 120 .* x3)[accuracy_order:end-accuracy_order], Inf) < 8000000*eps(T)
    A_mul_B!(res, D, x7); @test norm((res - 240 .* x4)[accuracy_order:end-accuracy_order], Inf) > 220000*eps(T)

    tmp = D * x7
    @test typeof(tmp) == typeof(res)
    @test norm(tmp - res, Inf) < eps(T)

    # higher derivative operators
    accuracy_order = 2
    @test_throws ArgumentError periodic_central_derivative_operator(4, accuracy_order, xmin, xmax, N)
    @test_throws ArgumentError periodic_central_derivative_operator(5, accuracy_order, xmin, xmax, N)
    @test_throws ArgumentError periodic_central_derivative_operator(6, accuracy_order, xmin, xmax, N)
end

let T = Float64
    xmin = -one(T)
    xmax = 2*one(T)
    N = 101
    x1 = linspace(xmin, xmax, N)

    x0 = ones(x1)
    x2 = x1 .* x1
    x3 = x2 .* x1
    x4 = x2 .* x2
    x5 = x2 .* x3
    x6 = x3 .* x3
    x7 = x4 .* x3
    x8 = x4 .* x4
    x9 = x4 .* x5

    res = zeros(x0)

    # first derivative operators
    accuracy_order = 2
    D = periodic_central_derivative_operator(1, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == false
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < eps(T)
    A_mul_B!(res, D, x1); @test norm((res - x0)[accuracy_order:end-accuracy_order], Inf) < 20*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x1)[accuracy_order:end-accuracy_order], Inf) < 70*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 3 .* x2)[accuracy_order:end-accuracy_order], Inf) > 140*eps(T)

    accuracy_order = 4
    D = periodic_central_derivative_operator(1, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == false
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 10*eps(T)
    A_mul_B!(res, D, x1); @test norm((res - x0)[accuracy_order:end-accuracy_order], Inf) < 40*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x1)[accuracy_order:end-accuracy_order], Inf) < 80*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 3 .* x2)[accuracy_order:end-accuracy_order], Inf) < 220*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 4 .* x3)[accuracy_order:end-accuracy_order], Inf) < 480*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 5 .* x4)[accuracy_order:end-accuracy_order], Inf) > 500*eps(T)

    accuracy_order = 6
    D = periodic_central_derivative_operator(1, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == false
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 10*eps(T)
    A_mul_B!(res, D, x1); @test norm((res - x0)[accuracy_order:end-accuracy_order], Inf) < 40*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x1)[accuracy_order:end-accuracy_order], Inf) < 120*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 3 .* x2)[accuracy_order:end-accuracy_order], Inf) < 300*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 4 .* x3)[accuracy_order:end-accuracy_order], Inf) < 420*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 5 .* x4)[accuracy_order:end-accuracy_order], Inf) < 1000*eps(T)
    A_mul_B!(res, D, x6); @test norm((res - 6 .* x5)[accuracy_order:end-accuracy_order], Inf) < 2200*eps(T)
    A_mul_B!(res, D, x7); @test norm((res - 7 .* x6)[accuracy_order:end-accuracy_order], Inf) > 2200*eps(T)

    tmp = D * x7
    @test typeof(tmp) == typeof(res)
    @test norm(tmp - res, Inf) < eps(T)

    # second derivative operators
    accuracy_order = 2
    D = periodic_central_derivative_operator(2, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == true
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < eps(T)
    A_mul_B!(res, D, x1); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 1400*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x0)[accuracy_order:end-accuracy_order], Inf) < 7000*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 6 .* x1)[accuracy_order:end-accuracy_order], Inf) > 7000*eps(T)

    accuracy_order = 4
    D = periodic_central_derivative_operator(2, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == true
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 400*eps(T)
    A_mul_B!(res, D, x1); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 5000*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x0)[accuracy_order:end-accuracy_order], Inf) < 11000*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 6 .* x1)[accuracy_order:end-accuracy_order], Inf) < 23000*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 12 .* x2)[accuracy_order:end-accuracy_order], Inf) < 40000*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 20 .* x3)[accuracy_order:end-accuracy_order], Inf) > 50000*eps(T)

    accuracy_order = 6
    D = periodic_central_derivative_operator(2, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == true
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 300*eps(T)
    A_mul_B!(res, D, x1); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 5000*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x0)[accuracy_order:end-accuracy_order], Inf) < 10020*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 6 .* x1)[accuracy_order:end-accuracy_order], Inf) < 30000*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 12 .* x2)[accuracy_order:end-accuracy_order], Inf) < 28000*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 20 .* x3)[accuracy_order:end-accuracy_order], Inf) < 100000*eps(T)
    A_mul_B!(res, D, x6); @test norm((res - 30 .* x4)[accuracy_order:end-accuracy_order], Inf) < 220000*eps(T)
    A_mul_B!(res, D, x7); @test norm((res - 42 .* x5)[accuracy_order:end-accuracy_order], Inf) > 220000*eps(T)

    tmp = D * x7
    @test typeof(tmp) == typeof(res)
    @test norm(tmp - res, Inf) < eps(T)

    # third derivative operators
    accuracy_order = 2
    D = periodic_central_derivative_operator(3, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == true # because this operator is zero!
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < eps(T)
    A_mul_B!(res, D, x1); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 1400*eps(T)
    A_mul_B!(res, D, x2); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 7000*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 6 .* x0)[accuracy_order:end-accuracy_order], Inf) > 7000*eps(T)

    accuracy_order = 4
    D = periodic_central_derivative_operator(3, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == false
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 400*eps(T)
    A_mul_B!(res, D, x1); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 20000*eps(T)
    A_mul_B!(res, D, x2); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 200000*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 6 .* x0)[accuracy_order:end-accuracy_order], Inf) < 300000*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 24 .* x1)[accuracy_order:end-accuracy_order], Inf) < 600000*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 60 .* x2)[accuracy_order:end-accuracy_order], Inf) > 600000*eps(T)

    accuracy_order = 6
    D = periodic_central_derivative_operator(3, accuracy_order, xmin, xmax, N)
    @test issymmetric(D) == false
    A_mul_B!(res, D, x0); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 300*eps(T)
    A_mul_B!(res, D, x1); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 70000*eps(T)
    A_mul_B!(res, D, x2); @test norm(res[accuracy_order:end-accuracy_order], Inf) < 200000*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 6 .* x0)[accuracy_order:end-accuracy_order], Inf) < 600000*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 24 .* x1)[accuracy_order:end-accuracy_order], Inf) < 1000000*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 60 .* x2)[accuracy_order:end-accuracy_order], Inf) < 4000000*eps(T)
    A_mul_B!(res, D, x6); @test norm((res - 120 .* x3)[accuracy_order:end-accuracy_order], Inf) < 8000000*eps(T)
    A_mul_B!(res, D, x7); @test norm((res - 240 .* x4)[accuracy_order:end-accuracy_order], Inf) > 220000*eps(T)

    tmp = D * x7
    @test typeof(tmp) == typeof(res)
    @test norm(tmp - res, Inf) < eps(T)

    # higher derivative operators
    accuracy_order = 2
    @test_throws ArgumentError periodic_central_derivative_operator(4, accuracy_order, xmin, xmax, N)
    @test_throws ArgumentError periodic_central_derivative_operator(5, accuracy_order, xmin, xmax, N)
    @test_throws ArgumentError periodic_central_derivative_operator(6, accuracy_order, xmin, xmax, N)
end


for T in (Float32, Float64), accuracy_order in 2:2:12, derivative_order in 1:3
    xmin = zero(T)
    xmax = one(T)
    N = 21
    D_central = periodic_central_derivative_operator(derivative_order, accuracy_order, xmin, xmax, N)
    D_general = periodic_derivative_operator(derivative_order, accuracy_order, xmin, xmax, N)
    @test norm(D_central.coefficients.lower_coef - D_general.coefficients.lower_coef) < 3*eps(T)
    @test abs(D_central.coefficients.central_coef - D_general.coefficients.central_coef) < 3*eps(T)
    @test norm(D_central.coefficients.upper_coef - D_general.coefficients.upper_coef) < 3*eps(T)
end


for T in (Float32, Float64), accuracy_order in 1:10, derivative_order in 1:3
    xmin = zero(T)
    xmax = 5*one(T)
    N = 51
    x = linspace(xmin, xmax, N)
    u = x.^5
    dest_serial = zeros(u)
    dest_threads= zeros(u)
    dest_full   = zeros(u)
    dest_sparse = zeros(u)
    D_serial = periodic_derivative_operator(derivative_order, accuracy_order, xmin, xmax, N)
    D_threads= periodic_derivative_operator(derivative_order, accuracy_order, xmin, xmax, N, Val{:threads}())
    D_full = full(D_serial)
    D_sparse = sparse(D_serial)
    A_mul_B!(dest_serial, D_serial, u)
    A_mul_B!(dest_threads, D_threads, u)
    A_mul_B!(dest_full, D_full, u)
    A_mul_B!(dest_sparse, D_sparse, u)
    for i in eachindex(u)
        @test dest_serial[i] ≈ dest_threads[i]
        @test isapprox(dest_serial[i], dest_full[i], rtol=5*sqrt(eps(T)))
        @test isapprox(dest_serial[i], dest_sparse[i], rtol=5*sqrt(eps(T)))
    end
end
