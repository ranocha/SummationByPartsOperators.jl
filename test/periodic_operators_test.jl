using Base.Test
using SummationByPartsOperators

let T = Float32
    xmin = -one(T)
    xmax = 2*one(T)
    N = 101
    x1 = linspace(xmin, xmax, N)
    Δx = (xmax - xmin) / (N-1)

    x0 = ones(x1)
    x2 = x1 .* x1
    x3 = x2 .* x1
    x4 = x2 .* x2
    x5 = x2 .* x3
    x6 = x3 .* x3
    x7 = x4 .* x3

    res = similar(x0)

    # first derivative operators
    accurcy_order = 2
    D = periodic_central_derivative_operator(Δx, 1, accurcy_order)
    A_mul_B!(res, D, x0); @test norm(res[accurcy_order:end-accurcy_order], Inf) < eps(T)
    A_mul_B!(res, D, x1); @test norm((res - x0)[accurcy_order:end-accurcy_order], Inf) < 20*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x1)[accurcy_order:end-accurcy_order], Inf) < 40*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 3 .* x2)[accurcy_order:end-accurcy_order], Inf) > 140*eps(T)

    accurcy_order = 4
    D = periodic_central_derivative_operator(Δx, 1, accurcy_order)
    A_mul_B!(res, D, x0); @test norm(res[accurcy_order:end-accurcy_order], Inf) < 10*eps(T)
    A_mul_B!(res, D, x1); @test norm((res - x0)[accurcy_order:end-accurcy_order], Inf) < 40*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x1)[accurcy_order:end-accurcy_order], Inf) < 60*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 3 .* x2)[accurcy_order:end-accurcy_order], Inf) < 140*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 4 .* x3)[accurcy_order:end-accurcy_order], Inf) < 300*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 5 .* x4)[accurcy_order:end-accurcy_order], Inf) > 500*eps(T)

    accurcy_order = 6
    D = periodic_central_derivative_operator(Δx, 1, accurcy_order)
    A_mul_B!(res, D, x0); @test norm(res[accurcy_order:end-accurcy_order], Inf) < 10*eps(T)
    A_mul_B!(res, D, x1); @test norm((res - x0)[accurcy_order:end-accurcy_order], Inf) < 40*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x1)[accurcy_order:end-accurcy_order], Inf) < 100*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 3 .* x2)[accurcy_order:end-accurcy_order], Inf) < 200*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 4 .* x3)[accurcy_order:end-accurcy_order], Inf) < 360*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 5 .* x4)[accurcy_order:end-accurcy_order], Inf) < 780*eps(T)
    A_mul_B!(res, D, x6); @test norm((res - 6 .* x5)[accurcy_order:end-accurcy_order], Inf) < 2200*eps(T)
    A_mul_B!(res, D, x7); @test norm((res - 7 .* x6)[accurcy_order:end-accurcy_order], Inf) > 2200*eps(T)

    tmp = D * x7
    @test typeof(tmp) == typeof(res)
    @test norm(tmp - res, Inf) < eps(T)
end

let T = Float64
    xmin = -one(T)
    xmax = 2*one(T)
    N = 101
    x1 = linspace(xmin, xmax, N)
    Δx = (xmax - xmin) / (N-1)

    x0 = ones(x1)
    x2 = x1 .* x1
    x3 = x2 .* x1
    x4 = x2 .* x2
    x5 = x2 .* x3
    x6 = x3 .* x3
    x7 = x4 .* x3
    x8 = x4 .* x4
    x9 = x4 .* x5

    res = similar(x0)

    # first derivative operators
    accurcy_order = 2
    D = periodic_central_derivative_operator(Δx, 1, accurcy_order)
    A_mul_B!(res, D, x0); @test norm(res[accurcy_order:end-accurcy_order], Inf) < eps(T)
    A_mul_B!(res, D, x1); @test norm((res - x0)[accurcy_order:end-accurcy_order], Inf) < 20*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x1)[accurcy_order:end-accurcy_order], Inf) < 70*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 3 .* x2)[accurcy_order:end-accurcy_order], Inf) > 140*eps(T)

    accurcy_order = 4
    D = periodic_central_derivative_operator(Δx, 1, accurcy_order)
    A_mul_B!(res, D, x0); @test norm(res[accurcy_order:end-accurcy_order], Inf) < 10*eps(T)
    A_mul_B!(res, D, x1); @test norm((res - x0)[accurcy_order:end-accurcy_order], Inf) < 40*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x1)[accurcy_order:end-accurcy_order], Inf) < 80*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 3 .* x2)[accurcy_order:end-accurcy_order], Inf) < 220*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 4 .* x3)[accurcy_order:end-accurcy_order], Inf) < 480*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 5 .* x4)[accurcy_order:end-accurcy_order], Inf) > 500*eps(T)

    accurcy_order = 6
    D = periodic_central_derivative_operator(Δx, 1, accurcy_order)
    A_mul_B!(res, D, x0); @test norm(res[accurcy_order:end-accurcy_order], Inf) < 10*eps(T)
    A_mul_B!(res, D, x1); @test norm((res - x0)[accurcy_order:end-accurcy_order], Inf) < 40*eps(T)
    A_mul_B!(res, D, x2); @test norm((res - 2 .* x1)[accurcy_order:end-accurcy_order], Inf) < 120*eps(T)
    A_mul_B!(res, D, x3); @test norm((res - 3 .* x2)[accurcy_order:end-accurcy_order], Inf) < 300*eps(T)
    A_mul_B!(res, D, x4); @test norm((res - 4 .* x3)[accurcy_order:end-accurcy_order], Inf) < 420*eps(T)
    A_mul_B!(res, D, x5); @test norm((res - 5 .* x4)[accurcy_order:end-accurcy_order], Inf) < 1000*eps(T)
    A_mul_B!(res, D, x6); @test norm((res - 6 .* x5)[accurcy_order:end-accurcy_order], Inf) < 2200*eps(T)
    A_mul_B!(res, D, x7); @test norm((res - 7 .* x6)[accurcy_order:end-accurcy_order], Inf) > 2200*eps(T)

    tmp = D * x7
    @test typeof(tmp) == typeof(res)
    @test norm(tmp - res, Inf) < eps(T)
end
