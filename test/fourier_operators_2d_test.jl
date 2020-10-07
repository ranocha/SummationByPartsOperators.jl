using Test, SummationByPartsOperators
using LinearAlgebra


function accuracy_test!(res, ufunc, dufunc, D, direction)
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    mul!(res, D, u, direction)
    maximum(abs, du-res) < 5*length(res)*eps(eltype(res))
end

# Accuracy Tests
for T in (Float32, Float64)
    xmin = -one(T)
    xmax = one(T)
    ymin = -one(T)
    ymax = one(T)

    for N in 2 .^ (3:6)
        D = fourier_derivative_operator(xmin, xmax, N, ymin, ymax, N)
        println(devnull, D)
        @test SummationByPartsOperators.derivative_order(D) == 1
        @test issymmetric(D) == false
        u = compute_coefficients(zero ∘ eltype, D)
        res = similar(u)
        for k in 0:(N÷2)-1
            ufunc = x->sinpi(k*x[1])
            dufunc = x->eltype(x)(k*π)*cospi(k*x[1])
            @test accuracy_test!(res, ufunc, dufunc, D, Val(:x))
            xplot, duplot = evaluate_coefficients(res, D)
            @test maximum(abs, duplot - dufunc.(xplot)) < 5N*eps(T)
            @test abs(integrate(u, D)) < N*eps(T)

            ufunc = x->cospi(k*x[1])
            dufunc = x->-eltype(x)(k*π)*sinpi(k*x[1])
            @test accuracy_test!(res, ufunc, dufunc, D, Val(:x))
            xplot, duplot = evaluate_coefficients(res, D)
            @test maximum(abs, duplot - dufunc.(xplot)) < 5N*eps(T)
            @test abs(integrate(u, D)) < N*eps(T)

            ufunc = x->sinpi(k*x[2])
            dufunc = x->eltype(x)(k*π)*cospi(k*x[2])
            @test accuracy_test!(res, ufunc, dufunc, D, Val(:y))
            xplot, duplot = evaluate_coefficients(res, D)
            @test maximum(abs, duplot - dufunc.(xplot)) < 5N*eps(T)
            @test abs(integrate(u, D)) < N*eps(T)

            ufunc = x->cospi(k*x[2])
            dufunc = x->-eltype(x)(k*π)*sinpi(k*x[2])
            @test accuracy_test!(res, ufunc, dufunc, D, Val(:y))
            xplot, duplot = evaluate_coefficients(res, D)
            @test maximum(abs, duplot - dufunc.(xplot)) < 5N*eps(T)
            @test abs(integrate(u, D)) < N*eps(T)
        end
    end
end

