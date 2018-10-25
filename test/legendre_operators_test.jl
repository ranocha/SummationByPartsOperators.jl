using Test, SummationByPartsOperators
using LinearAlgebra
using SpecialFunctions


function accuracy_test!(res, ufunc, dufunc, D)
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    mul!(res, D, u)
    maximum(abs, du-res) < 5*length(res)*eps(eltype(res))
end

# Accuracy Tests
for T in (Float32, Float64)
    xmin = -one(T)
    xmax = one(T)

    for N in 2 .^ (1:4)
        D = legendre_derivative_operator(xmin, xmax, N)
        println(devnull, D)
        @test SummationByPartsOperators.derivative_order(D) == 1
        @test issymmetric(D) == false
        u = compute_coefficients(zero, D)
        res = D*u
        for k in 1:N-1
            ufunc = x -> x^k/typeof(x)(gamma(k+1))
            dufunc = x -> x^(k-1)/typeof(x)(gamma(k))
            @test accuracy_test!(res, ufunc, dufunc, D)
            xplot, duplot = evaluate_coefficients(res, D)
            @test maximum(abs, duplot - dufunc.(xplot)) < 5N*eps(T)
            @test abs(integrate(u, D)) < N*eps(T)
        end
    end
end

# Modal Filtering
for T in (Float32, Float64), filter_type in (ExponentialFilter(),)
    xmin = -one(T)
    xmax = one(T)

    for N in 2 .^ (1:4)
        D = legendre_derivative_operator(xmin, xmax, N)
        filter! = ConstantFilter(D, filter_type)
        u = compute_coefficients(zero, D)
        res = D*u
        for k in 1:N-1
            compute_coefficients!(u, x->exp(sinpi(x)), D)
            norm2_u = integrate(u->u^2, u, D)
            filter!(u)
            @test integrate(u->u^2, u, D) <= norm2_u
        end
    end
end
