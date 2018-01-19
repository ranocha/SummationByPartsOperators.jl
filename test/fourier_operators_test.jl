using Base.Test, SummationByPartsOperators


function accuracy_test!(res, ufunc, dufunc, D)
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    A_mul_B!(res, D, u)
    maximum(abs, du-res) < 5N*eps(T)
end

# Accuracy tests.
for T in (Float32, Float64)
    xmin = -one(T)
    xmax = one(T)

    for N in 2 .^ (3:6)
        D = fourier_derivative_operator(xmin, xmax, N)
        println(DevNull, D)
        @test SummationByPartsOperators.derivative_order(D) == 1
        @test issymmetric(D) == false
        u = compute_coefficients(zero, D)
        res = D*u
        for k in 0:(N÷2)-1
            ufunc = x->sinpi(k*x)
            dufunc = x->typeof(x)(k*π)*cospi(k*x)
            @test accuracy_test!(res, ufunc, dufunc, D)
            xplot, duplot = evaluate_coefficients(res, D)
            @test maximum(abs, duplot - dufunc.(xplot)) < 5N*eps(T)

            ufunc = x->cospi(k*x)
            dufunc = x->-typeof(x)(k*π)*sinpi(k*x)
            @test accuracy_test!(res, ufunc, dufunc, D)
            xplot, duplot = evaluate_coefficients(res, D)
            @test maximum(abs, duplot - dufunc.(xplot)) < 5N*eps(T)
        end
    end
end
