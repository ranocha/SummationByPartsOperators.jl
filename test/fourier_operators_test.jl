using Base.Test, SummationByPartsOperators


function accuracy_test!(res, ufunc, dufunc, D)
    u = compute_coefficients(ufunc, D)
    du = compute_coefficients(dufunc, D)
    A_mul_B!(res, D, u)
    maximum(abs, du-res) < 5*length(res)*eps(eltype(res))
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
            @test abs(integrate(u, D)) < N*eps(T)

            ufunc = x->cospi(k*x)
            dufunc = x->-typeof(x)(k*π)*sinpi(k*x)
            @test accuracy_test!(res, ufunc, dufunc, D)
            xplot, duplot = evaluate_coefficients(res, D)
            @test maximum(abs, duplot - dufunc.(xplot)) < 5N*eps(T)
            @test abs(integrate(u, D)) < N*eps(T)
        end
    end
end


# Spectral Viscosity.
source_SV = (Tadmor1989(), MadayTadmor1989(), TadmorWaagan2012Standard(), TadmorWaagan2012Convergent())
source_SSV = (Tadmor1993(),)

for T in (Float32, Float64), source in source_SV
    xmin = -one(T)
    xmax = one(T)

    for N in 2 .^ (3:6)
        D = fourier_derivative_operator(xmin, xmax, N)
        println(DevNull, D)
        @test issymmetric(D) == false

        Di = dissipation_operator(source, D)
        println(DevNull, Di)
        @test issymmetric(Di) == true
        Di_full = full(Di)
        @test maximum(abs, Di_full-Di_full') < 80*eps(T)

        @test maximum(eigvals(Symmetric(Di_full))) < 10N*eps(T)
    end
end

for T in (Float32, Float64), source in source_SSV
    xmin = -one(T)
    xmax = one(T)

    for N in 2 .^ (3:6), order in 1:3
        D = fourier_derivative_operator(xmin, xmax, N)
        println(DevNull, D)
        @test issymmetric(D) == false

        Di = dissipation_operator(source, D, order=order)
        println(DevNull, Di)
        @test issymmetric(Di) == true
        Di_full = full(Di)
        @test maximum(abs, Di_full-Di_full') < 80*eps(T)

        @test maximum(eigvals(Symmetric(Di_full))) < 15N*eps(T)
    end
end
