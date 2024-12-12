using Test, SummationByPartsOperators
using LinearAlgebra


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

    for N in 2 .^ (3:6)
        D = fourier_derivative_operator(xmin, xmax, N)
        println(devnull, D)
        @test SummationByPartsOperators.derivative_order(D) == 1
        @test issymmetric(D) == false
        @test SummationByPartsOperators.xmin(D) ≈ xmin
        @test SummationByPartsOperators.xmax(D) ≈ xmax
        M = mass_matrix(D)
        @test isapprox(M * Matrix(D) + Matrix(D)' * M, zeros(T, N, N), atol = N*eps(T))
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

        # mass matrix scaling
        x1 = grid(D)
        M = @inferred mass_matrix(D)
        u = sinpi.(x1)
        v = copy(u)
        scale_by_mass_matrix!(v, D)
        @test_throws DimensionMismatch scale_by_mass_matrix!(@view(v[(begin + 1):(end - 1)]), D)
        @test v ≈ M * u
        scale_by_inverse_mass_matrix!(v, D)
        @test_throws DimensionMismatch scale_by_inverse_mass_matrix!(@view(v[(begin + 1):(end - 1)]), D)
        @test v ≈ u
    end
end


# Check Fourier polynomial/rational operators
for T in (Float32, Float64)
    xmin = -one(T)
    xmax = one(T)
    for N in (8, 9)
        D = @inferred fourier_derivative_operator(xmin, xmax, N)
        x = @inferred grid(D)
        u = @. sinpi(x) - cospi(x)^2 + exp(sinpi(x))
        println(devnull, D)

        # see e.g. Steven G. Johnson (2011) Notes on FFT based differentiation
        @test (Matrix(D^2) ≈ Matrix(D)^2) == isodd(N)

        @test Matrix(D^2) ≈ Matrix(D * D) ≈ Matrix((I * D) * D) ≈ Matrix(D * (D * I))
        @test issymmetric(I - D^2)
        @test !issymmetric(I + D)
        @test !issymmetric(D - I)

        @test @inferred(SummationByPartsOperators.xmin(D^2)) ≈ xmin
        @test @inferred(SummationByPartsOperators.xmax(D^2)) ≈ xmax

        poly = @inferred (I + 2D + 5*D^2) * (2I * D - D^3 * 5I) * (D*2 - D^2 * 5)
        @test poly.coef == (0.0, 0.0, 4.0, -2.0, -10.0, -45.0, 0.0, 125.0)
        println(devnull, poly)

        @test @inferred(I + one(T)/2*D) * u ≈ (u + D*u ./ 2)

        v = (I - D^2) * u
        @test inv(I - D^2) * v ≈ u

        @test @inferred(SummationByPartsOperators.xmin(inv(I - D^2))) ≈ xmin
        @test @inferred(SummationByPartsOperators.xmax(inv(I - D^2))) ≈ xmax

        v = @inferred(I - D^2) \ u
        @test D * v ≈ (D / (I - D^2)) * u

        rat = @inferred((I - D^2) / (I + D^4))
        println(devnull, rat)
        v = rat * u
        @test (I - D^2) \ (v + D^4 * v) ≈ u

        rat1 = @inferred((I - D^2) / (I + D^4))
        rat2 = @inferred((I + D^4) / (I - D^2))
        rat3 = @inferred((I - D^4) / (I - D^2))
        @test @inferred(rat2 + rat3) * u ≈ 2 * ((I - D^2) \ u)
        @test @inferred(rat1 * rat2) * u ≈ u

        @test integrate(u, D) ≈ sum(mass_matrix(D) * u)
        @test integrate(u->u^2, u, D) ≈ dot(u, mass_matrix(D), u)

        # combine rational operators and scalars
        @test @inferred(2 * rat1) * u ≈ rat1 * (2 * u)
        @test @inferred(rat1 * 2) * u ≈ rat1 * (2 * u)
        @test @inferred(2 / rat1) * u ≈ inv(rat1) * (2 * u)
        @test @inferred(rat1 / 2) * u ≈ rat1 * (u / 2)
        @test @inferred(2 \ rat1) * u ≈ rat1 * (u / 2)
        @test @inferred(rat1 \ 2) * u ≈ inv(rat1) * (2 * u)

        # combine rational operators and uniform scaling via *, /, \
        @test @inferred((2 * I) * rat1) == 2 * rat1
        @test @inferred(rat1 * (2 * I)) == 2 * rat1
        @test @inferred((2 * I) / rat1) == 2 / rat1
        @test @inferred(rat1 / (2 * I)) == rat1 / 2
        @test @inferred((2 * I) \ rat1) == rat1 / 2
        @test @inferred(rat1 \ (2 * I)) == 2 / rat1

        # combine rational operators and uniform scaling via +, -
        @test @inferred(I + rat1) == D^0 + rat1
        @test @inferred(rat1 + I) == D^0 + rat1
        @test @inferred(I - rat1) == D^0 - rat1
        @test @inferred(rat1 - I) == rat1 - D^0
        @test @inferred(-rat1) == @inferred(-1 * rat1)
    end
end


# (Super) Spectral Viscosity
source_SV = (Tadmor1989(), MadayTadmor1989(), TadmorWaagan2012Standard(), TadmorWaagan2012Convergent())
source_SSV = (Tadmor1993(),)

for T in (Float32, Float64), source in source_SV
    xmin = -one(T)
    xmax = one(T)

    for N in 2 .^ (3:6)
        D = fourier_derivative_operator(xmin, xmax, N)
        println(devnull, D)
        @test issymmetric(D) == false

        Di = dissipation_operator(source, D)
        println(devnull, Di)
        @test issymmetric(Di) == true
        Di_full = Matrix(Di)
        @test maximum(abs, Di_full-Di_full') < 80*eps(T)

        @test maximum(eigvals(Symmetric(Di_full))) < 10N*eps(T)
    end
end

for T in (Float32, Float64), source in source_SSV
    xmin = -one(T)
    xmax = one(T)

    for N in 2 .^ (3:6), order in 1:3
        D = fourier_derivative_operator(xmin, xmax, N)
        println(devnull, D)
        @test issymmetric(D) == false

        Di = dissipation_operator(source, D, order=order)
        println(devnull, Di)
        @test issymmetric(Di) == true
        Di_full = Matrix(Di)
        @test maximum(abs, Di_full-Di_full') < 80*eps(T)

        @test maximum(eigvals(Symmetric(Di_full))) < 15N*eps(T)
    end
end


# Modal Filtering
for T in (Float32, Float64), filter_type in (ExponentialFilter(),)
    xmin = -one(T)
    xmax = one(T)

    for N in 2 .^ (1:4)
        D = fourier_derivative_operator(xmin, xmax, N)
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
