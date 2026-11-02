using Test, SummationByPartsOperators

@testset "Burgers" for T in (Float32, Float64), split_form in (Val{true}(), Val{false}())
    xmin = T(-1)
    xmax = T(1)
    N = 2^6
    u0func = sinpi
    tspan = (zero(T), one(T))

    @testset "Fourier" begin
        let D = fourier_derivative_operator(xmin, xmax, N)
            Di = dissipation_operator(TadmorWaagan2012Standard(), D)
            semi = BurgersPeriodicSemidiscretization(D, Di, split_form)
            println(devnull, semi)
            ode = semidiscretize(u0func, semi, tspan)
            du = similar(ode.u0)
            semi(du, ode.u0, nothing, first(tspan))
        end
    end

    @testset "Periodic FD" begin
        for acc_order in 2:2:8
            D = periodic_derivative_operator(1, acc_order, xmin, xmax, N)
            Di = dissipation_operator(D)
            semi = BurgersPeriodicSemidiscretization(D, Di, split_form)
            println(devnull, semi)
            ode = semidiscretize(u0func, semi, tspan)
            du = similar(ode.u0)
            semi(du, ode.u0, nothing, first(tspan))
        end
    end

    @testset "Legendre" begin
        let D = legendre_derivative_operator(xmin, xmax, N)
            Di = nothing
            semi = BurgersNonperiodicSemidiscretization(D, Di, split_form, zero, zero)
            println(devnull, semi)
            ode = semidiscretize(u0func, semi, tspan)
            du = similar(ode.u0)
            semi(du, ode.u0, nothing, first(tspan))
        end
    end

    @testset "SBP FD" begin
        for acc_order in 2:2:8
            D = derivative_operator(MattssonSvärdNordström2004(), 1, acc_order, xmin, xmax,
                                    N)
            Di = dissipation_operator(D)
            semi = BurgersNonperiodicSemidiscretization(D, Di, split_form, zero, zero)
            println(devnull, semi)
            ode = semidiscretize(u0func, semi, tspan)
            du = similar(ode.u0)
            semi(du, ode.u0, nothing, first(tspan))

            # boundary values of another type than `T`
            semi_int = BurgersNonperiodicSemidiscretization(D, Di, split_form,
                                                            t -> 1, t -> -1)
            semi_T = BurgersNonperiodicSemidiscretization(D, Di, split_form,
                                                          t -> one(T), t -> -one(T))
            du_int = similar(ode.u0)
            semi_int(du_int, ode.u0, nothing, first(tspan))
            semi_T(du, ode.u0, nothing, first(tspan))
            @test du_int == du
        end
    end
end

@testset "Godunov flux with arguments of different types" begin
    godunov_flux_burgers = SummationByPartsOperators.godunov_flux_burgers
    for (i, x) in ((1, 0.5), (-1, 0.5), (1, -0.5), (-1, -0.5))
        @test godunov_flux_burgers(i, x) === godunov_flux_burgers(float(i), x)
        @test godunov_flux_burgers(x, i) === godunov_flux_burgers(x, float(i))
        @test godunov_flux_burgers(i, Float32(x)) ===
              godunov_flux_burgers(Float32(i), Float32(x))
    end
end
