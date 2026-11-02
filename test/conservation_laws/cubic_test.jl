using Test, SummationByPartsOperators
using OrdinaryDiffEqSSPRK, DiffEqCallbacks

@testset "Cubic" for T in (Float32, Float64), split_form in (Val{true}(), Val{false}())
    xmin = T(-1)
    xmax = T(1)
    N = 2^6
    u0func = sinpi
    tspan = (zero(T), one(T))

    # Fourier
    let D = fourier_derivative_operator(xmin, xmax, N)
        Di = dissipation_operator(TadmorWaagan2012Standard(), D)
        semi = CubicPeriodicSemidiscretization(D, Di, split_form)
        println(devnull, semi)
        ode = semidiscretize(u0func, semi, tspan)
        du = similar(ode.u0)
        semi(du, ode.u0, nothing, first(tspan))

        saving = SavingCallback(semi, saveat = range(tspan..., length = 10))
        sol = solve(ode, SSPRK104(), dt = T(1 / 5N), save_everystep = false,
                    callback = saving)
    end

    # Periodic FD
    for acc_order in 2:2:8
        D = periodic_derivative_operator(1, acc_order, xmin, xmax, N)
        Di = dissipation_operator(D)
        semi = CubicPeriodicSemidiscretization(D, Di, split_form)
        println(devnull, semi)
        ode = semidiscretize(u0func, semi, tspan)
        du = similar(ode.u0)
        semi(du, ode.u0, nothing, first(tspan))

        saving = SavingCallback(semi, saveat = range(tspan..., length = 10))
        sol = solve(ode, SSPRK104(), dt = T(1 / 5N), save_everystep = false,
                    callback = saving)
    end

    # Legendre
    let D = legendre_derivative_operator(xmin, xmax, N)
        Di = nothing
        semi = CubicNonperiodicSemidiscretization(D, Di, split_form, zero, zero)
        println(devnull, semi)
        ode = semidiscretize(u0func, semi, tspan)
        du = similar(ode.u0)
        semi(du, ode.u0, nothing, first(tspan))

        saving = SavingCallback(semi, saveat = range(tspan..., length = 10))
        sol = solve(ode, SSPRK104(), dt = T(1 / (1 + N^20)), save_everystep = false,
                    callback = saving)
    end

    # SBP FD
    for acc_order in 2:2:8
        D = derivative_operator(MattssonSvärdNordström2004(), 1, acc_order, xmin, xmax, N)
        Di = dissipation_operator(D)
        semi = CubicNonperiodicSemidiscretization(D, Di, split_form, zero, zero)
        println(devnull, semi)
        ode = semidiscretize(u0func, semi, tspan)
        du = similar(ode.u0)
        semi(du, ode.u0, nothing, first(tspan))

        saving = SavingCallback(semi, saveat = range(tspan..., length = 10))
        sol = solve(ode, SSPRK104(), dt = T(1 / 5N), save_everystep = false,
                    callback = saving)

        # boundary values of another type than `T`
        semi_int = CubicNonperiodicSemidiscretization(D, Di, split_form, t -> 1, t -> -1)
        semi_T = CubicNonperiodicSemidiscretization(D, Di, split_form,
                                                    t -> one(T), t -> -one(T))
        du_int = similar(ode.u0)
        semi_int(du_int, ode.u0, nothing, first(tspan))
        semi_T(du, ode.u0, nothing, first(tspan))
        @test du_int == du
    end
end

@testset "Godunov flux with arguments of different types" begin
    godunov_flux_cubic = SummationByPartsOperators.godunov_flux_cubic
    @test godunov_flux_cubic(2, 0.5) === godunov_flux_cubic(2.0, 0.5)
    @test godunov_flux_cubic(2, 0.5f0) === godunov_flux_cubic(2.0f0, 0.5f0)
end
