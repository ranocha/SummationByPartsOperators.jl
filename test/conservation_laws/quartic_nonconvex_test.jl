using Test, SummationByPartsOperators
using OrdinaryDiffEq, DiffEqCallbacks

for T in (Float32, Float64), split_form in (Val{true}(), Val{false}())
    xmin = T(-1)
    xmax = T(1)
    N = 2^6
    u0func = sinpi
    tspan = (zero(T), one(T)/50)

    # Fourier
    let D = fourier_derivative_operator(xmin, xmax, N)
        Di = dissipation_operator(TadmorWaagan2012Standard(), D)
        semidisc = QuarticNonconvexPeriodicSemidiscretisation(D, Di, split_form)
        println(devnull, semidisc)
        ode = semidiscretise(u0func, semidisc, tspan)
        du = similar(ode.u0)
        semidisc(du, ode.u0, nothing, first(tspan))

        saving = SavingCallback(semidisc, saveat=range(tspan..., length=10))
        sol = solve(ode, SSPRK104(), dt=1/20N, save_everystep=false, callback=saving)
    end

    # Periodic FD
    for acc_order in 2:2:8
        D = periodic_derivative_operator(1, acc_order, xmin, xmax, N)
        Di = dissipation_operator(D)
        semidisc = QuarticNonconvexPeriodicSemidiscretisation(D, Di, split_form)
        println(devnull, semidisc)
        ode = semidiscretise(u0func, semidisc, tspan)
        du = similar(ode.u0)
        semidisc(du, ode.u0, nothing, first(tspan))

        saving = SavingCallback(semidisc, saveat=range(tspan..., length=10))
        sol = solve(ode, SSPRK104(), dt=1/20N, save_everystep=false, callback=saving)
    end
end
