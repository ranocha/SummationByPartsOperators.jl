using Test, SummationByPartsOperators

for T in (Float32, Float64), split_form in (Val{true}(), Val{false}())
    xmin = T(-1)
    xmax = T(1)
    N = 2^6
    u0func = sinpi
    tspan = (zero(T), one(T))

    # Fourier
    let D = fourier_derivative_operator(xmin, xmax, N)
        Di = dissipation_operator(TadmorWaagan2012Standard(), D)
        semi = BurgersPeriodicSemidiscretization(D, Di, split_form)
        println(devnull, semi)
        ode = semidiscretize(u0func, semi, tspan)
        du = similar(ode.u0)
        semi(du, ode.u0, nothing, first(tspan))
    end

    # Periodic FD
    for acc_order in 2:2:8
        D = periodic_derivative_operator(1, acc_order, xmin, xmax, N)
        Di = dissipation_operator(D)
        semi = BurgersPeriodicSemidiscretization(D, Di, split_form)
        println(devnull, semi)
        ode = semidiscretize(u0func, semi, tspan)
        du = similar(ode.u0)
        semi(du, ode.u0, nothing, first(tspan))
    end

    # Legendre
    let D = legendre_derivative_operator(xmin, xmax, N)
        Di = nothing
        semi = BurgersNonperiodicSemidiscretization(D, Di, split_form, zero, zero)
        println(devnull, semi)
        ode = semidiscretize(u0func, semi, tspan)
        du = similar(ode.u0)
        semi(du, ode.u0, nothing, first(tspan))
    end

    # SBP FD
    for acc_order in 2:2:8
        D = derivative_operator(MattssonSvärdNordström2004(), 1, acc_order, xmin, xmax, N)
        Di = dissipation_operator(D)
        semi = BurgersNonperiodicSemidiscretization(D, Di, split_form, zero, zero)
        println(devnull, semi)
        ode = semidiscretize(u0func, semi, tspan)
        du = similar(ode.u0)
        semi(du, ode.u0, nothing, first(tspan))
    end
end
